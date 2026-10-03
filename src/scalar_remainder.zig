const common = @import("scalar_c23_common.zig");
const types = @import("scalar_c23_types.zig");
// Integer significands keep quotient selection exact across all supported exponent gaps.

const Reduction = enum {
    truncating,
    nearest_even,
};

const Format = struct {
    Float: type,
    Bits: type,
    Ieee: type,
    fraction_bits: u32,
    precision: u32,
    exponent_bias: i32,
    sign_mask: u64,
    exponent_mask: u64,
    fraction_mask: u64,
};

const h_format: Format = .{
    .Float = f16,
    .Bits = u16,
    .Ieee = common.H,
    .fraction_bits = 10,
    .precision = 11,
    .exponent_bias = 15,
    .sign_mask = 0x8000,
    .exponent_mask = 0x7c00,
    .fraction_mask = 0x03ff,
};

const f_format: Format = .{
    .Float = f32,
    .Bits = u32,
    .Ieee = common.F,
    .fraction_bits = 23,
    .precision = 24,
    .exponent_bias = 127,
    .sign_mask = 0x8000_0000,
    .exponent_mask = 0x7f80_0000,
    .fraction_mask = 0x007f_ffff,
};

const d_format: Format = .{
    .Float = f64,
    .Bits = u64,
    .Ieee = common.D,
    .fraction_bits = 52,
    .precision = 53,
    .exponent_bias = 1023,
    .sign_mask = 0x8000_0000_0000_0000,
    .exponent_mask = 0x7ff0_0000_0000_0000,
    .fraction_mask = 0x000f_ffff_ffff_ffff,
};

const Parts = struct {
    negative: bool,
    mantissa: u64,
    exponent: i32,
};

const Division = struct {
    remainder: u64,
    quotient_low: u8,
};

fn bitLength(value: u64) u32 {
    if (value == 0) return 0;
    return 64 - @as(u32, @intCast(@clz(value)));
}

fn decompose(comptime format: Format, value: format.Float) Parts {
    const bits = format.Ieee.toBits(value);
    const exponent_mask: format.Bits = @intCast(format.exponent_mask);
    const fraction_mask: format.Bits = @intCast(format.fraction_mask);
    const sign_mask: format.Bits = @intCast(format.sign_mask);
    const negative = bits & sign_mask != 0;
    const raw_exponent: u32 = @intCast((bits & exponent_mask) >> format.fraction_bits);
    const fraction: u64 = @intCast(bits & fraction_mask);

    if (raw_exponent == 0) {
        return .{
            .negative = negative,
            .mantissa = fraction,
            .exponent = 1 - format.exponent_bias - @as(i32, @intCast(format.fraction_bits)),
        };
    }

    return .{
        .negative = negative,
        .mantissa = fraction | (@as(u64, 1) << @intCast(format.fraction_bits)),
        .exponent = @as(i32, @intCast(raw_exponent)) - format.exponent_bias -
            @as(i32, @intCast(format.fraction_bits)),
    };
}

fn divideShifted(numerator: u64, denominator: u64, shift: u32) Division {
    stdAssert(denominator != 0);
    // Every divisor here is a significand of at most 53 bits.
    const chunk_max: u32 = @clz(denominator);
    stdAssert(chunk_max >= 11);

    // A power-of-two divisor has no significand bits below its width, so the
    // chunk loop can only reproduce a shift and a mask: the remainder is the
    // numerator truncated to the surviving width, and the quotient is that
    // numerator shifted by the same amount.
    if (denominator & (denominator - 1) == 0) {
        const width: u32 = bitLength(denominator) - 1;
        if (shift >= width) {
            // The quotient is the numerator shifted left by `shift - width`,
            // so its low three bits come from the numerator's low three.
            if (shift - width >= 3) return .{ .remainder = 0, .quotient_low = 0 };
            const kept: u8 = @intCast(numerator & 7);
            return .{
                .remainder = 0,
                .quotient_low = (kept << @as(u3, @intCast(shift - width))) & 7,
            };
        }
        const kept: u32 = width - shift;
        const mask = (@as(u64, 1) << @as(u6, @intCast(kept))) - 1;
        return .{
            .remainder = (numerator & mask) << @as(u6, @intCast(shift)),
            .quotient_low = @intCast((numerator >> @as(u6, @intCast(kept))) & 7),
        };
    }

    var remainder = numerator % denominator;
    var quotient_low: u8 = @intCast((numerator / denominator) & 7);
    var remaining_shift = shift;
    while (remaining_shift != 0) {
        if (remainder == 0) {
            if (remaining_shift >= 3) return .{ .remainder = 0, .quotient_low = 0 };
            quotient_low = (quotient_low << @as(u3, @intCast(remaining_shift))) & 7;
            return .{ .remainder = 0, .quotient_low = quotient_low };
        }

        const chunk: u32 = @min(remaining_shift, chunk_max);
        // remainder < denominator, so (remainder << chunk) fits in u64.
        const expanded = remainder << @as(u6, @intCast(chunk));
        const chunk_quotient = expanded / denominator;
        remainder = expanded % denominator;
        if (chunk >= 3) {
            quotient_low = @intCast(chunk_quotient & 7);
        } else {
            quotient_low = ((quotient_low << @as(u3, @intCast(chunk))) |
                @as(u8, @intCast(chunk_quotient))) & 7;
        }
        remaining_shift -= chunk;
    }

    return .{ .remainder = remainder, .quotient_low = quotient_low };
}

fn stdAssert(condition: bool) void {
    if (!condition) @panic("scalar remainder invariant failed");
}

fn packExact(
    comptime format: Format,
    magnitude: u64,
    exponent: i32,
    negative: bool,
) format.Float {
    const sign_mask: format.Bits = @intCast(format.sign_mask);
    const fraction_mask: format.Bits = @intCast(format.fraction_mask);
    const exponent_mask: format.Bits = @intCast(format.exponent_mask);
    const sign_bits: format.Bits = if (negative) sign_mask else 0;
    if (magnitude == 0) return format.Ieee.fromBits(sign_bits);

    const magnitude_bits = bitLength(magnitude);
    stdAssert(magnitude_bits <= format.precision);

    const top_exponent = exponent + @as(i32, @intCast(magnitude_bits)) - 1;
    const minimum_normal_exponent = 1 - format.exponent_bias;
    const minimum_subnormal_exponent = minimum_normal_exponent -
        @as(i32, @intCast(format.fraction_bits));

    if (top_exponent < minimum_normal_exponent) {
        const subnormal_shift = exponent - minimum_subnormal_exponent;
        stdAssert(subnormal_shift >= 0);
        const fraction: format.Bits = @intCast(magnitude << @intCast(subnormal_shift));
        stdAssert(fraction != 0);
        stdAssert(fraction & ~fraction_mask == 0);
        return format.Ieee.fromBits(sign_bits | fraction);
    }

    const significand_shift = format.precision - magnitude_bits;
    const significand: format.Bits = @intCast(magnitude << @intCast(significand_shift));
    const exponent_bits: format.Bits = @intCast(top_exponent + format.exponent_bias);
    const fraction = significand & fraction_mask;
    stdAssert(exponent_bits != 0);
    stdAssert(
        exponent_bits & ~@as(format.Bits, @intCast(exponent_mask >> format.fraction_bits)) == 0,
    );
    return format.Ieee.fromBits(sign_bits | (exponent_bits << format.fraction_bits) | fraction);
}

fn reduce(
    comptime format: Format,
    x: format.Float,
    y: format.Float,
    reduction: Reduction,
) types.RemquoResult(format.Float) {
    if (format.Ieee.firstNaN2(x, y)) |nan| {
        return .{ .value = nan, .quotient = 0 };
    }

    const x_bits = format.Ieee.toBits(x);
    const y_bits = format.Ieee.toBits(y);
    const x_class = format.Ieee.classifyBits(x_bits);
    const y_class = format.Ieee.classifyBits(y_bits);

    if (x_class == .infinity or y_class == .zero) {
        return .{ .value = format.Ieee.canonicalNaN(), .quotient = 0 };
    }
    if (y_class == .infinity or x_class == .zero) {
        return .{ .value = x, .quotient = 0 };
    }

    const x_parts = decompose(format, x);
    const y_parts = decompose(format, y);
    const quotient_negative = x_parts.negative != y_parts.negative;
    const exponent_delta = x_parts.exponent - y_parts.exponent;

    if (exponent_delta < 0) {
        // The higher-exponent divisor is normal: its shifted significand is
        // larger than the dividend's. Only a one-step gap can exceed half.
        if (reduction == .truncating or exponent_delta < -1 or
            x_parts.mantissa <= y_parts.mantissa)
        {
            return .{ .value = x, .quotient = 0 };
        }
        const magnitude = (y_parts.mantissa << 1) - x_parts.mantissa;
        return .{
            .value = packExact(format, magnitude, x_parts.exponent, !x_parts.negative),
            .quotient = if (quotient_negative) -1 else 1,
        };
    }

    const division = divideShifted(x_parts.mantissa, y_parts.mantissa, @intCast(exponent_delta));

    var quotient_low = division.quotient_low;
    var remainder_magnitude = division.remainder;
    var remainder_negative = x_parts.negative;
    if (reduction == .nearest_even) {
        const denominator = y_parts.mantissa;
        const twice_remainder = remainder_magnitude << 1;
        const increment = twice_remainder > denominator or
            (twice_remainder == denominator and quotient_low & 1 != 0);
        if (increment) {
            remainder_magnitude = denominator - remainder_magnitude;
            quotient_low = (quotient_low + 1) & 7;
            remainder_negative = !remainder_negative;
        }
    }
    var quotient: i32 = @intCast(quotient_low);
    if (quotient_negative) quotient = -quotient;
    const value = packExact(format, remainder_magnitude, y_parts.exponent, remainder_negative);
    return .{ .value = value, .quotient = quotient };
}

/// Returns the truncated quotient of `x` divided by `y`, that is the unique
/// value of the form `x - n * y` whose magnitude is below `y`'s and whose sign
/// is the sign of `x`.
///
/// An infinite `x` or a zero `y` returns the binary16 positive canonical quiet
/// NaN, and a NaN operand is returned with its quiet bit set and its sign and
/// payload kept, `x` first when both are NaN. A zero `x` and an infinite `y`
/// return `x` unchanged.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn fmodh(x: f16, y: f16) f16 {
    return reduce(h_format, x, y, .truncating).value;
}

/// Returns the truncated quotient of `x` divided by `y` at binary32 width, with
/// the same domain, NaN and signed-zero rules as `fmodh`.
pub fn fmodf(x: f32, y: f32) f32 {
    return reduce(f_format, x, y, .truncating).value;
}

/// Returns the truncated quotient of `x` divided by `y` at binary64 width, with
/// the same domain, NaN and signed-zero rules as `fmodh`.
pub fn fmodd(x: f64, y: f64) f64 {
    return reduce(d_format, x, y, .truncating).value;
}

/// Returns the IEEE remainder of `x` divided by `y`, the quotient rounded to
/// the nearest integer with ties to even, so the result's magnitude is at most
/// half of `y`'s.
///
/// An exact tie keeps the quotient even, and an increment flips the remainder's
/// sign, so the result carries the sign of `x` unless the rounding changed the
/// quotient by one. An infinite `x` or a zero `y` returns the binary16 positive
/// canonical quiet NaN, and a NaN operand is returned with its quiet bit set and
/// its sign and payload kept, `x` first when both are NaN.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn remainderh(x: f16, y: f16) f16 {
    return reduce(h_format, x, y, .nearest_even).value;
}

/// Returns the IEEE remainder of `x` divided by `y` at binary32 width, with the
/// same tie, sign, domain and NaN rules as `remainderh`.
pub fn remainderf(x: f32, y: f32) f32 {
    return reduce(f_format, x, y, .nearest_even).value;
}

/// Returns the IEEE remainder of `x` divided by `y` at binary64 width, with the
/// same tie, sign, domain and NaN rules as `remainderh`.
pub fn remainderd(x: f64, y: f64) f64 {
    return reduce(d_format, x, y, .nearest_even).value;
}

/// Returns the IEEE remainder of `x` divided by `y` and the low three bits of
/// the integer quotient, as `remainderh` defines `value`.
///
/// `quotient` is signed and lies in `[-7, 7]`, its sign being the sign of the
/// quotient, and it is set on every path: it is `0` when the operation is
/// invalid, which is when `value` is the positive canonical quiet NaN, and `0`
/// when a NaN operand is propagated.
pub fn remquoh(x: f16, y: f16) types.RemquoResult(f16) {
    return reduce(h_format, x, y, .nearest_even);
}

/// Returns the IEEE remainder and the low quotient bits at binary32 width, with
/// the same rules as `remquoh`.
pub fn remquof(x: f32, y: f32) types.RemquoResult(f32) {
    return reduce(f_format, x, y, .nearest_even);
}

/// Returns the IEEE remainder and the low quotient bits at binary64 width, with
/// the same rules as `remquoh`.
pub fn remquod(x: f64, y: f64) types.RemquoResult(f64) {
    return reduce(d_format, x, y, .nearest_even);
}
