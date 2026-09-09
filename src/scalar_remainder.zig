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

pub fn fmodh(x: f16, y: f16) f16 {
    return reduce(h_format, x, y, .truncating).value;
}

pub fn fmodf(x: f32, y: f32) f32 {
    return reduce(f_format, x, y, .truncating).value;
}

pub fn fmodd(x: f64, y: f64) f64 {
    return reduce(d_format, x, y, .truncating).value;
}

pub fn remainderh(x: f16, y: f16) f16 {
    return reduce(h_format, x, y, .nearest_even).value;
}

pub fn remainderf(x: f32, y: f32) f32 {
    return reduce(f_format, x, y, .nearest_even).value;
}

pub fn remainderd(x: f64, y: f64) f64 {
    return reduce(d_format, x, y, .nearest_even).value;
}

pub fn remquoh(x: f16, y: f16) types.RemquoResult(f16) {
    return reduce(h_format, x, y, .nearest_even);
}

pub fn remquof(x: f32, y: f32) types.RemquoResult(f32) {
    return reduce(f_format, x, y, .nearest_even);
}

pub fn remquod(x: f64, y: f64) types.RemquoResult(f64) {
    return reduce(d_format, x, y, .nearest_even);
}
