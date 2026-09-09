// Exact binary scaling for binary16, binary32, and binary64.
//
// The implementation operates on IEEE encodings so that gradual underflow,
// ties-to-even, signed special values, and extreme integer exponents do not
// depend on host floating-point environment or libm behavior. Callers use
// round-to-nearest, ties-to-even for this C23 surface.

const std = @import("std");

const HFormat = struct {
    const Float = f16;
    const Bits = u16;
    const sign_mask: Bits = 0x8000;
    const exponent_mask: Bits = 0x7c00;
    const fraction_mask: Bits = 0x03ff;
    const bias: i64 = 15;
};

const FFormat = struct {
    const Float = f32;
    const Bits = u32;
    const sign_mask: Bits = 0x8000_0000;
    const exponent_mask: Bits = 0x7f80_0000;
    const fraction_mask: Bits = 0x007f_ffff;
    const bias: i64 = 127;
};

const DFormat = struct {
    const Float = f64;
    const Bits = u64;
    const sign_mask: Bits = 0x8000_0000_0000_0000;
    const exponent_mask: Bits = 0x7ff0_0000_0000_0000;
    const fraction_mask: Bits = 0x000f_ffff_ffff_ffff;
    const bias: i64 = 1023;
};

fn addExponent(base: i64, delta: i64) i64 {
    if (delta > 0) {
        if (base > std.math.maxInt(i64) - delta) return std.math.maxInt(i64);
    } else if (delta < 0) {
        if (base < std.math.minInt(i64) - delta) return std.math.minInt(i64);
    }
    return base + delta;
}

fn scaleBits(comptime Format: type, value: Format.Float, exponent: i64) Format.Float {
    const Bits = Format.Bits;
    const Shift = std.math.Log2Int(Bits);
    const fraction_bits: u6 = @intCast(@ctz(Format.exponent_mask));
    const precision: i64 = @as(i64, fraction_bits) + 1;
    const hidden_bit: Bits = @as(Bits, 1) << fraction_bits;
    const quiet_bit: Bits = hidden_bit >> 1;
    const magnitude_mask: Bits = ~Format.sign_mask;
    const infinity_bits: Bits = Format.exponent_mask;
    const input_bits: Bits = @bitCast(value);
    const sign_bits: Bits = input_bits & Format.sign_mask;
    const magnitude_bits: Bits = input_bits & magnitude_mask;
    const exponent_bits: Bits = magnitude_bits & Format.exponent_mask;
    const fraction: Bits = magnitude_bits & Format.fraction_mask;

    if (exponent_bits == 0) {
        if (fraction == 0) return @bitCast(input_bits);
    } else if (exponent_bits == Format.exponent_mask) {
        if (fraction == 0) return @bitCast(input_bits);
        return @bitCast(input_bits | quiet_bit);
    }

    const minimum_normal_exponent: i64 = 1 - Format.bias;
    const minimum_subnormal_exponent: i64 = minimum_normal_exponent - @as(i64, fraction_bits);
    const maximum_normal_exponent: i64 = Format.bias;
    const minimum_normal_scale: i64 = minimum_subnormal_exponent;
    const maximum_normal_scale: i64 = maximum_normal_exponent - @as(i64, fraction_bits);

    var significand: Bits = fraction;
    var scale_exponent: i64 = minimum_subnormal_exponent;
    if (exponent_bits != 0) {
        significand |= hidden_bit;
        const unbiased_exponent: i64 =
            @as(i64, @intCast(exponent_bits >> fraction_bits)) - Format.bias;
        scale_exponent = unbiased_exponent - @as(i64, fraction_bits);
    } else {
        // Zero returned above, so 0 < significand < hidden_bit. The CLZ
        // difference is 1..fraction_bits and normalizes without losing bits.
        // The exponent decrement preserves the value and bounds scale_exponent
        // to [minimum_subnormal_exponent - fraction_bits, minimum_subnormal_exponent - 1].
        const normalization_shift: Shift = @intCast(@clz(significand) - @clz(hidden_bit));
        significand <<= normalization_shift;
        scale_exponent -= @as(i64, normalization_shift);
    }

    // Finite inputs bound scale_exponent to
    // [minimum_subnormal_exponent - fraction_bits, maximum_normal_scale].
    const result_scale: i64 = addExponent(scale_exponent, exponent);
    if (result_scale > maximum_normal_scale) {
        return @bitCast(sign_bits | infinity_bits);
    }

    if (result_scale >= minimum_normal_scale) {
        const result_exponent: Bits = @intCast(
            result_scale + @as(i64, fraction_bits) + Format.bias,
        );
        const result_bits: Bits = sign_bits |
            (result_exponent << fraction_bits) |
            (significand & Format.fraction_mask);
        return @bitCast(result_bits);
    }

    const shift: i64 = minimum_subnormal_exponent - result_scale;
    var rounded_significand: Bits = 0;
    if (shift < precision) {
        const shift_bits: Shift = @intCast(shift);
        rounded_significand = significand >> shift_bits;
        const remainder_mask: Bits = (@as(Bits, 1) << shift_bits) - 1;
        const remainder: Bits = significand & remainder_mask;
        const halfway_shift: Shift = @intCast(shift - 1);
        const halfway: Bits = @as(Bits, 1) << halfway_shift;
        if (remainder > halfway or
            (remainder == halfway and (rounded_significand & 1) != 0))
        {
            rounded_significand += 1;
        }
    } else if (shift == precision) {
        if (significand > hidden_bit) rounded_significand = 1;
    }

    if (rounded_significand >= hidden_bit) {
        return @bitCast(sign_bits | hidden_bit);
    }
    return @bitCast(sign_bits | rounded_significand);
}

/// Scales a binary16 value by an exact power of two.
pub fn ldexph(value: f16, exponent: i32) f16 {
    @setFloatMode(.strict);
    return scaleBits(HFormat, value, @as(i64, exponent));
}

/// Scales a binary32 value by an exact power of two.
pub fn ldexpf(value: f32, exponent: i32) f32 {
    @setFloatMode(.strict);
    return scaleBits(FFormat, value, @as(i64, exponent));
}

/// Scales a binary64 value by an exact power of two.
pub fn ldexpd(value: f64, exponent: i32) f64 {
    @setFloatMode(.strict);
    return scaleBits(DFormat, value, @as(i64, exponent));
}

/// Scales a binary16 value by an exact power of two.
pub fn scalbnh(value: f16, exponent: i32) f16 {
    @setFloatMode(.strict);
    return scaleBits(HFormat, value, @as(i64, exponent));
}

/// Scales a binary32 value by an exact power of two.
pub fn scalbnf(value: f32, exponent: i32) f32 {
    @setFloatMode(.strict);
    return scaleBits(FFormat, value, @as(i64, exponent));
}

/// Scales a binary64 value by an exact power of two.
pub fn scalbnd(value: f64, exponent: i32) f64 {
    @setFloatMode(.strict);
    return scaleBits(DFormat, value, @as(i64, exponent));
}

/// Scales a binary16 value by an exact long integer power of two.
pub fn scalblnh(value: f16, exponent: i64) f16 {
    @setFloatMode(.strict);
    return scaleBits(HFormat, value, exponent);
}

/// Scales a binary32 value by an exact long integer power of two.
pub fn scalblnf(value: f32, exponent: i64) f32 {
    @setFloatMode(.strict);
    return scaleBits(FFormat, value, exponent);
}

/// Scales a binary64 value by an exact long integer power of two.
pub fn scalblnd(value: f64, exponent: i64) f64 {
    @setFloatMode(.strict);
    return scaleBits(DFormat, value, exponent);
}
