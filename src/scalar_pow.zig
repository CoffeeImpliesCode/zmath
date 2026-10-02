// Correctly rounded power functions for binary16, binary32, and binary64.
//
// The finite kernels are derived directly from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and are used under the MIT License.
// The wrappers construct exceptional results from target-format bits so that
// NaN payloads, signs, and C23 identity and pole rules are deterministic.

const std = @import("std");

extern fn zmath_core_powh_bits(base: u16, exponent: u16) callconv(.c) u16;
extern fn cr_powfFiniteBits(base_bits: u32, exponent: f32) callconv(.c) f32;
extern fn cr_powFinite(base: f64, exponent: f64) callconv(.c) f64;

const sign_mask_h: u16 = 0x8000;
const sign_mask_f: u32 = 0x8000_0000;
const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_h: u16 = 0x7fff;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_h: u16 = 0x7c00;
const infinity_bits_f: u32 = 0x7f80_0000;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_h: u16 = 0x0200;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits_h: u16 = 0x7e00;
const canonical_nan_bits_f: u32 = 0x7fc0_0000;
const canonical_nan_bits_d: u64 = 0x7ff8_0000_0000_0000;

fn classifyFinitePowerF(base_magnitude: u32, exponent_bits: u32, exponent_magnitude: u32) ?u32 {
    if (exponent_magnitude -% 0x2f80_0001 < 0x1f7f_ffff) return null;
    if (exponent_magnitude <= 0x2f80_0000) return 0x3f80_0000;
    const base_less_than_one = base_magnitude < 0x3f80_0000;
    const exponent_negative = (exponent_bits & sign_mask_f) != 0;
    if (base_less_than_one) {
        if (exponent_negative) return infinity_bits_f;
        if (exponent_magnitude >= 0x4f80_0000) return 0;
        return null;
    }
    return if (exponent_negative) 0 else infinity_bits_f;
}

fn classifyFinitePowerD(base_magnitude: u64, exponent_bits: u64) ?u64 {
    const exponent_magnitude = exponent_bits & magnitude_mask_d;
    if (exponent_magnitude <= 0x3bf0_0000_0000_0000) return 0x3ff0_0000_0000_0000;
    if (exponent_magnitude < 0x43e0_0000_0000_0000) return null;
    const base_less_than_one = base_magnitude < 0x3ff0_0000_0000_0000;
    const exponent_negative = (exponent_bits & sign_mask_d) != 0;
    if (base_less_than_one) {
        if (exponent_negative) return infinity_bits_d;
        if (exponent_magnitude >= 0x43f0_0000_0000_0000) return 0;
        return null;
    }
    return if (exponent_negative) 0 else infinity_bits_d;
}
const IntegerClass = enum {
    non_integer,
    even_integer,
    odd_integer,
};

fn classifyInteger(
    comptime Bits: type,
    comptime fraction_bits: comptime_int,
    comptime exponent_bias: comptime_int,
    comptime magnitude_mask: Bits,
    bits: Bits,
) IntegerClass {
    const magnitude = bits & magnitude_mask;
    if (magnitude == 0) return .even_integer;

    const exponent_field = magnitude >> fraction_bits;
    if (exponent_field < exponent_bias) return .non_integer;

    const unit_exponent = exponent_bias + fraction_bits;
    if (exponent_field > unit_exponent) return .even_integer;

    const fraction_mask = (@as(Bits, 1) << fraction_bits) - 1;
    const significand = (@as(Bits, 1) << fraction_bits) | (magnitude & fraction_mask);
    const shift: std.math.Log2Int(Bits) = @intCast(unit_exponent - exponent_field);
    const truncated_mask = (@as(Bits, 1) << shift) - 1;
    if ((significand & truncated_mask) != 0) return .non_integer;
    return if (((significand >> shift) & 1) != 0)
        .odd_integer
    else
        .even_integer;
}

fn classifyIntegerH(bits: u16) IntegerClass {
    return classifyInteger(u16, 10, 15, magnitude_mask_h, bits);
}

fn classifyIntegerF(bits: u32) IntegerClass {
    return classifyInteger(u32, 23, 127, magnitude_mask_f, bits);
}

fn classifyIntegerD(bits: u64) IntegerClass {
    return classifyInteger(u64, 52, 1023, magnitude_mask_d, bits);
}

/// Computes the correctly rounded power of two binary16 values.
pub fn powh(base: f16, exponent: f16) f16 {
    @setFloatMode(.strict);
    const base_bits: u16 = @bitCast(base);
    const exponent_bits: u16 = @bitCast(exponent);
    const base_magnitude = base_bits & magnitude_mask_h;
    const exponent_magnitude = exponent_bits & magnitude_mask_h;

    if (exponent_magnitude == 0) {
        if (base_magnitude > infinity_bits_h and (base_bits & quiet_nan_bit_h) == 0) {
            return @bitCast(base_bits | quiet_nan_bit_h);
        }
        return @bitCast(@as(u16, 0x3c00));
    }
    if (base_magnitude == 0x3c00) {
        if (exponent_magnitude > infinity_bits_h) {
            if ((base_bits & sign_mask_h) != 0 or (exponent_bits & quiet_nan_bit_h) == 0) {
                return @bitCast(exponent_bits | quiet_nan_bit_h);
            }
            return @bitCast(@as(u16, 0x3c00));
        }
        const exponent_class = classifyIntegerH(exponent_bits);
        if ((base_bits & sign_mask_h) != 0 and exponent_class == .non_integer) {
            return @bitCast(canonical_nan_bits_h);
        }
        if ((base_bits & sign_mask_h) != 0 and exponent_class == .odd_integer) {
            return @bitCast(@as(u16, 0xbc00));
        }
        return @bitCast(@as(u16, 0x3c00));
    }
    if (base_magnitude > infinity_bits_h) return @bitCast(base_bits | quiet_nan_bit_h);
    if (exponent_magnitude > infinity_bits_h) return @bitCast(exponent_bits | quiet_nan_bit_h);

    if (exponent_magnitude == infinity_bits_h) {
        const base_less_than_one = base_magnitude < 0x3c00;
        const exponent_negative = (exponent_bits & sign_mask_h) != 0;
        const result_bits: u16 = if (base_less_than_one == exponent_negative) 0x7c00 else 0x0000;
        return @bitCast(result_bits);
    }
    if (base_magnitude == 0) {
        const exponent_negative = (exponent_bits & sign_mask_h) != 0;
        const exponent_class = classifyIntegerH(exponent_bits);
        const negative_result = (base_bits & sign_mask_h) != 0 and exponent_class == .odd_integer;
        if (exponent_negative) {
            const result_bits: u16 = if (negative_result) 0xfc00 else 0x7c00;
            return @bitCast(result_bits);
        }
        const result_bits: u16 = if (negative_result) 0x8000 else 0x0000;
        return @bitCast(result_bits);
    }
    if (base_magnitude == infinity_bits_h) {
        const exponent_class = classifyIntegerH(exponent_bits);
        const negative_result =
            (base_bits & sign_mask_h) != 0 and exponent_class == .odd_integer;
        const result_sign = if (negative_result) sign_mask_h else 0;
        return @bitCast(if ((exponent_bits & sign_mask_h) != 0) result_sign else result_sign | infinity_bits_h);
    }
    if ((base_bits & sign_mask_h) != 0) {
        const exponent_class = classifyIntegerH(exponent_bits);
        if (exponent_class == .non_integer) return @bitCast(canonical_nan_bits_h);
        return @bitCast(zmath_core_powh_bits(base_bits & magnitude_mask_h, exponent_bits) |
            if (exponent_class == .odd_integer) sign_mask_h else 0);
    }
    return @bitCast(zmath_core_powh_bits(base_bits, exponent_bits));
}

noinline fn powfNonRegular(
    base_bits: u32,
    exponent_bits: u32,
    base_magnitude: u32,
    exponent_magnitude: u32,
    exponent: f32,
) f32 {
    if (exponent_magnitude == 0) {
        if (base_magnitude > infinity_bits_f and (base_bits & quiet_nan_bit_f) == 0) {
            return @bitCast(base_bits | quiet_nan_bit_f);
        }
        return @bitCast(@as(u32, 0x3f80_0000));
    }
    if (base_magnitude == 0x3f80_0000) {
        if (exponent_magnitude > infinity_bits_f) {
            if ((base_bits & sign_mask_f) != 0 or (exponent_bits & quiet_nan_bit_f) == 0) {
                return @bitCast(exponent_bits | quiet_nan_bit_f);
            }
            return @bitCast(@as(u32, 0x3f80_0000));
        }
        const exponent_class = classifyIntegerF(exponent_bits);
        if ((base_bits & sign_mask_f) != 0 and exponent_class == .non_integer) {
            return @bitCast(canonical_nan_bits_f);
        }
        if ((base_bits & sign_mask_f) != 0 and exponent_class == .odd_integer) {
            return @bitCast(@as(u32, 0xbf80_0000));
        }
        return @bitCast(@as(u32, 0x3f80_0000));
    }
    if (base_magnitude > infinity_bits_f) return @bitCast(base_bits | quiet_nan_bit_f);
    if (exponent_magnitude > infinity_bits_f) return @bitCast(exponent_bits | quiet_nan_bit_f);

    if (exponent_magnitude == infinity_bits_f) {
        const base_less_than_one = base_magnitude < 0x3f80_0000;
        const exponent_negative = (exponent_bits & sign_mask_f) != 0;
        const result_bits: u32 =
            if (base_less_than_one == exponent_negative) 0x7f80_0000 else 0x0000_0000;
        return @bitCast(result_bits);
    }
    if (base_magnitude == 0) {
        const exponent_negative = (exponent_bits & sign_mask_f) != 0;
        const exponent_class = classifyIntegerF(exponent_bits);
        const negative_result = (base_bits & sign_mask_f) != 0 and exponent_class == .odd_integer;
        if (exponent_negative) {
            const result_bits: u32 = if (negative_result) 0xff80_0000 else 0x7f80_0000;
            return @bitCast(result_bits);
        }
        const result_bits: u32 = if (negative_result) 0x8000_0000 else 0x0000_0000;
        return @bitCast(result_bits);
    }
    if (base_magnitude == infinity_bits_f) {
        const exponent_class = classifyIntegerF(exponent_bits);
        const negative_result = (base_bits & sign_mask_f) != 0 and exponent_class == .odd_integer;
        const result_sign = if (negative_result) sign_mask_f else 0;
        return @bitCast(if ((exponent_bits & sign_mask_f) != 0) result_sign else result_sign | infinity_bits_f);
    }
    if ((base_bits & sign_mask_f) != 0) {
        const exponent_class = classifyIntegerF(exponent_bits);
        if (exponent_class == .non_integer) return @bitCast(canonical_nan_bits_f);
        const shortcut_bits = classifyFinitePowerF(base_magnitude, exponent_bits, exponent_magnitude);
        const magnitude_result: u32 = shortcut_bits orelse
            @as(u32, @bitCast(cr_powfFiniteBits(base_magnitude, exponent)));
        return @bitCast(magnitude_result | if (exponent_class == .odd_integer) sign_mask_f else 0);
    }
    if (classifyFinitePowerF(base_magnitude, exponent_bits, exponent_magnitude)) |result_bits| {
        return @bitCast(result_bits);
    }
    return cr_powfFiniteBits(base_magnitude, exponent);
}

/// Computes the correctly rounded power of two binary32 values.
pub fn powf(base: f32, exponent: f32) f32 {
    @setFloatMode(.strict);
    const base_bits: u32 = @bitCast(base);
    const exponent_bits: u32 = @bitCast(exponent);
    const exponent_magnitude = exponent_bits & magnitude_mask_f;

    const positive_finite_nonzero_base =
        base_bits -% 1 < infinity_bits_f - 1;
    // This is exactly the exponent range where the classifier needs
    // the finite kernel. All identities and extreme exponents stay in
    // the complete nonregular dispatch below.
    const regular_exponent = exponent_magnitude -% 0x2f80_0001 < 0x1f7f_ffff;
    if (positive_finite_nonzero_base and regular_exponent and base_bits != 0x3f80_0000) {
        return cr_powfFiniteBits(base_bits, exponent);
    }

    return powfNonRegular(
        base_bits,
        exponent_bits,
        base_bits & magnitude_mask_f,
        exponent_magnitude,
        exponent,
    );
}

noinline fn powdNonRegular(
    base_bits: u64,
    exponent_bits: u64,
    base_magnitude: u64,
    exponent_magnitude: u64,
    base: f64,
    exponent: f64,
) f64 {
    if (exponent_magnitude == 0) {
        if (base_magnitude > infinity_bits_d and (base_bits & quiet_nan_bit_d) == 0) {
            return @bitCast(base_bits | quiet_nan_bit_d);
        }
        return @bitCast(@as(u64, 0x3ff0_0000_0000_0000));
    }
    if (base_magnitude == 0x3ff0_0000_0000_0000) {
        if (exponent_magnitude > infinity_bits_d) {
            if ((base_bits & sign_mask_d) != 0 or (exponent_bits & quiet_nan_bit_d) == 0) {
                return @bitCast(exponent_bits | quiet_nan_bit_d);
            }
            return @bitCast(@as(u64, 0x3ff0_0000_0000_0000));
        }
        const exponent_class = classifyIntegerD(exponent_bits);
        if ((base_bits & sign_mask_d) != 0 and exponent_class == .non_integer) {
            return @bitCast(canonical_nan_bits_d);
        }
        if ((base_bits & sign_mask_d) != 0 and exponent_class == .odd_integer) {
            return @bitCast(@as(u64, 0xbff0_0000_0000_0000));
        }
        return @bitCast(@as(u64, 0x3ff0_0000_0000_0000));
    }
    if (base_magnitude > infinity_bits_d) return @bitCast(base_bits | quiet_nan_bit_d);
    if (exponent_magnitude > infinity_bits_d) return @bitCast(exponent_bits | quiet_nan_bit_d);

    if (exponent_magnitude == infinity_bits_d) {
        const base_less_than_one = base_magnitude < 0x3ff0_0000_0000_0000;
        const exponent_negative = (exponent_bits & sign_mask_d) != 0;
        const result_bits: u64 =
            if (base_less_than_one == exponent_negative) 0x7ff0_0000_0000_0000 else 0x0000_0000_0000_0000;
        return @bitCast(result_bits);
    }
    if (base_magnitude == 0) {
        const exponent_negative = (exponent_bits & sign_mask_d) != 0;
        const exponent_class = classifyIntegerD(exponent_bits);
        const negative_result = (base_bits & sign_mask_d) != 0 and exponent_class == .odd_integer;
        if (exponent_negative) {
            const result_bits: u64 =
                if (negative_result) 0xfff0_0000_0000_0000 else 0x7ff0_0000_0000_0000;
            return @bitCast(result_bits);
        }
        const result_bits: u64 =
            if (negative_result) 0x8000_0000_0000_0000 else 0x0000_0000_0000_0000;
        return @bitCast(result_bits);
    }
    if (base_magnitude == infinity_bits_d) {
        const exponent_class = classifyIntegerD(exponent_bits);
        const negative_result = (base_bits & sign_mask_d) != 0 and exponent_class == .odd_integer;
        const result_sign = if (negative_result) sign_mask_d else 0;
        return @bitCast(if ((exponent_bits & sign_mask_d) != 0) result_sign else result_sign | infinity_bits_d);
    }
    if ((base_bits & sign_mask_d) != 0) {
        const exponent_class = classifyIntegerD(exponent_bits);
        if (exponent_class == .non_integer) return @bitCast(canonical_nan_bits_d);
        const shortcut_bits = classifyFinitePowerD(base_magnitude, exponent_bits);
        const magnitude_result: u64 = shortcut_bits orelse
            @as(u64, @bitCast(cr_powFinite(@bitCast(base_magnitude), exponent)));
        return @bitCast(magnitude_result | if (exponent_class == .odd_integer) sign_mask_d else 0);
    }
    if (classifyFinitePowerD(base_magnitude, exponent_bits)) |result_bits| {
        return @bitCast(result_bits);
    }
    return cr_powFinite(base, exponent);
}

/// Computes the correctly rounded power of two binary64 values.
pub fn powd(base: f64, exponent: f64) f64 {
    @setFloatMode(.strict);
    const base_bits: u64 = @bitCast(base);
    const exponent_bits: u64 = @bitCast(exponent);
    const exponent_magnitude = exponent_bits & magnitude_mask_d;

    const positive_finite_nonzero_base =
        base_bits -% 1 < infinity_bits_d - 1;
    // This is exactly the exponent range where the classifier always needs
    // the finite kernel. All identities and extreme exponents stay in the
    // complete nonregular dispatch below.
    const regular_exponent =
        exponent_magnitude -% 0x3bf0_0000_0000_0001 < 0x07ef_ffff_ffff_ffff;
    if (positive_finite_nonzero_base and
        regular_exponent and
        base_bits != 0x3ff0_0000_0000_0000)
    {
        return cr_powFinite(base, exponent);
    }

    return powdNonRegular(
        base_bits,
        exponent_bits,
        base_bits & magnitude_mask_d,
        exponent_magnitude,
        base,
        exponent,
    );
}
