// Correctly rounded power functions for binary16, binary32, and binary64.
//
// The finite kernels are derived directly from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and are used under the MIT License.
// The wrappers construct exceptional results from target-format bits so that
// NaN payloads, signs, and C23 identity and pole rules are deterministic.

extern fn cr_powf16(base: f16, exponent: f16) callconv(.c) f16;
extern fn cr_powf(base: f32, exponent: f32) callconv(.c) f32;
extern fn cr_pow(base: f64, exponent: f64) callconv(.c) f64;

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

fn isIntegerH(bits: u16) bool {
    const magnitude = bits & magnitude_mask_h;
    const exponent_field = (magnitude >> 10) & 0x1f;
    const fraction = magnitude & 0x03ff;
    if (exponent_field == 0) return magnitude == 0;
    const exponent: i32 = @as(i32, exponent_field) - 15;
    if (exponent < 0) return false;
    if (exponent >= 10) return true;
    const mask: u16 = (@as(u16, 1) << @as(u4, @intCast(10 - exponent))) - 1;
    return (fraction & mask) == 0;
}

fn isOddIntegerH(bits: u16) bool {
    const magnitude = bits & magnitude_mask_h;
    const exponent_field = (magnitude >> 10) & 0x1f;
    const fraction = magnitude & 0x03ff;
    if (exponent_field == 0) return false;
    const exponent: i32 = @as(i32, exponent_field) - 15;
    if (exponent < 0 or exponent > 10) return false;
    if (!isIntegerH(bits)) return false;
    const significand = @as(u16, 0x0400) | fraction;
    return ((significand >> @as(u4, @intCast(10 - exponent))) & 1) != 0;
}

fn isIntegerF(bits: u32) bool {
    const magnitude = bits & magnitude_mask_f;
    const exponent_field = (magnitude >> 23) & 0xff;
    const fraction = magnitude & 0x007f_ffff;
    if (exponent_field == 0) return magnitude == 0;
    const exponent: i32 = @as(i32, @intCast(exponent_field)) - 127;
    if (exponent < 0) return false;
    if (exponent >= 23) return true;
    const mask: u32 = (@as(u32, 1) << @as(u5, @intCast(23 - exponent))) - 1;
    return (fraction & mask) == 0;
}

fn isOddIntegerF(bits: u32) bool {
    const magnitude = bits & magnitude_mask_f;
    const exponent_field = (magnitude >> 23) & 0xff;
    const fraction = magnitude & 0x007f_ffff;
    if (exponent_field == 0) return false;
    const exponent: i32 = @as(i32, @intCast(exponent_field)) - 127;
    if (exponent < 0 or exponent > 23) return false;
    if (!isIntegerF(bits)) return false;
    const significand = @as(u32, 0x0080_0000) | fraction;
    return ((significand >> @as(u5, @intCast(23 - exponent))) & 1) != 0;
}

fn isIntegerD(bits: u64) bool {
    const magnitude = bits & magnitude_mask_d;
    const exponent_field = (magnitude >> 52) & 0x7ff;
    const fraction = magnitude & 0x000f_ffff_ffff_ffff;
    if (exponent_field == 0) return magnitude == 0;
    const exponent: i32 = @as(i32, @intCast(exponent_field)) - 1023;
    if (exponent < 0) return false;
    if (exponent >= 52) return true;
    const mask: u64 = (@as(u64, 1) << @as(u6, @intCast(52 - exponent))) - 1;
    return (fraction & mask) == 0;
}

fn isOddIntegerD(bits: u64) bool {
    const magnitude = bits & magnitude_mask_d;
    const exponent_field = (magnitude >> 52) & 0x7ff;
    const fraction = magnitude & 0x000f_ffff_ffff_ffff;
    if (exponent_field == 0) return false;
    const exponent: i32 = @as(i32, @intCast(exponent_field)) - 1023;
    if (exponent < 0 or exponent > 52) return false;
    if (!isIntegerD(bits)) return false;
    const significand = @as(u64, 0x0010_0000_0000_0000) | fraction;
    return ((significand >> @as(u6, @intCast(52 - exponent))) & 1) != 0;
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
        if ((base_bits & sign_mask_h) != 0 and !isIntegerH(exponent_bits)) {
            return @bitCast(canonical_nan_bits_h);
        }
        if ((base_bits & sign_mask_h) != 0 and isOddIntegerH(exponent_bits)) {
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
        const negative_result = (base_bits & sign_mask_h) != 0 and isOddIntegerH(exponent_bits);
        if (exponent_negative) {
            const result_bits: u16 = if (negative_result) 0xfc00 else 0x7c00;
            return @bitCast(result_bits);
        }
        const result_bits: u16 = if (negative_result) 0x8000 else 0x0000;
        return @bitCast(result_bits);
    }
    if (base_magnitude == infinity_bits_h) {
        const negative_result = (base_bits & sign_mask_h) != 0 and isOddIntegerH(exponent_bits);
        const result_sign = if (negative_result) sign_mask_h else 0;
        return @bitCast(if ((exponent_bits & sign_mask_h) != 0) result_sign else result_sign | infinity_bits_h);
    }
    if ((base_bits & sign_mask_h) != 0) {
        if (!isIntegerH(exponent_bits)) return @bitCast(canonical_nan_bits_h);
        return @bitCast(@as(u16, @bitCast(cr_powf16(@bitCast(base_bits & magnitude_mask_h), exponent))) |
            if (isOddIntegerH(exponent_bits)) sign_mask_h else 0);
    }
    return cr_powf16(base, exponent);
}

/// Computes the correctly rounded power of two binary32 values.
pub fn powf(base: f32, exponent: f32) f32 {
    @setFloatMode(.strict);
    const base_bits: u32 = @bitCast(base);
    const exponent_bits: u32 = @bitCast(exponent);
    const base_magnitude = base_bits & magnitude_mask_f;
    const exponent_magnitude = exponent_bits & magnitude_mask_f;

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
        if ((base_bits & sign_mask_f) != 0 and !isIntegerF(exponent_bits)) {
            return @bitCast(canonical_nan_bits_f);
        }
        if ((base_bits & sign_mask_f) != 0 and isOddIntegerF(exponent_bits)) {
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
        const negative_result = (base_bits & sign_mask_f) != 0 and isOddIntegerF(exponent_bits);
        if (exponent_negative) {
            const result_bits: u32 = if (negative_result) 0xff80_0000 else 0x7f80_0000;
            return @bitCast(result_bits);
        }
        const result_bits: u32 = if (negative_result) 0x8000_0000 else 0x0000_0000;
        return @bitCast(result_bits);
    }
    if (base_magnitude == infinity_bits_f) {
        const negative_result = (base_bits & sign_mask_f) != 0 and isOddIntegerF(exponent_bits);
        const result_sign = if (negative_result) sign_mask_f else 0;
        return @bitCast(if ((exponent_bits & sign_mask_f) != 0) result_sign else result_sign | infinity_bits_f);
    }
    if ((base_bits & sign_mask_f) != 0) {
        if (!isIntegerF(exponent_bits)) return @bitCast(canonical_nan_bits_f);
        const magnitude_result: u32 = @bitCast(cr_powf(@bitCast(base_bits & magnitude_mask_f), exponent));
        return @bitCast(magnitude_result | if (isOddIntegerF(exponent_bits)) sign_mask_f else 0);
    }
    return cr_powf(base, exponent);
}

/// Computes the correctly rounded power of two binary64 values.
pub fn powd(base: f64, exponent: f64) f64 {
    @setFloatMode(.strict);
    const base_bits: u64 = @bitCast(base);
    const exponent_bits: u64 = @bitCast(exponent);
    const base_magnitude = base_bits & magnitude_mask_d;
    const exponent_magnitude = exponent_bits & magnitude_mask_d;

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
        if ((base_bits & sign_mask_d) != 0 and !isIntegerD(exponent_bits)) {
            return @bitCast(canonical_nan_bits_d);
        }
        if ((base_bits & sign_mask_d) != 0 and isOddIntegerD(exponent_bits)) {
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
        const negative_result = (base_bits & sign_mask_d) != 0 and isOddIntegerD(exponent_bits);
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
        const negative_result = (base_bits & sign_mask_d) != 0 and isOddIntegerD(exponent_bits);
        const result_sign = if (negative_result) sign_mask_d else 0;
        return @bitCast(if ((exponent_bits & sign_mask_d) != 0) result_sign else result_sign | infinity_bits_d);
    }
    if ((base_bits & sign_mask_d) != 0) {
        if (!isIntegerD(exponent_bits)) return @bitCast(canonical_nan_bits_d);
        const magnitude_result: u64 = @bitCast(cr_pow(@bitCast(base_bits & magnitude_mask_d), exponent));
        return @bitCast(magnitude_result | if (isOddIntegerD(exponent_bits)) sign_mask_d else 0);
    }
    return cr_pow(base, exponent);
}
