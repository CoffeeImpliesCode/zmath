// Correctly rounded true gamma function for binary16, binary32, and binary64.
//
// The finite kernels are derived directly from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and are used under the MIT License.
// The wrappers construct exceptional results from target-format bits so that
// NaN payloads, signs, and poles do not depend on host library conventions.

extern fn cr_tgammaf16(value: f16) callconv(.c) f16;
extern fn cr_tgammaf(value: f32) callconv(.c) f32;
extern fn cr_tgamma(value: f64) callconv(.c) f64;

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

fn isNegativeIntegerH(bits: u16) bool {
    const exponent_field: u16 = (bits >> 10) & 0x1f;
    const fraction: u16 = bits & 0x03ff;
    if ((bits & sign_mask_h) == 0 or exponent_field == 0) return false;
    const unbiased_exponent: i32 = @as(i32, exponent_field) - 15;
    if (unbiased_exponent < 0) return false;
    if (unbiased_exponent >= 10) return true;
    const fractional_mask: u16 = (@as(u16, 1) << @as(u4, @intCast(10 - unbiased_exponent))) - 1;
    return (fraction & fractional_mask) == 0;
}

fn isNegativeIntegerF(bits: u32) bool {
    const exponent_field: u32 = (bits >> 23) & 0xff;
    const fraction: u32 = bits & 0x007f_ffff;
    if ((bits & sign_mask_f) == 0 or exponent_field == 0) return false;
    const unbiased_exponent: i32 = @as(i32, @intCast(exponent_field)) - 127;
    if (unbiased_exponent < 0) return false;
    if (unbiased_exponent >= 23) return true;
    const fractional_mask: u32 = (@as(u32, 1) << @as(u5, @intCast(23 - unbiased_exponent))) - 1;
    return (fraction & fractional_mask) == 0;
}

fn isNegativeIntegerD(bits: u64) bool {
    const exponent_field: u64 = (bits >> 52) & 0x7ff;
    const fraction: u64 = bits & 0x000f_ffff_ffff_ffff;
    if ((bits & sign_mask_d) == 0 or exponent_field == 0) return false;
    const unbiased_exponent: i32 = @as(i32, @intCast(exponent_field)) - 1023;
    if (unbiased_exponent < 0) return false;
    if (unbiased_exponent >= 52) return true;
    const fractional_mask: u64 = (@as(u64, 1) << @as(u6, @intCast(52 - unbiased_exponent))) - 1;
    return (fraction & fractional_mask) == 0;
}

/// Computes the true gamma function of a binary16 value.
pub fn tgammah(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    const magnitude: u16 = bits & magnitude_mask_h;
    if (magnitude > infinity_bits_h) return @bitCast(bits | quiet_nan_bit_h);
    if (magnitude == infinity_bits_h) {
        if ((bits & sign_mask_h) != 0) return @bitCast(canonical_nan_bits_h);
        return value;
    }
    if (magnitude == 0) {
        const result_bits: u16 = if ((bits & sign_mask_h) != 0) 0xfc00 else 0x7c00;
        return @bitCast(result_bits);
    }
    if (isNegativeIntegerH(bits)) return @bitCast(canonical_nan_bits_h);
    return cr_tgammaf16(value);
}

/// Computes the true gamma function of a binary32 value.
pub fn tgammaf(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const magnitude: u32 = bits & magnitude_mask_f;
    if (magnitude > infinity_bits_f) return @bitCast(bits | quiet_nan_bit_f);
    if (magnitude == infinity_bits_f) {
        if ((bits & sign_mask_f) != 0) return @bitCast(canonical_nan_bits_f);
        return value;
    }
    if (magnitude == 0) {
        const result_bits: u32 = if ((bits & sign_mask_f) != 0) 0xff80_0000 else 0x7f80_0000;
        return @bitCast(result_bits);
    }
    if (isNegativeIntegerF(bits)) return @bitCast(canonical_nan_bits_f);
    return cr_tgammaf(value);
}

/// Computes the true gamma function of a binary64 value.
pub fn tgammad(value: f64) f64 {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    const magnitude: u64 = bits & magnitude_mask_d;
    if (magnitude > infinity_bits_d) return @bitCast(bits | quiet_nan_bit_d);
    if (magnitude == infinity_bits_d) {
        if ((bits & sign_mask_d) != 0) return @bitCast(canonical_nan_bits_d);
        return value;
    }
    if (magnitude == 0) {
        const result_bits: u64 =
            if ((bits & sign_mask_d) != 0) 0xfff0_0000_0000_0000 else 0x7ff0_0000_0000_0000;
        return @bitCast(result_bits);
    }
    if (isNegativeIntegerD(bits)) return @bitCast(canonical_nan_bits_d);
    return cr_tgamma(value);
}
