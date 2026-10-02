// Correctly rounded base-10 logarithm of one plus a value.
//
// The finite kernels are derived directly from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and are used under the MIT License.
// The wrappers construct exceptional results from target-format bits so that
// NaN payloads, signs, domains, and poles do not depend on host conventions.

const tiny = @import("scalar_tiny.zig");

extern fn zmath_core_log10p1h_bits(value: u16) callconv(.c) u16;
extern fn zmath_core_log10p1f_finite(value: f32) callconv(.c) f32;
extern fn zmath_core_log10p1_finite(value: f64) callconv(.c) f64;

const sign_mask_h: u16 = 0x8000;
const sign_mask_f: u32 = 0x8000_0000;
const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_h: u16 = 0x7fff;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_h: u16 = 0x7c00;
const infinity_bits_f: u32 = 0x7f80_0000;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const negative_one_bits_h: u16 = 0xbc00;
const negative_one_bits_f: u32 = 0xbf80_0000;
const negative_one_bits_d: u64 = 0xbff0_0000_0000_0000;
const quiet_nan_bit_h: u16 = 0x0200;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits_h: u16 = 0x7e00;
const canonical_nan_bits_f: u32 = 0x7fc0_0000;
const canonical_nan_bits_d: u64 = 0x7ff8_0000_0000_0000;

/// Computes the correctly rounded base-10 logarithm of one plus a binary16 value.
/// The caller must keep IEEE round-to-nearest, ties-to-even active.
pub fn log10p1h(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;
    if (magnitude > infinity_bits_h) return @bitCast(bits | quiet_nan_bit_h);
    if (magnitude == infinity_bits_h) {
        if ((bits & sign_mask_h) != 0) return @bitCast(canonical_nan_bits_h);
        return value;
    }
    if (bits == negative_one_bits_h) return @bitCast(@as(u16, 0xfc00));
    if ((bits & sign_mask_h) != 0 and bits > negative_one_bits_h) {
        return @bitCast(canonical_nan_bits_h);
    }
    if (magnitude == 0) return value;
    return @bitCast(zmath_core_log10p1h_bits(bits));
}

/// Computes the correctly rounded base-10 logarithm of one plus a binary32 value.
/// The caller must keep IEEE round-to-nearest, ties-to-even active.
pub fn log10p1f(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask_f;
    if (magnitude > infinity_bits_f) return @bitCast(bits | quiet_nan_bit_f);
    if (magnitude == infinity_bits_f) {
        if ((bits & sign_mask_f) != 0) return @bitCast(canonical_nan_bits_f);
        return value;
    }
    if (bits == negative_one_bits_f) return @bitCast(@as(u32, 0xff80_0000));
    if ((bits & sign_mask_f) != 0 and bits > negative_one_bits_f) {
        return @bitCast(canonical_nan_bits_f);
    }
    if (magnitude == 0) return value;
    return zmath_core_log10p1f_finite(value);
}

/// Computes the correctly rounded base-10 logarithm of one plus a binary64 value.
/// The caller must keep IEEE round-to-nearest, ties-to-even active.
pub fn log10p1d(value: f64) f64 {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    const magnitude = bits & magnitude_mask_d;
    if (magnitude >= infinity_bits_d or bits >= negative_one_bits_d) {
        @branchHint(.unlikely);
        if (magnitude > infinity_bits_d) return @bitCast(bits | quiet_nan_bit_d);
        if (bits == negative_one_bits_d) return @bitCast(@as(u64, 0xfff0_0000_0000_0000));
        if (bits > negative_one_bits_d) return @bitCast(canonical_nan_bits_d);
        return value;
    }
    if (magnitude < tiny.magnitude_limit) return tiny.scale(.inv_ln10, value);
    return zmath_core_log10p1_finite(value);
}
