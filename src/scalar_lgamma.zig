// Correctly rounded logarithm of the absolute gamma function for binary16, binary32, and
// binary64.
//
// The finite kernels are derived directly from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and are used under the MIT License.
// The wrappers construct exceptional results from target-format bits so that NaN payloads,
// signs, and poles do not depend on host conventions.

const c23 = @import("scalar_c23_types.zig");

/// The rounded logarithm of |gamma(x)| and the sign of gamma(x).
///
/// The sign is -1, 0, or 1. Callers must keep IEEE round-to-nearest,
/// ties-to-even active. The implementation does not inspect or change fenv.
pub const LgammaResult = c23.LgammaResult;

extern fn zmath_core_lgammah(value: f16, sign: *c_int) callconv(.c) f16;
extern fn zmath_core_lgammaf(value: f32, sign: *c_int) callconv(.c) f32;
extern fn zmath_core_lgammad(value: f64, sign: *c_int) callconv(.c) f64;

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

fn isNegativeIntegerH(bits: u16) bool {
    if ((bits & sign_mask_h) == 0) return false;

    const exponent_field: u16 = (bits >> 10) & 0x1f;
    if (exponent_field == 0 or exponent_field == 0x1f) return false;

    const unbiased_exponent: i32 = @as(i32, exponent_field) - 15;
    if (unbiased_exponent < 0) return false;
    if (unbiased_exponent >= 10) return true;

    const fractional_mask: u16 =
        (@as(u16, 1) << @as(u4, @intCast(10 - unbiased_exponent))) - 1;
    return (bits & fractional_mask) == 0;
}

fn isNegativeIntegerF(bits: u32) bool {
    if ((bits & sign_mask_f) == 0) return false;

    const exponent_field: u32 = (bits >> 23) & 0xff;
    if (exponent_field == 0 or exponent_field == 0xff) return false;

    const unbiased_exponent: i32 = @as(i32, @intCast(exponent_field)) - 127;
    if (unbiased_exponent < 0) return false;
    if (unbiased_exponent >= 23) return true;

    const fractional_mask: u32 =
        (@as(u32, 1) << @as(u5, @intCast(23 - unbiased_exponent))) - 1;
    return (bits & fractional_mask) == 0;
}

fn isNegativeIntegerD(bits: u64) bool {
    if ((bits & sign_mask_d) == 0) return false;

    const exponent_field: u64 = (bits >> 52) & 0x7ff;
    if (exponent_field == 0 or exponent_field == 0x7ff) return false;

    const unbiased_exponent: i32 = @as(i32, @intCast(exponent_field)) - 1023;
    if (unbiased_exponent < 0) return false;
    if (unbiased_exponent >= 52) return true;

    const fractional_mask: u64 =
        (@as(u64, 1) << @as(u6, @intCast(52 - unbiased_exponent))) - 1;
    return (bits & fractional_mask) == 0;
}

/// Computes the correctly rounded binary16 lgamma result.
pub fn lgammah(value: f16) c23.LgammaResult(f16) {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;

    if (magnitude > infinity_bits_h) {
        return .{ .value = @bitCast(bits | quiet_nan_bit_h), .sign = 1 };
    }
    if (magnitude == infinity_bits_h) {
        return .{ .value = @bitCast(infinity_bits_h), .sign = 1 };
    }
    if (magnitude == 0) {
        const sign: i2 = if ((bits & sign_mask_h) == 0) 1 else -1;
        return .{ .value = @bitCast(infinity_bits_h), .sign = sign };
    }
    if (isNegativeIntegerH(bits)) {
        return .{ .value = @bitCast(infinity_bits_h), .sign = -1 };
    }

    var sign_c: c_int = 1;
    const result = zmath_core_lgammah(value, &sign_c);
    const sign: i2 = @intCast(sign_c);
    return .{ .value = result, .sign = sign };
}

/// Computes the correctly rounded binary32 lgamma result.
pub fn lgammaf(value: f32) c23.LgammaResult(f32) {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask_f;

    if (magnitude > infinity_bits_f) {
        return .{ .value = @bitCast(bits | quiet_nan_bit_f), .sign = 1 };
    }
    if (magnitude == infinity_bits_f) {
        return .{ .value = @bitCast(infinity_bits_f), .sign = 1 };
    }
    if (magnitude == 0) {
        const sign: i2 = if ((bits & sign_mask_f) == 0) 1 else -1;
        return .{ .value = @bitCast(infinity_bits_f), .sign = sign };
    }
    if (isNegativeIntegerF(bits)) {
        return .{ .value = @bitCast(infinity_bits_f), .sign = -1 };
    }

    var sign_c: c_int = 1;
    const result = zmath_core_lgammaf(value, &sign_c);
    const sign: i2 = @intCast(sign_c);
    return .{ .value = result, .sign = sign };
}

/// Computes the correctly rounded binary64 lgamma result.
pub fn lgammad(value: f64) c23.LgammaResult(f64) {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    const magnitude = bits & magnitude_mask_d;

    if (magnitude > infinity_bits_d) {
        return .{ .value = @bitCast(bits | quiet_nan_bit_d), .sign = 1 };
    }
    if (magnitude == infinity_bits_d) {
        return .{ .value = @bitCast(infinity_bits_d), .sign = 1 };
    }
    if (magnitude == 0) {
        const sign: i2 = if ((bits & sign_mask_d) == 0) 1 else -1;
        return .{ .value = @bitCast(infinity_bits_d), .sign = sign };
    }
    if (isNegativeIntegerD(bits)) {
        return .{ .value = @bitCast(infinity_bits_d), .sign = -1 };
    }

    var sign_c: c_int = 1;
    const result = zmath_core_lgammad(value, &sign_c);
    const sign: i2 = @intCast(sign_c);
    return .{ .value = result, .sign = sign };
}
