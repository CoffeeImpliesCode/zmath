// Correctly rounded logarithm of the absolute gamma function for binary16, binary32, and
// binary64.
//
// The finite kernels are derived directly from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and are used under the MIT License.
// The kernels own the infinity, zero, and integer-pole results and their signs.
// The wrappers keep only the NaN path in Zig so that a NaN payload and sign
// never depend on host conventions.

const c23 = @import("scalar_c23_types.zig");

/// The rounded logarithm of |gamma(x)| and the sign of gamma(x).
///
/// The sign is -1, 0, or 1. Callers must keep IEEE round-to-nearest,
/// ties-to-even active. The implementation does not inspect or change fenv.
pub const LgammaResult = c23.LgammaResult;

extern fn zmath_core_lgammah_bits(value: u16, sign: *c_int) callconv(.c) u16;
extern fn zmath_core_lgammaf(value: f32, sign: *c_int) callconv(.c) f32;
extern fn zmath_core_lgammad(value: f64, sign: *c_int) callconv(.c) f64;

const magnitude_mask_h: u16 = 0x7fff;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_h: u16 = 0x7c00;
const infinity_bits_f: u32 = 0x7f80_0000;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_h: u16 = 0x0200;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;

/// Computes the correctly rounded binary16 lgamma result.
pub fn lgammah(value: f16) c23.LgammaResult(f16) {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    if (bits & magnitude_mask_h > infinity_bits_h) {
        return .{ .value = @bitCast(bits | quiet_nan_bit_h), .sign = 1 };
    }

    var sign_c: c_int = 1;
    const result = zmath_core_lgammah_bits(bits, &sign_c);
    return .{ .value = @bitCast(result), .sign = @intCast(sign_c) };
}

/// Computes the correctly rounded binary32 lgamma result.
pub fn lgammaf(value: f32) c23.LgammaResult(f32) {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    if (bits & magnitude_mask_f > infinity_bits_f) {
        return .{ .value = @bitCast(bits | quiet_nan_bit_f), .sign = 1 };
    }

    var sign_c: c_int = 1;
    const result = zmath_core_lgammaf(value, &sign_c);
    return .{ .value = result, .sign = @intCast(sign_c) };
}

/// Computes the correctly rounded binary64 lgamma result.
pub fn lgammad(value: f64) c23.LgammaResult(f64) {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    if (bits & magnitude_mask_d > infinity_bits_d) {
        return .{ .value = @bitCast(bits | quiet_nan_bit_d), .sign = 1 };
    }

    var sign_c: c_int = 1;
    const result = zmath_core_lgammad(value, &sign_c);
    return .{ .value = result, .sign = @intCast(sign_c) };
}
