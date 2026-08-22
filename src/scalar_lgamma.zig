const std = @import("std");

/// The rounded logarithm of |gamma(x)| and the sign of gamma(x).
///
/// The sign is -1, 0, or 1. Callers must keep IEEE round-to-nearest,
/// ties-to-even active. The implementation does not inspect or change fenv.
pub fn LgammaResult(comptime T: type) type {
    return struct {
        value: T,
        sign: i2,
    };
}

extern fn zmath_core_lgammah(value: f16, sign: *c_int) f16;
extern fn zmath_core_lgammaf(value: f32, sign: *c_int) f32;
extern fn zmath_core_lgammad(value: f64, sign: *c_int) f64;

fn poleSignH(value: f16, sign: i2) i2 {
    const bits: u16 = @bitCast(value);
    if ((bits & 0x7fff) == 0) return if ((bits & 0x8000) != 0) -1 else 1;
    if (value < 0 and @trunc(value) == value) return -1;
    return sign;
}

fn poleSignF(value: f32, sign: i2) i2 {
    const bits: u32 = @bitCast(value);
    if ((bits & 0x7fff_ffff) == 0) return if ((bits & 0x8000_0000) != 0) -1 else 1;
    if (value < 0 and @trunc(value) == value) return -1;
    return sign;
}

fn poleSignD(value: f64, sign: i2) i2 {
    const bits: u64 = @bitCast(value);
    if ((bits & 0x7fff_ffff_ffff_ffff) == 0) {
        return if ((bits & 0x8000_0000_0000_0000) != 0) -1 else 1;
    }
    if (value < 0 and @trunc(value) == value) return -1;
    return sign;
}

/// Correctly rounded binary16 lgamma.
pub fn lgammah(value: f16) LgammaResult(f16) {
    @setFloatMode(.strict);
    var sign: c_int = 1;
    const result = zmath_core_lgammah(value, &sign);
    const sign_i2: i2 = @intCast(sign);
    return .{ .value = result, .sign = poleSignH(value, sign_i2) };
}

/// Correctly rounded binary32 lgamma.
pub fn lgammaf(value: f32) LgammaResult(f32) {
    @setFloatMode(.strict);
    var sign: c_int = 1;
    const result = zmath_core_lgammaf(value, &sign);
    const sign_i2: i2 = @intCast(sign);
    return .{ .value = result, .sign = poleSignF(value, sign_i2) };
}

/// Correctly rounded binary64 lgamma.
pub fn lgammad(value: f64) LgammaResult(f64) {
    @setFloatMode(.strict);
    var sign: c_int = 1;
    const result = zmath_core_lgammad(value, &sign);
    const sign_i2: i2 = @intCast(sign);
    return .{ .value = result, .sign = poleSignD(value, sign_i2) };
}

test "LgammaResult has the required shape" {
    comptime {
        if (@typeInfo(LgammaResult(f32)) != .@"struct") @compileError("not a struct");
        if (@TypeOf(@as(LgammaResult(f32), undefined).value) != f32) @compileError("bad value");
        if (@TypeOf(@as(LgammaResult(f32), undefined).sign) != i2) @compileError("bad sign");
    }
    try std.testing.expect(true);
}
