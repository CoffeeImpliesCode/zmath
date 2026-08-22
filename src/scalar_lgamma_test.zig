const std = @import("std");
const zmath = @import("root.zig");
const scalar = @import("scalar.zig");

fn expectResult(
    comptime T: type,
    expected_bits: std.meta.Int(.unsigned, @bitSizeOf(T)),
    expected_sign: i2,
    actual: zmath.LgammaResult(T),
) !void {
    try std.testing.expectEqual(expected_bits, @as(@TypeOf(expected_bits), @bitCast(actual.value)));
    try std.testing.expectEqual(expected_sign, actual.sign);
}

test "lgammah returns pinned CORE-MATH value and sign" {
    try expectResult(f16, 0x7c00, 1, scalar.lgammah(@bitCast(@as(u16, 0x0000))));
    try expectResult(f16, 0x7c00, -1, scalar.lgammah(@bitCast(@as(u16, 0x8000))));
    try expectResult(f16, 0x3894, 1, scalar.lgammah(@as(f16, 0.5)));
    try expectResult(f16, 0x3d10, -1, scalar.lgammah(@as(f16, -0.5)));
    try expectResult(f16, 0x0000, 1, scalar.lgammah(@as(f16, 1.0)));
    try expectResult(f16, 0x3ae1, 1, scalar.lgammah(@as(f16, -1.5)));
    try expectResult(f16, 0x1c98, -1, scalar.lgammah(@as(f16, -12.0)));
    try expectResult(f16, 0xff23, 1, scalar.lgammah(@bitCast(@as(u16, 0xfd23))));
}

test "lgammaf returns pinned CORE-MATH value and sign" {
    try expectResult(f32, 0x7f80_0000, 1, scalar.lgammaf(@bitCast(@as(u32, 0x0000_0000))));
    try expectResult(f32, 0x7f80_0000, -1, scalar.lgammaf(@bitCast(@as(u32, 0x8000_0000))));
    try expectResult(f32, 0x3f12_8682, 1, scalar.lgammaf(@as(f32, 0.5)));
    try expectResult(f32, 0x3fa1_fc4d, -1, scalar.lgammaf(@as(f32, -0.5)));
    try expectResult(f32, 0x0000_0000, 1, scalar.lgammaf(@as(f32, 2.0)));
    try expectResult(f32, 0x3f5c_2c0b, 1, scalar.lgammaf(@as(f32, -1.5)));
    try expectResult(f32, 0x7fe1_2345, 1, scalar.lgammaf(@bitCast(@as(u32, 0x7fa1_2345))));
    try expectResult(f32, 0xffe1_2345, 1, scalar.lgammaf(@bitCast(@as(u32, 0xffa1_2345))));
}

test "lgammad returns pinned CORE-MATH value and sign" {
    try expectResult(f64, 0x7ff0_0000_0000_0000, 1, scalar.lgammad(@bitCast(@as(u64, 0))));
    try expectResult(f64, 0x7ff0_0000_0000_0000, -1, scalar.lgammad(@bitCast(@as(u64, 0x8000_0000_0000_0000))));
    try expectResult(f64, 0x3fe2_50d0_48e7_a1bd, 1, scalar.lgammad(@as(f64, 0.5)));
    try expectResult(f64, 0x3ff4_3f89_a3f0_edd6, -1, scalar.lgammad(@as(f64, -0.5)));
    try expectResult(f64, 0x0000_0000_0000_0000, 1, scalar.lgammad(@as(f64, 1.0)));
    try expectResult(f64, 0x3feb_8581_5182_0f86, 1, scalar.lgammad(@as(f64, -1.5)));
    try expectResult(f64, 0x7ff9_2345_6789_abcd, 1, scalar.lgammad(@bitCast(@as(u64, 0x7ff1_2345_6789_abcd))));
    try expectResult(f64, 0xfff9_2345_6789_abcd, 1, scalar.lgammad(@bitCast(@as(u64, 0xfff1_2345_6789_abcd))));
}

test "root lgamma dispatch preserves the result type" {
    const half = zmath.lgamma(@as(f16, 0.5));
    const single = zmath.lgamma(@as(f32, -0.5));
    const double = zmath.lgamma(@as(f64, 3.0));

    try std.testing.expect(@TypeOf(half) == zmath.LgammaResult(f16));
    try std.testing.expect(@TypeOf(single) == zmath.LgammaResult(f32));
    try std.testing.expect(@TypeOf(double) == zmath.LgammaResult(f64));
    try expectResult(f16, 0x3894, 1, half);
    try expectResult(f32, 0x3fa1_fc4d, -1, single);
    try expectResult(f64, 0x3fe6_2e42_fefa_39ef, 1, double);
}
