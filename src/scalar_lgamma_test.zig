const std = @import("std");
const c23 = @import("scalar_c23_types.zig");
const scalar = @import("scalar.zig");

fn expectLgammaH(input_bits: u16, expected_bits: u16, expected_sign: i2) !void {
    const result = scalar.lgammah(@bitCast(input_bits));
    try std.testing.expectEqual(expected_bits, @as(u16, @bitCast(result.value)));
    try std.testing.expectEqual(expected_sign, result.sign);
}

fn expectLgammaF(input_bits: u32, expected_bits: u32, expected_sign: i2) !void {
    const result = scalar.lgammaf(@bitCast(input_bits));
    try std.testing.expectEqual(expected_bits, @as(u32, @bitCast(result.value)));
    try std.testing.expectEqual(expected_sign, result.sign);
}

fn expectLgammaD(input_bits: u64, expected_bits: u64, expected_sign: i2) !void {
    const result = scalar.lgammad(@bitCast(input_bits));
    try std.testing.expectEqual(expected_bits, @as(u64, @bitCast(result.value)));
    try std.testing.expectEqual(expected_sign, result.sign);
}

test "lgamma exposes the required result shape" {
    comptime {
        const result_type = c23.LgammaResult(f32);
        const initialized_result: result_type = .{ .value = 0, .sign = 0 };
        std.debug.assert(@typeInfo(result_type) == .@"struct");
        std.debug.assert(@TypeOf(initialized_result.value) == f32);
        std.debug.assert(@TypeOf(initialized_result.sign) == i2);
        std.debug.assert(@TypeOf(scalar.lgammah(@as(f16, 1.0))) == c23.LgammaResult(f16));
        std.debug.assert(@TypeOf(scalar.lgammaf(@as(f32, 1.0))) == c23.LgammaResult(f32));
        std.debug.assert(@TypeOf(scalar.lgammad(@as(f64, 1.0))) == c23.LgammaResult(f64));
    }
}

test "lgamma preserves signed-zero pole signs" {
    try expectLgammaH(0x0000, 0x7c00, 1);
    try expectLgammaH(0x8000, 0x7c00, -1);
    try expectLgammaF(0x0000_0000, 0x7f80_0000, 1);
    try expectLgammaF(0x8000_0000, 0x7f80_0000, -1);
    try expectLgammaD(0x0000_0000_0000_0000, 0x7ff0_0000_0000_0000, 1);
    try expectLgammaD(0x8000_0000_0000_0000, 0x7ff0_0000_0000_0000, -1);
}

test "lgamma returns exact finite positive and negative values" {
    try expectLgammaH(0x3800, 0x3894, 1);
    try expectLgammaH(0xb800, 0x3d10, -1);
    try expectLgammaH(0xbe00, 0x3ae1, 1);

    try expectLgammaF(0x3f00_0000, 0x3f12_8682, 1);
    try expectLgammaF(0xbf00_0000, 0x3fa1_fc4d, -1);
    try expectLgammaF(0xbfc0_0000, 0x3f5c_2c0b, 1);

    try expectLgammaD(0x3fe0_0000_0000_0000, 0x3fe2_50d0_48e7_a1bd, 1);
    try expectLgammaD(0xbfe0_0000_0000_0000, 0x3ff4_3f89_a3f0_edd6, -1);
    try expectLgammaD(0xbff8_0000_0000_0000, 0x3feb_8581_5182_0f86, 1);
}

test "lgamma returns positive infinity at negative integer poles" {
    try expectLgammaH(0xbc00, 0x7c00, -1);
    try expectLgammaH(0xc000, 0x7c00, -1);
    try expectLgammaH(0xca00, 0x7c00, -1);
    try expectLgammaF(0xbf80_0000, 0x7f80_0000, -1);
    try expectLgammaF(0xc000_0000, 0x7f80_0000, -1);
    try expectLgammaF(0xc140_0000, 0x7f80_0000, -1);
    try expectLgammaD(0xbff0_0000_0000_0000, 0x7ff0_0000_0000_0000, -1);
    try expectLgammaD(0xc000_0000_0000_0000, 0x7ff0_0000_0000_0000, -1);
    try expectLgammaD(0xc028_0000_0000_0000, 0x7ff0_0000_0000_0000, -1);
}

test "lgamma maps both infinities to positive infinity" {
    try expectLgammaH(0x7c00, 0x7c00, 1);
    try expectLgammaH(0xfc00, 0x7c00, 1);
    try expectLgammaF(0x7f80_0000, 0x7f80_0000, 1);
    try expectLgammaF(0xff80_0000, 0x7f80_0000, 1);
    try expectLgammaD(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000, 1);
    try expectLgammaD(0xfff0_0000_0000_0000, 0x7ff0_0000_0000_0000, 1);
}

test "lgamma quiets NaNs while preserving first-input sign and payload" {
    try expectLgammaH(0x7e55, 0x7e55, 1);
    try expectLgammaH(0xfe55, 0xfe55, 1);
    try expectLgammaH(0x7c55, 0x7e55, 1);
    try expectLgammaH(0xfc55, 0xfe55, 1);

    try expectLgammaF(0x7fc1_2345, 0x7fc1_2345, 1);
    try expectLgammaF(0xffc1_2345, 0xffc1_2345, 1);
    try expectLgammaF(0x7f81_2345, 0x7fc1_2345, 1);
    try expectLgammaF(0xff81_2345, 0xffc1_2345, 1);

    try expectLgammaD(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc, 1);
    try expectLgammaD(0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc, 1);
    try expectLgammaD(0x7ff1_2345_6789_abcd, 0x7ff9_2345_6789_abcd, 1);
    try expectLgammaD(0xfff1_2345_6789_abcd, 0xfff9_2345_6789_abcd, 1);
}
