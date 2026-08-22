const std = @import("std");
const scalar = @import("scalar.zig");

fn expectRsqrthBits(input: u16, expected: u16) !void {
    const value: f16 = @bitCast(input);
    try std.testing.expectEqual(expected, @as(u16, @bitCast(scalar.rsqrth(value))));
}

fn expectRsqrtfBits(input: u32, expected: u32) !void {
    const value: f32 = @bitCast(input);
    try std.testing.expectEqual(expected, @as(u32, @bitCast(scalar.rsqrtf(value))));
}

fn expectRsqrtdBits(input: u64, expected: u64) !void {
    const value: f64 = @bitCast(input);
    try std.testing.expectEqual(expected, @as(u64, @bitCast(scalar.rsqrtd(value))));
}

test "Rsqrt returns exact ordinary values" {
    try expectRsqrthBits(0x3c00, 0x3c00);
    try expectRsqrthBits(0x4400, 0x3800);
    try expectRsqrthBits(0x4000, 0x39a8);

    try expectRsqrtfBits(0x3f80_0000, 0x3f80_0000);
    try expectRsqrtfBits(0x4080_0000, 0x3f00_0000);
    try expectRsqrtfBits(0x4000_0000, 0x3f35_04f3);

    try expectRsqrtdBits(0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectRsqrtdBits(0x4010_0000_0000_0000, 0x3fe0_0000_0000_0000);
    try expectRsqrtdBits(0x4000_0000_0000_0000, 0x3fe6_a09e_667f_3bcd);
}

test "Rsqrt covers pinned CORE-MATH hard cases" {
    try expectRsqrtfBits(0x002f_7e2a, 0x5f52_22e0);
    try expectRsqrtfBits(0x0055_b7bd, 0x5f1c_6a5a);
    try expectRsqrtfBits(0x00bd_f8a8, 0x5ed2_22e0);

    try expectRsqrtdBits(0x3ffa_6a9c_c15a_bcce, 0x3fe8_e77a_118a_3095);
    try expectRsqrtdBits(0x400c_562b_8574_53dd, 0x3fe1_00b9_26df_6e73);
    try expectRsqrtdBits(0x000f_ffff_ffff_ffff, 0x5fe0_0000_0000_0001);
}

test "Rsqrt handles subnormal inputs" {
    try expectRsqrthBits(0x0001, 0x6c00);
    try expectRsqrtfBits(0x0000_0001, 0x64b5_04f3);
    try expectRsqrtdBits(0x0000_0000_0000_0001, 0x6180_0000_0000_0000);
}

test "Rsqrt handles signed zero and infinity" {
    try expectRsqrthBits(0x0000, 0x7c00);
    try expectRsqrthBits(0x8000, 0xfc00);
    try expectRsqrthBits(0x7c00, 0x0000);

    try expectRsqrtfBits(0x0000_0000, 0x7f80_0000);
    try expectRsqrtfBits(0x8000_0000, 0xff80_0000);
    try expectRsqrtfBits(0x7f80_0000, 0x0000_0000);

    try expectRsqrtdBits(0x0000_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectRsqrtdBits(0x8000_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectRsqrtdBits(0x7ff0_0000_0000_0000, 0x0000_0000_0000_0000);
}

test "Rsqrt returns stable quiet NaNs" {
    try expectRsqrthBits(0xbc00, 0x7e00);
    try expectRsqrthBits(0xfc00, 0x7e00);
    try expectRsqrthBits(0x7d55, 0x7f55);
    try expectRsqrthBits(0xfd55, 0xff55);

    try expectRsqrtfBits(0xbf80_0000, 0x7fc0_0000);
    try expectRsqrtfBits(0xff80_0000, 0x7fc0_0000);
    try expectRsqrtfBits(0x7fa1_2345, 0x7fe1_2345);
    try expectRsqrtfBits(0xffa1_2345, 0xffe1_2345);

    try expectRsqrtdBits(0xbff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectRsqrtdBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectRsqrtdBits(0x7ff1_2345_6789_abcd, 0x7ff9_2345_6789_abcd);
    try expectRsqrtdBits(0xfff1_2345_6789_abcd, 0xfff9_2345_6789_abcd);
}
