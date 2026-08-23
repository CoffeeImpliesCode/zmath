const std = @import("std");
const ieee = @import("ieee.zig");
const scalar = @import("scalar.zig");

fn expectExp2m1hBits(input_bits: u16, expected_bits: u16) !void {
    const input = ieee.fromBitsh(input_bits);
    try std.testing.expectEqual(expected_bits, ieee.toBitsh(scalar.exp2m1h(input)));
}

fn expectExp2m1fBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u32, @bitCast(scalar.exp2m1f(input))));
}

fn expectExp2m1dBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u64, @bitCast(scalar.exp2m1d(input))));
}

test "exp2m1 exports nearest-even scalar signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.exp2m1h) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.exp2m1f) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.exp2m1d) == fn (f64) f64);
    }
}

test "exp2m1 preserves signed zeros" {
    try expectExp2m1hBits(0x0000, 0x0000);
    try expectExp2m1hBits(0x8000, 0x8000);
    try expectExp2m1fBits(0x0000_0000, 0x0000_0000);
    try expectExp2m1fBits(0x8000_0000, 0x8000_0000);
    try expectExp2m1dBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectExp2m1dBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
}

test "exp2m1 maps infinities to exact limits" {
    try expectExp2m1hBits(0x7c00, 0x7c00);
    try expectExp2m1hBits(0xfc00, 0xbc00);
    try expectExp2m1fBits(0x7f80_0000, 0x7f80_0000);
    try expectExp2m1fBits(0xff80_0000, 0xbf80_0000);
    try expectExp2m1dBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectExp2m1dBits(0xfff0_0000_0000_0000, 0xbff0_0000_0000_0000);
}

test "exp2m1 quiets NaNs while preserving sign and payload" {
    try expectExp2m1hBits(0x7c01, 0x7e01);
    try expectExp2m1hBits(0xfd55, 0xff55);
    try expectExp2m1fBits(0x7f81_2345, 0x7fc1_2345);
    try expectExp2m1fBits(0xff81_2345, 0xffc1_2345);
    try expectExp2m1dBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectExp2m1dBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "exp2m1 returns exact integral powers" {
    try expectExp2m1hBits(0xbc00, 0xb800);
    try expectExp2m1hBits(0x3c00, 0x3c00);
    try expectExp2m1hBits(0x4000, 0x4200);
    try expectExp2m1fBits(0xbf80_0000, 0xbf00_0000);
    try expectExp2m1fBits(0x3f80_0000, 0x3f80_0000);
    try expectExp2m1fBits(0x4000_0000, 0x4040_0000);
    try expectExp2m1dBits(0xbff0_0000_0000_0000, 0xbfe0_0000_0000_0000);
    try expectExp2m1dBits(0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectExp2m1dBits(0x4000_0000_0000_0000, 0x4008_0000_0000_0000);
}

test "exp2m1h rounds pinned exceptional inputs" {
    try expectExp2m1hBits(0x118c, 0x0fb1);
    try expectExp2m1hBits(0x21bc, 0x1ffb);
    try expectExp2m1hBits(0x8da9, 0x8bd9);
    try expectExp2m1hBits(0x973f, 0x9505);
    try expectExp2m1hBits(0x8461, 0x8309);
}

test "exp2m1f rounds pinned exceptional inputs" {
    try expectExp2m1fBits(0xb3d8_5005, 0xb395_efbb);
    try expectExp2m1fBits(0x3338_428d, 0x32ff_7045);
    try expectExp2m1fBits(0x388b_ca4f, 0x3841_cb81);
}

test "exp2m1d applies each pinned correction path" {
    try expectExp2m1dBits(0xbc7a_1682_6a8e_825d, 0xbc72_1530_a306_cc85);
    try expectExp2m1dBits(0xbfbf_6ec7_3d39_48c3, 0xbfb4_e2d8_b0ce_ad45);
    try expectExp2m1dBits(0xbffd_a226_1125_3866, 0xbfe7_22e3_006b_8e0b);
}
