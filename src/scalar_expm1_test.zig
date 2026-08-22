const std = @import("std");
const ieee = @import("ieee.zig");
const scalar = @import("scalar.zig");

fn expectExpm1hBits(input_bits: u16, expected_bits: u16) !void {
    const input = ieee.fromBitsh(input_bits);
    try std.testing.expectEqual(expected_bits, ieee.toBitsh(scalar.expm1h(input)));
}

fn expectExpm1fBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u32, @bitCast(scalar.expm1f(input))));
}

fn expectExpm1dBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u64, @bitCast(scalar.expm1d(input))));
}

test "expm1 exports use nearest-even scalar signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.expm1h) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.expm1f) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.expm1d) == fn (f64) f64);
    }
}

test "expm1 preserves signed zeros" {
    try expectExpm1hBits(0x0000, 0x0000);
    try expectExpm1hBits(0x8000, 0x8000);
    try expectExpm1fBits(0x0000_0000, 0x0000_0000);
    try expectExpm1fBits(0x8000_0000, 0x8000_0000);
    try expectExpm1dBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectExpm1dBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
}

test "expm1 maps infinities to exact limits" {
    try expectExpm1hBits(0x7c00, 0x7c00);
    try expectExpm1hBits(0xfc00, 0xbc00);
    try expectExpm1fBits(0x7f80_0000, 0x7f80_0000);
    try expectExpm1fBits(0xff80_0000, 0xbf80_0000);
    try expectExpm1dBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectExpm1dBits(0xfff0_0000_0000_0000, 0xbff0_0000_0000_0000);
}

test "expm1 quiets NaNs while preserving sign and payload" {
    try expectExpm1hBits(0x7e55, 0x7e55);
    try expectExpm1hBits(0xfd55, 0xff55);
    try expectExpm1fBits(0x7fc1_2345, 0x7fc1_2345);
    try expectExpm1fBits(0xff81_2345, 0xffc1_2345);
    try expectExpm1dBits(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectExpm1dBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "expm1 rounds representative CORE-MATH values" {
    try expectExpm1hBits(0x11a8, 0x11a9);
    try expectExpm1hBits(0x3c00, 0x3ee0);
    try expectExpm1hBits(0xbc00, 0xb90f);
    try expectExpm1hBits(0x498b, 0x7bf7);
    try expectExpm1hBits(0x498c, 0x7c00);

    try expectExpm1fBits(0x3f80_0000, 0x3fdb_f0a9);
    try expectExpm1fBits(0xbf80_0000, 0xbf21_d2a7);
    try expectExpm1fBits(0x3e20_0000, 0x3e2d_2d63);
    try expectExpm1fBits(0xbe20_0000, 0xbe14_205b);
    try expectExpm1fBits(0x42b1_7217, 0x7f7f_ff84);
    try expectExpm1fBits(0x42b1_7218, 0x7f80_0000);
    try expectExpm1fBits(0xc18a_a123, 0xbf80_0000);

    try expectExpm1dBits(0x3ff0_0000_0000_0000, 0x3ffb_7e15_1628_aed3);
    try expectExpm1dBits(0xbff0_0000_0000_0000, 0xbfe4_3a54_e4e9_8864);
    try expectExpm1dBits(0x3fd0_0000_0000_0000, 0x3fd2_2d78_f0fa_061a);
    try expectExpm1dBits(0xbfd0_0000_0000_0000, 0xbfcc_5041_854d_f7d4);
    try expectExpm1dBits(0x4086_2e42_fefa_39ef, 0x7fef_ffff_ffff_ff2a);
}

test "expm1d rounds accurate-small hard cases" {
    try expectExpm1dBits(0x3f30_0288_1d0c_becf, 0x3f30_0308_483d_bda3);
    try expectExpm1dBits(0x3d80_00bf_fb80_30a7, 0x3d80_00bf_fb80_40a8);
    try expectExpm1dBits(0xbe30_b5d6_cc46_b3f8, 0xbe30_b5d6_cbbb_15f7);
}

test "expm1d applies the database correction" {
    try expectExpm1dBits(0x3fd8_172a_0e02_f90e, 0x3fdd_404e_9760_1d65);
    try expectExpm1dBits(0xbff0_a54d_8778_3d6f, 0xbfe4_b188_7d4d_477c);
}
