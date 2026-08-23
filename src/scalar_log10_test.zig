const std = @import("std");
const scalar = @import("scalar.zig");

fn expectLog10hBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.log10h(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectLog10fBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.log10f(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectLog10dBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.log10d(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "log10 exposes nearest-even H, F, and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.log10h) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.log10f) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.log10d) == fn (f64) f64);
    }
}

test "log10 handles exact special values" {
    try expectLog10hBits(0x0000, 0xfc00);
    try expectLog10hBits(0x8000, 0xfc00);
    try expectLog10hBits(0xbc00, 0x7e00);
    try expectLog10hBits(0xfc00, 0x7e00);
    try expectLog10hBits(0x7c00, 0x7c00);

    try expectLog10fBits(0x0000_0000, 0xff80_0000);
    try expectLog10fBits(0x8000_0000, 0xff80_0000);
    try expectLog10fBits(0xbf80_0000, 0x7fc0_0000);
    try expectLog10fBits(0xff80_0000, 0x7fc0_0000);
    try expectLog10fBits(0x7f80_0000, 0x7f80_0000);

    try expectLog10dBits(0x0000_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectLog10dBits(0x8000_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectLog10dBits(0xbff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectLog10dBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectLog10dBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
}

test "log10 quiets NaNs while preserving sign and payload" {
    try expectLog10hBits(0x7c55, 0x7e55);
    try expectLog10hBits(0xfc55, 0xfe55);
    try expectLog10hBits(0x7e55, 0x7e55);
    try expectLog10hBits(0xfe55, 0xfe55);

    try expectLog10fBits(0x7f81_2345, 0x7fc1_2345);
    try expectLog10fBits(0xff81_2345, 0xffc1_2345);
    try expectLog10fBits(0x7fc1_2345, 0x7fc1_2345);
    try expectLog10fBits(0xffc1_2345, 0xffc1_2345);

    try expectLog10dBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectLog10dBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
    try expectLog10dBits(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectLog10dBits(0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "log10 returns exact decimal powers and rounds log10 two" {
    try expectLog10hBits(0x3c00, 0x0000);
    try expectLog10hBits(0x4000, 0x34d1);
    try expectLog10hBits(0x4900, 0x3c00);
    try expectLog10hBits(0x5640, 0x4000);
    try expectLog10hBits(0x63d0, 0x4200);

    try expectLog10fBits(0x3f80_0000, 0x0000_0000);
    try expectLog10fBits(0x4000_0000, 0x3e9a_209b);
    try expectLog10fBits(0x4120_0000, 0x3f80_0000);
    try expectLog10fBits(0x42c8_0000, 0x4000_0000);

    try expectLog10dBits(0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectLog10dBits(0x4000_0000_0000_0000, 0x3fd3_4413_509f_79ff);
    try expectLog10dBits(0x4024_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectLog10dBits(0x4059_0000_0000_0000, 0x4000_0000_0000_0000);
}

test "log10 covers subnormal and normal boundaries" {
    try expectLog10hBits(0x0001, 0xc73a);
    try expectLog10hBits(0x03ff, 0xc437);
    try expectLog10hBits(0x0400, 0xc437);

    try expectLog10fBits(0x0000_0001, 0xc233_69f4);
    try expectLog10fBits(0x007f_ffff, 0xc217_b818);
    try expectLog10fBits(0x0080_0000, 0xc217_b818);

    try expectLog10dBits(0x0000_0000_0000_0001, 0xc074_34e6_420f_4374);
    try expectLog10dBits(0x000f_ffff_ffff_ffff, 0xc073_3a71_46f7_2a42);
    try expectLog10dBits(0x0010_0000_0000_0000, 0xc073_3a71_46f7_2a42);
}
