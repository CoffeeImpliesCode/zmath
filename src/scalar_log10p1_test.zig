const std = @import("std");
const scalar = @import("scalar.zig");

fn expectLog10p1hBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u16, @bitCast(scalar.log10p1h(input))));
}

fn expectLog10p1fBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u32, @bitCast(scalar.log10p1f(input))));
}

fn expectLog10p1dBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u64, @bitCast(scalar.log10p1d(input))));
}

test "log10p1 exposes nearest-even H, F, and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.log10p1h) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.log10p1f) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.log10p1d) == fn (f64) f64);
    }
}

test "log10p1 preserves signed zero and handles domain limits" {
    try expectLog10p1hBits(0x0000, 0x0000);
    try expectLog10p1hBits(0x8000, 0x8000);
    try expectLog10p1hBits(0xbc00, 0xfc00);
    try expectLog10p1hBits(0xbc01, 0x7e00);
    try expectLog10p1hBits(0x7c00, 0x7c00);

    try expectLog10p1fBits(0x0000_0000, 0x0000_0000);
    try expectLog10p1fBits(0x8000_0000, 0x8000_0000);
    try expectLog10p1fBits(0xbf80_0000, 0xff80_0000);
    try expectLog10p1fBits(0xbf80_0001, 0x7fc0_0000);
    try expectLog10p1fBits(0x7f80_0000, 0x7f80_0000);

    try expectLog10p1dBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectLog10p1dBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectLog10p1dBits(0xbff0_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectLog10p1dBits(0xbff0_0000_0000_0001, 0x7ff8_0000_0000_0000);
    try expectLog10p1dBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
}

test "log10p1 quiets NaNs while preserving sign and payload" {
    try expectLog10p1hBits(0x7c55, 0x7e55);
    try expectLog10p1hBits(0xfc55, 0xfe55);
    try expectLog10p1fBits(0x7f81_2345, 0x7fc1_2345);
    try expectLog10p1fBits(0xff81_2345, 0xffc1_2345);
    try expectLog10p1dBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectLog10p1dBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "log10p1 returns exact decimal powers minus one" {
    try expectLog10p1hBits(0x4880, 0x3c00);
    try expectLog10p1hBits(0x5630, 0x4000);
    try expectLog10p1hBits(0x63ce, 0x4200);

    try expectLog10p1fBits(0x4110_0000, 0x3f80_0000);
    try expectLog10p1fBits(0x42c6_0000, 0x4000_0000);
    try expectLog10p1fBits(0x4479_c000, 0x4040_0000);

    try expectLog10p1dBits(0x4022_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectLog10p1dBits(0x4058_c000_0000_0000, 0x4000_0000_0000_0000);
    try expectLog10p1dBits(0x408f_3800_0000_0000, 0x4008_0000_0000_0000);
}

test "log10p1 rounds representative values" {
    try expectLog10p1hBits(0x3c00, 0x34d1);
    try expectLog10p1fBits(0x3f80_0000, 0x3e9a_209b);
    try expectLog10p1dBits(0x3ff0_0000_0000_0000, 0x3fd3_4413_509f_79ff);
}

test "log10p1 rounds pinned CORE-MATH hard cases" {
    try expectLog10p1hBits(0x9d08, 0x9861);
    try expectLog10p1hBits(0x89dc, 0x8517);
    try expectLog10p1hBits(0x213c, 0x1c86);
    try expectLog10p1hBits(0x0b80, 0x0683);

    try expectLog10p1fBits(0x399a_7c00, 0x3906_29e5);
    try expectLog10p1fBits(0xb051_e173, 0xafb6_4ccf);
    try expectLog10p1fBits(0xbb0e_a6d9, 0xba78_14e2);
    try expectLog10p1fBits(0xbd86_ffb9, 0xbcf2_9a9b);

    try expectLog10p1dBits(0xbf9e_0648_eff3_dad8, 0xbf8a_780e_0727_d64a);
    try expectLog10p1dBits(0xbf98_dd15_7e27_ade5, 0xbf85_dcde_5324_70b8);
    try expectLog10p1dBits(0xbf97_39b8_4607_1578, 0xbf84_67da_c57c_b711);
}
