const std = @import("std");
const ieee = @import("ieee.zig");
const scalar = @import("scalar.zig");

fn expectExp10m1hBits(input_bits: u16, expected_bits: u16) !void {
    const input = ieee.fromBitsh(input_bits);
    try std.testing.expectEqual(expected_bits, ieee.toBitsh(scalar.exp10m1h(input)));
}

fn expectExp10m1fBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u32, @bitCast(scalar.exp10m1f(input))));
}

fn expectExp10m1dBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u64, @bitCast(scalar.exp10m1d(input))));
}

test "exp10m1 exports use nearest-even scalar signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.exp10m1h) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.exp10m1f) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.exp10m1d) == fn (f64) f64);
    }
}

test "exp10m1 preserves signed zeros" {
    try expectExp10m1hBits(0x0000, 0x0000);
    try expectExp10m1hBits(0x8000, 0x8000);
    try expectExp10m1fBits(0x0000_0000, 0x0000_0000);
    try expectExp10m1fBits(0x8000_0000, 0x8000_0000);
    try expectExp10m1dBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectExp10m1dBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
}

test "exp10m1 maps infinities to exact limits" {
    try expectExp10m1hBits(0x7c00, 0x7c00);
    try expectExp10m1hBits(0xfc00, 0xbc00);
    try expectExp10m1fBits(0x7f80_0000, 0x7f80_0000);
    try expectExp10m1fBits(0xff80_0000, 0xbf80_0000);
    try expectExp10m1dBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectExp10m1dBits(0xfff0_0000_0000_0000, 0xbff0_0000_0000_0000);
}

test "exp10m1 quiets NaNs while preserving sign and payload" {
    try expectExp10m1hBits(0x7c01, 0x7e01);
    try expectExp10m1hBits(0xfd55, 0xff55);
    try expectExp10m1fBits(0x7f81_2345, 0x7fc1_2345);
    try expectExp10m1fBits(0xff81_2345, 0xffc1_2345);
    try expectExp10m1dBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectExp10m1dBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "exp10m1 returns exact integral powers" {
    try expectExp10m1hBits(0x3c00, 0x4880);
    try expectExp10m1hBits(0x4000, 0x5630);
    try expectExp10m1fBits(0x3f80_0000, 0x4110_0000);
    try expectExp10m1fBits(0x4000_0000, 0x42c6_0000);
    try expectExp10m1fBits(0x4040_0000, 0x4479_c000);
    try expectExp10m1fBits(0x4080_0000, 0x461c_3c00);
    try expectExp10m1fBits(0x40a0_0000, 0x47c3_4f80);
    try expectExp10m1fBits(0x40c0_0000, 0x4974_23f0);
    try expectExp10m1fBits(0x40e0_0000, 0x4b18_967f);
    try expectExp10m1dBits(0x3ff0_0000_0000_0000, 0x4022_0000_0000_0000);
    try expectExp10m1dBits(0x4000_0000_0000_0000, 0x4058_c000_0000_0000);
    try expectExp10m1dBits(0x4008_0000_0000_0000, 0x408f_3800_0000_0000);
    try expectExp10m1dBits(0x4010_0000_0000_0000, 0x40c3_8780_0000_0000);
    try expectExp10m1dBits(0x4014_0000_0000_0000, 0x40f8_69f0_0000_0000);
    try expectExp10m1dBits(0x4018_0000_0000_0000, 0x412e_847e_0000_0000);
    try expectExp10m1dBits(0x401c_0000_0000_0000, 0x4163_12cf_e000_0000);
    try expectExp10m1dBits(0x4020_0000_0000_0000, 0x4197_d783_fc00_0000);
    try expectExp10m1dBits(0x4022_0000_0000_0000, 0x41cd_cd64_ff80_0000);
    try expectExp10m1dBits(0x4024_0000_0000_0000, 0x4202_a05f_1ff8_0000);
    try expectExp10m1dBits(0x4026_0000_0000_0000, 0x4237_4876_e7ff_0000);
    try expectExp10m1dBits(0x4028_0000_0000_0000, 0x426d_1a94_a1ff_e000);
    try expectExp10m1dBits(0x402a_0000_0000_0000, 0x42a2_309c_e53f_fe00);
    try expectExp10m1dBits(0x402c_0000_0000_0000, 0x42d6_bcc4_1e8f_ffc0);
    try expectExp10m1dBits(0x402e_0000_0000_0000, 0x430c_6bf5_2633_fff8);
}

test "exp10m1h rounds all pinned exceptional inputs" {
    try expectExp10m1hBits(0x8757, 0x8c39);
    try expectExp10m1hBits(0x9788, 0x9c53);
    try expectExp10m1hBits(0xa1fe, 0xa6ce);
    try expectExp10m1hBits(0x2147, 0x2626);
    try expectExp10m1hBits(0x2394, 0x2870);
    try expectExp10m1hBits(0x23fc, 0x28ae);
    try expectExp10m1hBits(0xa4a8, 0xa940);
    try expectExp10m1hBits(0xa4f6, 0xa996);
    try expectExp10m1hBits(0x286c, 0x2d4c);
    try expectExp10m1hBits(0x24e6, 0x29c4);
    try expectExp10m1hBits(0x2ce8, 0x322d);
    try expectExp10m1hBits(0x2d71, 0x32eb);
    try expectExp10m1hBits(0x2615, 0x2b32);
    try expectExp10m1hBits(0x2633, 0x2b57);
    try expectExp10m1hBits(0x2ffa, 0x3551);
}

test "exp10m1f rounds pinned exceptional inputs" {
    try expectExp10m1fBits(0x2c99_4b7b, 0x2d30_7cba);
    try expectExp10m1fBits(0xb6fa_215b, 0xb78f_fc3f);
}

test "exp10m1d applies pinned correction triples" {
    try expectExp10m1dBits(0xbfa7_dcbb_bb8f_ad9c, 0xbfba_0cb7_77fb_6deb);
    try expectExp10m1dBits(0xc024_657d_88b3_9382, 0xbfef_ffff_fff7_4aeb);
    try expectExp10m1dBits(0xbc7f_b30f_6a48_e36d, 0xbc92_3f67_77b3_aab6);
}

test "exp10m1d applies pinned subnormal corrections" {
    try expectExp10m1dBits(0x0000_86c7_3059_343c, 0x0001_3656_8740_cb55);
    try expectExp10m1dBits(0x8000_86c7_3059_343c, 0x8001_3656_8740_cb55);
    try expectExp10m1dBits(0x0001_3a7b_70d0_248c, 0x0002_d41f_3b97_2fc7);
    try expectExp10m1dBits(0x8001_3a7b_70d0_248c, 0x8002_d41f_3b97_2fc7);
}
