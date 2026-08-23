const std = @import("std");
const scalar = @import("scalar.zig");

fn expectExp10hBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u16, @bitCast(scalar.exp10h(input))));
}

fn expectExp10fBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u32, @bitCast(scalar.exp10f(input))));
}

fn expectExp10dBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u64, @bitCast(scalar.exp10d(input))));
}

test "exp10 exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.exp10h) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.exp10f) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.exp10d) == fn (f64) f64);
    }
}

test "exp10 preserves signed zeros as one" {
    try expectExp10hBits(0x0000, 0x3c00);
    try expectExp10hBits(0x8000, 0x3c00);
    try expectExp10fBits(0x0000_0000, 0x3f80_0000);
    try expectExp10fBits(0x8000_0000, 0x3f80_0000);
    try expectExp10dBits(0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectExp10dBits(0x8000_0000_0000_0000, 0x3ff0_0000_0000_0000);
}

test "exp10 maps infinities to exact limits" {
    try expectExp10hBits(0x7c00, 0x7c00);
    try expectExp10hBits(0xfc00, 0x0000);
    try expectExp10fBits(0x7f80_0000, 0x7f80_0000);
    try expectExp10fBits(0xff80_0000, 0x0000_0000);
    try expectExp10dBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectExp10dBits(0xfff0_0000_0000_0000, 0x0000_0000_0000_0000);
}

test "exp10 quiets NaNs while preserving sign and payload" {
    try expectExp10hBits(0x7c01, 0x7e01);
    try expectExp10hBits(0xfd55, 0xff55);
    try expectExp10fBits(0x7f81_2345, 0x7fc1_2345);
    try expectExp10fBits(0xff81_2345, 0xffc1_2345);
    try expectExp10dBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectExp10dBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "exp10 returns exact integral powers" {
    try expectExp10hBits(0x3c00, 0x4900);
    try expectExp10hBits(0xbc00, 0x2e66);
    try expectExp10fBits(0x3f80_0000, 0x4120_0000);
    try expectExp10fBits(0xbf80_0000, 0x3dcc_cccd);
    try expectExp10dBits(0x3ff0_0000_0000_0000, 0x4024_0000_0000_0000);
    try expectExp10dBits(0xbff0_0000_0000_0000, 0x3fb9_9999_9999_999a);
}

test "exp10d covers CORE-MATH database hard cases" {
    try expectExp10dBits(0x3f48_21e0_f2af_b970, 0x3ff0_06f3_a519_f98a);
    try expectExp10dBits(0x3fae_7b52_5705_edef, 0x3ff2_59cc_ebe3_17bf);
    try expectExp10dBits(0x4038_9063_309f_3004, 0x4508_3fe3_1fe6_a9d6);
    try expectExp10dBits(0xbf04_5ddb_1038_2e3f, 0x3fef_ff44_6cc3_068d);
    try expectExp10dBits(0xbfd5_6ff3_0582_2f26, 0x3fdd_9862_d491_3979);
    try expectExp10dBits(0xc05d_a5b1_0d86_89fd, 0x2750_a267_4040_67fb);
}

test "exp10h visits every finite encoding" {
    var raw_bits: u32 = 0;
    while (raw_bits < 0x1_0000) : (raw_bits += 1) {
        const input_bits: u16 = @intCast(raw_bits);
        if ((input_bits & 0x7fff) >= 0x7c00) continue;
        const output_bits: u16 = @bitCast(scalar.exp10h(@bitCast(input_bits)));
        try std.testing.expect((output_bits & 0x7c00) != 0x7c00 or (output_bits & 0x03ff) == 0);
        try std.testing.expect((output_bits & 0x8000) == 0);
    }
}
