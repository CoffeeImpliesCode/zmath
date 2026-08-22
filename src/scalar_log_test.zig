const std = @import("std");
const scalar = @import("scalar.zig");

fn expectLoghBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.logh(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "logh returns exact special values" {
    try expectLoghBits(0x0000, 0xfc00);
    try expectLoghBits(0x8000, 0xfc00);
    try expectLoghBits(0x3c00, 0x0000);
    try expectLoghBits(0x7c00, 0x7c00);
    try expectLoghBits(0xbc00, 0x7e00);
    try expectLoghBits(0xfc00, 0x7e00);
}

test "logh quiets NaNs while preserving sign and payload" {
    try expectLoghBits(0x7e55, 0x7e55);
    try expectLoghBits(0x7c55, 0x7e55);
    try expectLoghBits(0xfe55, 0xfe55);
    try expectLoghBits(0xfc55, 0xfe55);
}

fn expectLogdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.logd(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "logd has the nearest-even binary64 signature" {
    comptime {
        std.debug.assert(@TypeOf(scalar.logd) == fn (f64) f64);
    }
}

test "logd returns exact special values" {
    try expectLogdBits(0x0000_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectLogdBits(0x8000_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectLogdBits(0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectLogdBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectLogdBits(0xbff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectLogdBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
}

test "logd quiets NaNs while preserving sign and payload" {
    try expectLogdBits(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectLogdBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectLogdBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "logd handles normal and subnormal boundaries" {
    try expectLogdBits(0x0000_0000_0000_0001, 0xc087_4385_446d_71c3);
    try expectLogdBits(0x0010_0000_0000_0000, 0xc086_232b_dd7a_bcd2);
    try expectLogdBits(0x3fe0_0000_0000_0000, 0xbfe6_2e42_fefa_39ef);
    try expectLogdBits(0x4000_0000_0000_0000, 0x3fe6_2e42_fefa_39ef);
    try expectLogdBits(0x7fef_ffff_ffff_ffff, 0x4086_2e42_fefa_39ef);
}

test "logd rounds values adjacent to one" {
    try expectLogdBits(0x3fef_ffff_ffff_ffff, 0xbca0_0000_0000_0000);
    try expectLogdBits(0x3ff0_0000_0000_0001, 0x3caf_ffff_ffff_ffff);
}

test "logd resolves CORE-MATH hard cases" {
    try expectLogdBits(0x3ffa_6ae5_1423_26b5, 0x3fe0_0bcc_31eb_ded7);
    try expectLogdBits(0x3ffa_6ba9_5b15_20f9, 0x3fe0_0cb9_ef61_45e1);
    try expectLogdBits(0x400f_fcd1_f1e0_c16f, 0x3ff6_2cab_e3b0_07f3);
    try expectLogdBits(0x47c7_1f7c_59ed_e8e0, 0x4055_c0be_a603_6ab8);
}

fn expectLogfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.logf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "logf uses the nearest-even binary32 API" {
    comptime std.debug.assert(@TypeOf(scalar.logf) == fn (f32) f32);
}

test "logf maps exact special values" {
    try expectLogfBits(0x0000_0000, 0xff80_0000);
    try expectLogfBits(0x8000_0000, 0xff80_0000);
    try expectLogfBits(0x7f80_0000, 0x7f80_0000);
    try expectLogfBits(0xff80_0000, 0x7fc0_0000);
    try expectLogfBits(0x7fc1_2345, 0x7fc1_2345);
    try expectLogfBits(0x7fa1_2345, 0x7fe1_2345);
    try expectLogfBits(0xffc1_2345, 0xffc1_2345);
    try expectLogfBits(0xffa1_2345, 0xffe1_2345);
    try expectLogfBits(0xbf80_0000, 0x7fc0_0000);
}

test "logf rounds subnormal and normal reductions exactly" {
    const cases = [_]struct { input_bits: u32, expected_bits: u32 }{
        .{ .input_bits = 0x0000_0001, .expected_bits = 0xc2ce_8ed0 },
        .{ .input_bits = 0x0000_0002, .expected_bits = 0xc2cd_2bec },
        .{ .input_bits = 0x0000_0003, .expected_bits = 0xc2cc_5c53 },
        .{ .input_bits = 0x007f_ffff, .expected_bits = 0xc2ae_ac50 },
        .{ .input_bits = 0x0080_0000, .expected_bits = 0xc2ae_ac50 },
        .{ .input_bits = 0x0100_0000, .expected_bits = 0xc2ad_496b },
        .{ .input_bits = 0x3a80_0000, .expected_bits = 0xc0dd_ce9e },
        .{ .input_bits = 0x3d_cccc_cd, .expected_bits = 0xc013_5d8e },
        .{ .input_bits = 0x3e_aaaa_ab, .expected_bits = 0xbf8c_9f54 },
        .{ .input_bits = 0x3f00_0000, .expected_bits = 0xbf31_7218 },
        .{ .input_bits = 0x3f35_04f3, .expected_bits = 0xbeb1_7219 },
        .{ .input_bits = 0x3f7f_ffff, .expected_bits = 0xb380_0000 },
        .{ .input_bits = 0x3f80_0000, .expected_bits = 0x0000_0000 },
        .{ .input_bits = 0x3f80_0001, .expected_bits = 0x33ff_ffff },
        .{ .input_bits = 0x3f80_0002, .expected_bits = 0x347f_fffe },
        .{ .input_bits = 0x3f80_0100, .expected_bits = 0x37ff_ff00 },
        .{ .input_bits = 0x3f90_0000, .expected_bits = 0x3df1_383b },
        .{ .input_bits = 0x3fb5_04f3, .expected_bits = 0x3eb1_7217 },
        .{ .input_bits = 0x3fc0_0000, .expected_bits = 0x3ecf_991f },
        .{ .input_bits = 0x3fe0_0000, .expected_bits = 0x3f0f_42fb },
        .{ .input_bits = 0x3fff_fffe, .expected_bits = 0x3f31_7216 },
        .{ .input_bits = 0x3fff_ffff, .expected_bits = 0x3f31_7217 },
        .{ .input_bits = 0x4000_0000, .expected_bits = 0x3f31_7218 },
        .{ .input_bits = 0x4000_0001, .expected_bits = 0x3f31_721a },
        .{ .input_bits = 0x4040_0000, .expected_bits = 0x3f8c_9f54 },
        .{ .input_bits = 0x4120_0000, .expected_bits = 0x4013_5d8e },
        .{ .input_bits = 0x42c8_0000, .expected_bits = 0x4093_5d8e },
        .{ .input_bits = 0x7f7f_ffff, .expected_bits = 0x42b1_7218 },
    };

    for (cases) |case| {
        try expectLogfBits(case.input_bits, case.expected_bits);
    }
}

test "logf rounds accurate fallback and low-correction vectors" {
    const cases = [_]struct { input_bits: u32, expected_bits: u32 }{
        // This input is outside the near-one accurate-path shortcut.
        .{ .input_bits = 0x3f80_2000, .expected_bits = 0x3a7f_e005 },
        // These inputs exercise the final low-28-bit rounding correction.
        .{ .input_bits = 0x1f11_6ab8, .expected_bits = 0xc233_b53e },
        .{ .input_bits = 0x3c41_3d3a, .expected_bits = 0xc08e_158f },
        .{ .input_bits = 0x4117_8feb, .expected_bits = 0x400f_e5e7 },
        .{ .input_bits = 0x4c5d_65a5, .expected_bits = 0x418f_034b },
        .{ .input_bits = 0x4d60_4ebe, .expected_bits = 0x419a_352c },
        .{ .input_bits = 0x65d8_90d3, .expected_bits = 0x4254_d1f9 },
        .{ .input_bits = 0x66a8_c860, .expected_bits = 0x4259_5e46 },
        .{ .input_bits = 0x6f31_a8ec, .expected_bits = 0x4284_5a89 },
    };

    for (cases) |case| {
        try expectLogfBits(case.input_bits, case.expected_bits);
    }
}
