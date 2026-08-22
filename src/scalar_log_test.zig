const std = @import("std");
const scalar = @import("scalar.zig");

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
