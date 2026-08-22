const std = @import("std");
const scalar = @import("scalar.zig");

fn expectTanhhBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.tanhh(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectTanhfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.tanhf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectTanhdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.tanhd(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "tanh exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.tanhh) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.tanhf) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.tanhd) == fn (f64) f64);
    }
}

test "tanh preserves signed zero and rounds representative values" {
    try expectTanhhBits(0x0000, 0x0000);
    try expectTanhhBits(0x8000, 0x8000);
    try expectTanhhBits(0x3800, 0x3765);
    try expectTanhhBits(0x3c00, 0x3a18);
    try expectTanhhBits(0x4000, 0x3bb6);

    try expectTanhfBits(0x0000_0000, 0x0000_0000);
    try expectTanhfBits(0x8000_0000, 0x8000_0000);
    try expectTanhfBits(0x3f00_0000, 0x3eec_9a9f);
    try expectTanhfBits(0x3f80_0000, 0x3f42_f7d6);
    try expectTanhfBits(0x4000_0000, 0x3f76_ca83);

    try expectTanhdBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectTanhdBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectTanhdBits(0x3fe0_0000_0000_0000, 0x3fdd_9353_d756_8af3);
    try expectTanhdBits(0x3ff0_0000_0000_0000, 0x3fe8_5efa_b514_f394);
    try expectTanhdBits(0x4000_0000_0000_0000, 0x3fee_d950_5e1b_c3d4);
}

test "tanh handles infinities and quiets NaNs with stable payloads" {
    try expectTanhhBits(0x7c00, 0x3c00);
    try expectTanhhBits(0xfc00, 0xbc00);
    try expectTanhhBits(0x7c55, 0x7e55);
    try expectTanhhBits(0xfc55, 0xfe55);

    try expectTanhfBits(0x7f80_0000, 0x3f80_0000);
    try expectTanhfBits(0xff80_0000, 0xbf80_0000);
    try expectTanhfBits(0x7f81_2345, 0x7fc1_2345);
    try expectTanhfBits(0xff81_2345, 0xffc1_2345);

    try expectTanhdBits(0x7ff0_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectTanhdBits(0xfff0_0000_0000_0000, 0xbff0_0000_0000_0000);
    try expectTanhdBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectTanhdBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "tanh covers CORE-MATH correction cases" {
    try expectTanhhBits(0xbfd5, 0xbbb0);

    try expectTanhdBits(0x3fca_c343_b179_fec4, 0x3fca_6124_99c5_3078);
    try expectTanhdBits(0x3fd0_0764_a988_bf73, 0x3fcf_6764_84c0_703b);
    try expectTanhdBits(0xbfd1_7d1e_8a63_711f, 0xbfd1_10e9_6a6c_2d96);
}

test "tanhh is odd for every finite binary16 input" {
    var magnitude: u16 = 0;
    while (magnitude < 0x7c00) : (magnitude += 1) {
        const positive_bits: u16 = @bitCast(scalar.tanhh(@bitCast(magnitude)));
        const negative_bits: u16 = @bitCast(scalar.tanhh(@bitCast(magnitude | 0x8000)));
        try std.testing.expectEqual(positive_bits | 0x8000, negative_bits);
    }
}
