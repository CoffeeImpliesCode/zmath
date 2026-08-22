const std = @import("std");
const scalar = @import("scalar.zig");

fn expectAtanhhBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.atanhh(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAtanhfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.atanhf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAtanhdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.atanhd(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "atanh exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.atanhh) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.atanhf) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.atanhd) == fn (f64) f64);
    }
}

test "atanh preserves signed zero and rounds representative values" {
    try expectAtanhhBits(0x0000, 0x0000);
    try expectAtanhhBits(0x8000, 0x8000);
    try expectAtanhhBits(0x3400, 0x3416);
    try expectAtanhhBits(0x3800, 0x3865);
    try expectAtanhhBits(0x3a00, 0x3bc9);

    try expectAtanhfBits(0x0000_0000, 0x0000_0000);
    try expectAtanhfBits(0x8000_0000, 0x8000_0000);
    try expectAtanhfBits(0x3e80_0000, 0x3e82_c578);
    try expectAtanhfBits(0x3f00_0000, 0x3f0c_9f54);
    try expectAtanhfBits(0x3f40_0000, 0x3f79_1395);

    try expectAtanhdBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectAtanhdBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectAtanhdBits(0x3fd0_0000_0000_0000, 0x3fd0_58ae_fa81_1452);
    try expectAtanhdBits(0x3fe0_0000_0000_0000, 0x3fe1_93ea_7aad_030b);
    try expectAtanhdBits(0x3fe8_0000_0000_0000, 0x3fef_2272_ae32_5a57);
}

test "atanh handles endpoints invalid inputs and NaNs" {
    try expectAtanhhBits(0x3c00, 0x7c00);
    try expectAtanhhBits(0xbc00, 0xfc00);
    try expectAtanhhBits(0x3c01, 0x7e00);
    try expectAtanhhBits(0x7c00, 0x7e00);
    try expectAtanhhBits(0x7c55, 0x7e55);
    try expectAtanhhBits(0xfc55, 0xfe55);

    try expectAtanhfBits(0x3f80_0000, 0x7f80_0000);
    try expectAtanhfBits(0xbf80_0000, 0xff80_0000);
    try expectAtanhfBits(0x3f80_0001, 0x7fc0_0000);
    try expectAtanhfBits(0x7f80_0000, 0x7fc0_0000);
    try expectAtanhfBits(0x7f81_2345, 0x7fc1_2345);
    try expectAtanhfBits(0xff81_2345, 0xffc1_2345);

    try expectAtanhdBits(0x3ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectAtanhdBits(0xbff0_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectAtanhdBits(0x3ff0_0000_0000_0001, 0x7ff8_0000_0000_0000);
    try expectAtanhdBits(0x7ff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectAtanhdBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectAtanhdBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "atanh covers tiny and CORE-MATH correction inputs" {
    try expectAtanhhBits(0x0001, 0x0001);
    try expectAtanhhBits(0x8001, 0x8001);
    try expectAtanhhBits(0x2d3c, 0x2d3f);

    try expectAtanhfBits(0x0000_0001, 0x0000_0001);
    try expectAtanhfBits(0x8000_0001, 0x8000_0001);

    try expectAtanhdBits(0x0000_0000_0000_0001, 0x0000_0000_0000_0001);
    try expectAtanhdBits(0x8000_0000_0000_0001, 0x8000_0000_0000_0001);
    try expectAtanhdBits(0x3fd2_dbb7_b1c9_1363, 0x3fd3_6f33_d51c_264d);
    try expectAtanhdBits(0xbfdc_493d_c899_e4a5, 0xbfde_611a_a58a_b608);
    try expectAtanhdBits(0x3ead_ffff_ffff_faba, 0x3eae_0000_0000_0384);
}
