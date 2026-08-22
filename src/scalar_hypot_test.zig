const std = @import("std");
const scalar = @import("scalar.zig");

fn expectHypothBits(x_bits: u16, y_bits: u16, expected_bits: u16) !void {
    const x: f16 = @bitCast(x_bits);
    const y: f16 = @bitCast(y_bits);
    const actual_bits: u16 = @bitCast(scalar.hypoth(x, y));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectHypotfBits(x_bits: u32, y_bits: u32, expected_bits: u32) !void {
    const x: f32 = @bitCast(x_bits);
    const y: f32 = @bitCast(y_bits);
    const actual_bits: u32 = @bitCast(scalar.hypotf(x, y));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectHypotdBits(x_bits: u64, y_bits: u64, expected_bits: u64) !void {
    const x: f64 = @bitCast(x_bits);
    const y: f64 = @bitCast(y_bits);
    const actual_bits: u64 = @bitCast(scalar.hypotd(x, y));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "hypot exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.hypoth) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(scalar.hypotf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(scalar.hypotd) == fn (f64, f64) f64);
    }
}

test "hypot computes exact triples and ignores input signs" {
    try expectHypothBits(0x4200, 0x4400, 0x4500);
    try expectHypothBits(0xc200, 0x4400, 0x4500);
    try expectHypothBits(0x4400, 0xc200, 0x4500);

    try expectHypotfBits(0x4040_0000, 0x4080_0000, 0x40a0_0000);
    try expectHypotfBits(0xc040_0000, 0x4080_0000, 0x40a0_0000);
    try expectHypotfBits(0x4080_0000, 0xc040_0000, 0x40a0_0000);

    try expectHypotdBits(0x4008_0000_0000_0000, 0x4010_0000_0000_0000, 0x4014_0000_0000_0000);
    try expectHypotdBits(0xc008_0000_0000_0000, 0x4010_0000_0000_0000, 0x4014_0000_0000_0000);
    try expectHypotdBits(0x4010_0000_0000_0000, 0xc008_0000_0000_0000, 0x4014_0000_0000_0000);
}

test "hypot returns positive zero and the nonzero magnitude" {
    try expectHypothBits(0x0000, 0x8000, 0x0000);
    try expectHypothBits(0x8000, 0xbc00, 0x3c00);
    try expectHypotfBits(0x0000_0000, 0x8000_0000, 0x0000_0000);
    try expectHypotfBits(0x8000_0000, 0xbf80_0000, 0x3f80_0000);
    try expectHypotdBits(0x0000_0000_0000_0000, 0x8000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectHypotdBits(0x8000_0000_0000_0000, 0xbff0_0000_0000_0000, 0x3ff0_0000_0000_0000);
}

test "hypot handles subnormal rounding and finite range extremes" {
    try expectHypothBits(0x0001, 0x0001, 0x0001);
    try expectHypothBits(0x7bff, 0x0001, 0x7bff);
    try expectHypothBits(0x7bff, 0x7bff, 0x7c00);

    try expectHypotfBits(0x0000_0001, 0x0000_0001, 0x0000_0001);
    try expectHypotfBits(0x7f7f_ffff, 0x0000_0001, 0x7f7f_ffff);
    try expectHypotfBits(0x7f7f_ffff, 0x7f7f_ffff, 0x7f80_0000);
    try expectHypotfBits(0x7f7f_ffff, 0x797f_ffff, 0x7f7f_ffff);

    try expectHypotdBits(0x0000_0000_0000_0001, 0x0000_0000_0000_0001, 0x0000_0000_0000_0001);
    try expectHypotdBits(0x0010_0000_0000_0000, 0x0000_0000_0000_0001, 0x0010_0000_0000_0000);
    try expectHypotdBits(0x0010_0000_0000_0000, 0x0000_0000_0200_0000, 0x0010_0000_0000_0000);
    try expectHypotdBits(0x7fef_ffff_ffff_ffff, 0x0000_0000_0000_0001, 0x7fef_ffff_ffff_ffff);
    try expectHypotdBits(0x7fef_ffff_ffff_ffff, 0x7fef_ffff_ffff_ffff, 0x7ff0_0000_0000_0000);
}

test "hypot applies stable NaN selection and infinity precedence" {
    try expectHypothBits(0xfe55, 0x3c00, 0xfe55);
    try expectHypothBits(0x3c00, 0x7c55, 0x7e55);
    try expectHypothBits(0xfc00, 0x7e55, 0x7c00);
    try expectHypothBits(0x7c55, 0xfc00, 0x7e55);
    try expectHypothBits(0x7e55, 0x7c55, 0x7e55);

    try expectHypotfBits(0xffc1_2345, 0x3f80_0000, 0xffc1_2345);
    try expectHypotfBits(0x3f80_0000, 0x7f81_2345, 0x7fc1_2345);
    try expectHypotfBits(0xff80_0000, 0x7fc1_2345, 0x7f80_0000);
    try expectHypotfBits(0x7f81_2345, 0xff80_0000, 0x7fc1_2345);
    try expectHypotfBits(0x7fc1_2345, 0x7f81_2345, 0x7fc1_2345);

    try expectHypotdBits(0xfff8_1234_5678_9abc, 0x3ff0_0000_0000_0000, 0xfff8_1234_5678_9abc);
    try expectHypotdBits(0x3ff0_0000_0000_0000, 0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectHypotdBits(0xfff0_0000_0000_0000, 0x7ff8_1234_5678_9abc, 0x7ff0_0000_0000_0000);
    try expectHypotdBits(0x7ff0_1234_5678_9abc, 0xfff0_0000_0000_0000, 0x7ff8_1234_5678_9abc);
    try expectHypotdBits(0x7ff8_1234_5678_9abc, 0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
}

test "hypotf covers CORE-MATH refinement and scaling cases" {
    try expectHypotfBits(0x25dc_7285, 0x9fed_8f3c, 0x25dc_7285);
    try expectHypotfBits(0x2f01_daa6, 0x2936_535f, 0x2f01_daa6);
    try expectHypotfBits(0x7071_77fb, 0xf582_2659, 0x7582_265d);
    try expectHypotfBits(0x805a_c171, 0x00b5_98a3, 0x00cb_0303);
    try expectHypotfBits(0x0049_ac62, 0x8034_9e5e, 0x005a_88d1);
    try expectHypotfBits(0x611d_64c0, 0xe65d_2f61, 0x665d_2f64);
}

test "hypotd covers CORE-MATH subnormal and hard rounding cases" {
    try expectHypotdBits(0x000e_779e_1bd4_ad3f, 0x0006_d568_c1f6_fbab, 0x0010_0000_0000_0000);
    try expectHypotdBits(0x4330_0000_1cc7_84cb, 0x433f_dd67_bceb_35f0, 0x4341_d400_aba3_c4da);
    try expectHypotdBits(0x3d26_0000_0000_0000, 0x3b82_c2fc_5954_56a7, 0x3d26_0000_0000_0001);
    try expectHypotdBits(0x3f18_7de2_9ce1_0f34, 0x3ff0_0000_02d4_13cd, 0x3ff0_0000_03ff_ffff);
    try expectHypotdBits(0x0006_03e5_2daf_0bfd, 0x800a_622d_0a9a_433b, 0x000b_ffff_b1b0_6483);
    try expectHypotdBits(0x3ffa_308e_1455_f447, 0x3ff9_d931_a83e_f879, 0x4002_660d_009d_54f9);
}
