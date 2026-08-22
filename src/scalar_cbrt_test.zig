const std = @import("std");
const scalar = @import("scalar.zig");

fn expectCbrthBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.cbrth(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectCbrtfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.cbrtf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectCbrtdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.cbrtd(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "cbrt exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.cbrth) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.cbrtf) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.cbrtd) == fn (f64) f64);
    }
}

test "cbrt preserves signed zero and exact signed cubes" {
    try expectCbrthBits(0x0000, 0x0000);
    try expectCbrthBits(0x8000, 0x8000);
    try expectCbrthBits(0x3c00, 0x3c00);
    try expectCbrthBits(0x4800, 0x4000);
    try expectCbrthBits(0xcec0, 0xc200);

    try expectCbrtfBits(0x0000_0000, 0x0000_0000);
    try expectCbrtfBits(0x8000_0000, 0x8000_0000);
    try expectCbrtfBits(0x3f80_0000, 0x3f80_0000);
    try expectCbrtfBits(0x4100_0000, 0x4000_0000);
    try expectCbrtfBits(0xc1d8_0000, 0xc040_0000);

    try expectCbrtdBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectCbrtdBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectCbrtdBits(0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectCbrtdBits(0x4020_0000_0000_0000, 0x4000_0000_0000_0000);
    try expectCbrtdBits(0xc03b_0000_0000_0000, 0xc008_0000_0000_0000);
}

test "cbrt handles finite boundaries and rounds cbrt two" {
    try expectCbrthBits(0x0001, 0x1c00);
    try expectCbrthBits(0x03ff, 0x290a);
    try expectCbrthBits(0x0400, 0x290a);
    try expectCbrthBits(0x4000, 0x3d0a);
    try expectCbrthBits(0x7bff, 0x510a);

    try expectCbrtfBits(0x0000_0001, 0x26a1_4518);
    try expectCbrtfBits(0x007f_ffff, 0x2a7f_ffff);
    try expectCbrtfBits(0x0080_0000, 0x2a80_0000);
    try expectCbrtfBits(0x4000_0000, 0x3fa1_4518);
    try expectCbrtfBits(0x7f7f_ffff, 0x54cb_2ff5);

    try expectCbrtdBits(0x0000_0000_0000_0001, 0x2990_0000_0000_0000);
    try expectCbrtdBits(0x000f_ffff_ffff_ffff, 0x2aa4_28a2_f98d_728a);
    try expectCbrtdBits(0x0010_0000_0000_0000, 0x2aa4_28a2_f98d_728b);
    try expectCbrtdBits(0x4000_0000_0000_0000, 0x3ff4_28a2_f98d_728b);
    try expectCbrtdBits(0x7fef_ffff_ffff_ffff, 0x5544_28a2_f98d_728b);
}

test "cbrt preserves infinities and quiets NaNs with stable payloads" {
    try expectCbrthBits(0x7c00, 0x7c00);
    try expectCbrthBits(0xfc00, 0xfc00);
    try expectCbrthBits(0x7c55, 0x7e55);
    try expectCbrthBits(0xfc55, 0xfe55);

    try expectCbrtfBits(0x7f80_0000, 0x7f80_0000);
    try expectCbrtfBits(0xff80_0000, 0xff80_0000);
    try expectCbrtfBits(0x7f81_2345, 0x7fc1_2345);
    try expectCbrtfBits(0xff81_2345, 0xffc1_2345);

    try expectCbrtdBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectCbrtdBits(0xfff0_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectCbrtdBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectCbrtdBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "cbrtd covers CORE-MATH nearest-even hard cases" {
    try expectCbrtdBits(0x4009_b782_23aa_307c, 0x3ff7_9d15_d0e8_d59c);
    try expectCbrtdBits(0xc009_b782_23aa_307c, 0xbff7_9d15_d0e8_d59c);
    try expectCbrtdBits(0x401a_202b_fc89_ddff, 0x3ffd_e87a_a837_820f);
    try expectCbrtdBits(0xc01a_202b_fc89_ddff, 0xbffd_e87a_a837_820f);
}

test "cbrth is odd for every finite binary16 input" {
    var magnitude: u16 = 0;
    while (magnitude < 0x7c00) : (magnitude += 1) {
        const positive_bits: u16 = @bitCast(scalar.cbrth(@bitCast(magnitude)));
        const negative_input = magnitude | 0x8000;
        const negative_bits: u16 = @bitCast(scalar.cbrth(@bitCast(negative_input)));
        try std.testing.expectEqual(positive_bits | 0x8000, negative_bits);
    }
}
