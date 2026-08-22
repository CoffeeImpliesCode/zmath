const std = @import("std");
const scalar = @import("scalar.zig");

fn expectBits(input_bits: u64, expected_bits: u64) !void {
    const actual_bits: u64 = @bitCast(scalar.cosd(@as(f64, @bitCast(input_bits))));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "cosd exports a nearest-even binary64 function" {
    comptime std.debug.assert(@TypeOf(scalar.cosd) == fn (f64) f64);
}

test "cosd maps signed zero and tiny values to positive one" {
    try expectBits(0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectBits(0x8000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectBits(0x0000_0000_0000_0001, 0x3ff0_0000_0000_0000);
    try expectBits(0x800f_ffff_ffff_ffff, 0x3ff0_0000_0000_0000);
}

test "cosd quiets NaNs and maps infinities to canonical positive NaN" {
    try expectBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
    try expectBits(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectBits(0x7ff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
}

test "cosd reduces ordinary quadrant and huge finite arguments" {
    try expectBits(0x3ff0_0000_0000_0000, 0x3fe1_4a28_0fb5_068c);
    try expectBits(0xbff0_0000_0000_0000, 0x3fe1_4a28_0fb5_068c);
    try expectBits(0x3ff9_21fb_5444_2d18, 0x3c91_a626_3314_5c07);
    try expectBits(0x4009_21fb_5444_2d18, 0xbff0_0000_0000_0000);
    try expectBits(0x4019_21fb_5444_2d18, 0x3ff0_0000_0000_0000);
    try expectBits(0x7fe0_0000_0000_0000, 0xbfea_719f_26c2_32bf);
}

test "cosd handles accurate fallback hard cases" {
    try expectBits(0x3e88_0000_0000_0009, 0x3fef_ffff_ffff_ff70);
    try expectBits(0x3e98_0000_0000_0024, 0x3fef_ffff_ffff_fdc0);
    try expectBits(0x3ea8_0000_0000_0090, 0x3fef_ffff_ffff_f700);
    try expectBits(0x3eb2_0000_0000_00f3, 0x3fef_ffff_ffff_ebc0);
    try expectBits(0x3eb8_0000_0000_0240, 0x3fef_ffff_ffff_dc00);
}
