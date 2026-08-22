const std = @import("std");
const scalar = @import("scalar.zig");

fn expectBits(input_bits: u64, expected_bits: u64) !void {
    const actual_bits: u64 = @bitCast(scalar.sind(@as(f64, @bitCast(input_bits))));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "sind exports a nearest-even binary64 function" {
    comptime std.debug.assert(@TypeOf(scalar.sind) == fn (f64) f64);
}

test "sind preserves signed zero and quiets NaNs" {
    try expectBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
    try expectBits(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
}

test "sind maps infinities to the canonical positive NaN" {
    try expectBits(0x7ff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
}

test "sind reduces ordinary and quadrant-boundary arguments" {
    try expectBits(0x3ff0_0000_0000_0000, 0x3fea_ed54_8f09_0cee);
    try expectBits(0xbff0_0000_0000_0000, 0xbfea_ed54_8f09_0cee);
    try expectBits(0x3ff9_21fb_5444_2d18, 0x3ff0_0000_0000_0000);
    try expectBits(0x4009_21fb_5444_2d18, 0x3ca1_a626_3314_5c07);
    try expectBits(0x4019_21fb_5444_2d18, 0xbcb1_a626_3314_5c07);
}

test "sind handles accurate fallback hard cases" {
    try expectBits(0x3ebe_0000_0000_01c2, 0x3ebd_ffff_ffff_f02e);
    try expectBits(0x7506_ac5b_262c_a1ff, 0x3ff0_0000_0000_0000);
    try expectBits(0xf506_ac5b_262c_a1ff, 0xbff0_0000_0000_0000);
}

test "sind uses identity rounding through the subnormal boundary" {
    try expectBits(0x3e57_1374_4912_3ef6, 0x3e57_1374_4912_3ef6);
    try expectBits(0x3e57_1374_4912_3ef7, 0x3e57_1374_4912_3ef6);
    try expectBits(0x0000_0000_0000_0001, 0x0000_0000_0000_0001);
    try expectBits(0x000f_ffff_ffff_ffff, 0x000f_ffff_ffff_ffff);
    try expectBits(0x8010_0000_0000_0000, 0x8010_0000_0000_0000);
}
