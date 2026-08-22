const std = @import("std");
const scalar = @import("scalar.zig");

fn expectCoshBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.cosh(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "cosh exposes the nearest-even binary16 signature" {
    comptime std.debug.assert(@TypeOf(scalar.cosh) == fn (f16) f16);
}

test "cosh maps both signed zeros and tiny values to one" {
    try expectCoshBits(0x0000, 0x3c00);
    try expectCoshBits(0x8000, 0x3c00);
    try expectCoshBits(0x0001, 0x3c00);
    try expectCoshBits(0x8001, 0x3c00);
}

test "cosh is even for ordinary and reduced finite values" {
    try expectCoshBits(0x3c00, 0x3853);
    try expectCoshBits(0xbc00, 0x3853);
    try expectCoshBits(0x4248, 0xbc00);
    try expectCoshBits(0xc248, 0xbc00);
    try expectCoshBits(0x7bff, 0xb30f);
    try expectCoshBits(0xfbff, 0xb30f);
}

test "cosh covers the two finite underflow results" {
    try expectCoshBits(0x598c, 0x80fd);
    try expectCoshBits(0xd98c, 0x80fd);
    try expectCoshBits(0x6029, 0x02f7);
    try expectCoshBits(0xe029, 0x02f7);
}

test "cosh returns a canonical NaN for both infinities" {
    try expectCoshBits(0x7c00, 0x7e00);
    try expectCoshBits(0xfc00, 0x7e00);
}

test "cosh quiets NaNs while preserving sign and payload" {
    try expectCoshBits(0x7e55, 0x7e55);
    try expectCoshBits(0x7c55, 0x7e55);
    try expectCoshBits(0xfe55, 0xfe55);
    try expectCoshBits(0xfc55, 0xfe55);
}

fn expectCosfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.cosf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "cosf exposes the nearest-even binary32 signature" {
    comptime std.debug.assert(@TypeOf(scalar.cosf) == fn (f32) f32);
}

test "cosf maps both signed zero values to one" {
    try expectCosfBits(0x0000_0000, 0x3f80_0000);
    try expectCosfBits(0x8000_0000, 0x3f80_0000);
}

test "cosf quiets NaNs while preserving sign and payload" {
    try expectCosfBits(0x7fc1_2345, 0x7fc1_2345);
    try expectCosfBits(0x7f81_2345, 0x7fc1_2345);
    try expectCosfBits(0xffc1_2345, 0xffc1_2345);
    try expectCosfBits(0xff81_2345, 0xffc1_2345);
}

test "cosf maps both infinities to the canonical quiet NaN" {
    try expectCosfBits(0x7f80_0000, 0x7fc0_0000);
    try expectCosfBits(0xff80_0000, 0x7fc0_0000);
}

test "cosf handles tiny and ordinary finite values" {
    try expectCosfBits(0x0000_0001, 0x3f80_0000);
    try expectCosfBits(0x8000_0001, 0x3f80_0000);
    try expectCosfBits(0x397f_ffff, 0x3f80_0000);
    try expectCosfBits(0x3980_0001, 0x3f7f_ffff);
    try expectCosfBits(0x3f80_0000, 0x3f0a_5140);
    try expectCosfBits(0xbf80_0000, 0x3f0a_5140);
}

test "cosf covers every CORE-MATH hard case" {
    try expectCosfBits(0x4096_cbe4, 0x324c_de2e);
    try expectCosfBits(0xc096_cbe4, 0x324c_de2e);
    try expectCosfBits(0x5922_aa80, 0x3f08_aebf);
    try expectCosfBits(0x5aa4_542c, 0x3efa_40a4);
    try expectCosfBits(0x5f18_b878, 0x3f7f_14bb);
    try expectCosfBits(0x6115_cb11, 0x3f78_142f);
}

test "cosf handles huge reduction and the largest finite values" {
    try expectCosfBits(0x4c80_0001, 0xbe91_8e4a);
    try expectCosfBits(0x7f7f_ffff, 0x3f5a_5f96);
    try expectCosfBits(0xff7f_ffff, 0x3f5a_5f96);
}

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
