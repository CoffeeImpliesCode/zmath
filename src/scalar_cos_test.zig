const std = @import("std");
const scalar = @import("scalar.zig");

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
