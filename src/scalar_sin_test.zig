const std = @import("std");
const scalar = @import("scalar.zig");

fn expectSinhBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.sinh(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "sinh exports the pinned binary16 sine signature" {
    comptime std.debug.assert(@TypeOf(scalar.sinh) == fn (f16) f16);
}

test "sinh preserves signed zero and ordinary finite values" {
    try expectSinhBits(0x0000, 0x0000);
    try expectSinhBits(0x8000, 0x8000);
    try expectSinhBits(0x0001, 0x0001);
    try expectSinhBits(0x8001, 0x8001);
    try expectSinhBits(0x3c00, 0x3abb);
    try expectSinhBits(0xbc00, 0xbabb);
}

test "sinh returns a canonical NaN for both infinities" {
    try expectSinhBits(0x7c00, 0x7e00);
    try expectSinhBits(0xfc00, 0x7e00);
}

test "sinh quiets NaNs while preserving sign and payload" {
    try expectSinhBits(0x7e55, 0x7e55);
    try expectSinhBits(0x7c55, 0x7e55);
    try expectSinhBits(0xfe55, 0xfe55);
    try expectSinhBits(0xfc55, 0xfe55);
}

test "sinh returns the two pinned exceptional values" {
    try expectSinhBits(0x658c, 0x07e7);
    try expectSinhBits(0xe58c, 0x87e7);
    try expectSinhBits(0x6429, 0x85ed);
    try expectSinhBits(0xe429, 0x05ed);
}

fn expectSinfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.sinf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "sinf exposes the nearest-even binary32 signature" {
    comptime std.debug.assert(@TypeOf(scalar.sinf) == fn (f32) f32);
}

test "sinf preserves signed zero" {
    try expectSinfBits(0x0000_0000, 0x0000_0000);
    try expectSinfBits(0x8000_0000, 0x8000_0000);
}

test "sinf quiets NaNs while preserving sign and payload" {
    try expectSinfBits(0x7fc1_2345, 0x7fc1_2345);
    try expectSinfBits(0xff81_2345, 0xffc1_2345);
}

test "sinf maps both infinities to the canonical quiet NaN" {
    try expectSinfBits(0x7f80_0000, 0x7fc0_0000);
    try expectSinfBits(0xff80_0000, 0x7fc0_0000);
}

test "sinf handles small and medium reduction boundaries" {
    try expectSinfBits(0x0000_0001, 0x0000_0001);
    try expectSinfBits(0x8000_0001, 0x8000_0001);
    try expectSinfBits(0x397f_ffff, 0x397f_ffff);
    try expectSinfBits(0x3980_0000, 0x3980_0000);
    try expectSinfBits(0x3f80_0000, 0x3f57_6aa4);
    try expectSinfBits(0xbf80_0000, 0xbf57_6aa4);
}

test "sinf covers CORE-MATH hard cases" {
    try expectSinfBits(0x4619_9998, 0xbeb1_fa5d);
    try expectSinfBits(0x3f3a_dc51, 0x3f2a_b445);
    try expectSinfBits(0x3fa7_832a, 0x3f77_41b6);
    try expectSinfBits(0x4116_cbe4, 0xb2cc_de2e);
}

test "sinf handles huge reduction and the largest finite inputs" {
    try expectSinfBits(0x4c80_0001, 0xbf75_6fe7);
    try expectSinfBits(0x7f7f_ffff, 0xbf05_99b3);
    try expectSinfBits(0xff7f_ffff, 0x3f05_99b3);
}
