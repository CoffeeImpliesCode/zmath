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
