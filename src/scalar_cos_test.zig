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
