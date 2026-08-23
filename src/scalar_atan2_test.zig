const std = @import("std");
const scalar = @import("scalar.zig");

fn expectAtan2hBits(y_bits: u16, x_bits: u16, expected_bits: u16) !void {
    const y: f16 = @bitCast(y_bits);
    const x: f16 = @bitCast(x_bits);
    const actual_bits: u16 = @bitCast(scalar.atan2h(y, x));
    try std.testing.expectEqual(expected_bits, actual_bits);
}
fn expectAtan2fBits(y_bits: u32, x_bits: u32, expected_bits: u32) !void {
    const y: f32 = @bitCast(y_bits);
    const x: f32 = @bitCast(x_bits);
    const actual_bits: u32 = @bitCast(scalar.atan2f(y, x));
    try std.testing.expectEqual(expected_bits, actual_bits);
}
fn expectAtan2dBits(y_bits: u64, x_bits: u64, expected_bits: u64) !void {
    const y: f64 = @bitCast(y_bits);
    const x: f64 = @bitCast(x_bits);
    const actual_bits: u64 = @bitCast(scalar.atan2d(y, x));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "atan2 family handles special cases" {
    try expectAtan2hBits(0x7c01, 0x0000, 0x7e01);
    try expectAtan2fBits(0x7f800001, 0x00000000, 0x7fc00001);
    try expectAtan2dBits(0x7ff0000000000001, 0, 0x7ff8000000000001);
    try expectAtan2hBits(0x7c00, 0x7c00, 0x40b6);
    try expectAtan2hBits(0xfc00, 0x7c00, 0xc0b6);
    try expectAtan2hBits(0x7c00, 0xfc00, 0x40b6);
    try expectAtan2hBits(0xfc00, 0xfc00, 0xc0b6);
    try expectAtan2hBits(0x7c00, 0x3c00, 0x3e48);
    try expectAtan2hBits(0xfc00, 0x3c00, 0xbe48);
    try expectAtan2hBits(0x3c00, 0x7c00, 0x0000);
    try expectAtan2hBits(0x3c00, 0xfc00, 0x4248);
    try expectAtan2hBits(0xbc00, 0xfc00, 0xc248);
    try expectAtan2hBits(0x0000, 0x0000, 0x0000);
    try expectAtan2hBits(0x8000, 0x0000, 0x8000);
    try expectAtan2hBits(0x0000, 0x8000, 0x4248);
    try expectAtan2hBits(0x8000, 0x8000, 0xc248);
    try expectAtan2hBits(0x0000, 0x3c00, 0x0000);
    try expectAtan2hBits(0x0000, 0xbc00, 0x4248);
    try expectAtan2hBits(0x8000, 0x3c00, 0x8000);
    try expectAtan2hBits(0x8000, 0xbc00, 0xc248);
    try expectAtan2hBits(0x3c00, 0x0000, 0x3e48);
    try expectAtan2hBits(0xbc00, 0x0000, 0xbe48);
}
