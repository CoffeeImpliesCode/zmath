const std = @import("std");
const scalar = @import("scalar.zig");

fn expectSqrtHBits(input: u16, expected: u16) !void {
    const value: f16 = @bitCast(input);
    try std.testing.expectEqual(expected, @as(u16, @bitCast(scalar.sqrth(value))));
}

fn expectSqrtFBits(input: u32, expected: u32) !void {
    const value: f32 = @bitCast(input);
    try std.testing.expectEqual(expected, @as(u32, @bitCast(scalar.sqrtf(value))));
}

fn expectSqrtDBits(input: u64, expected: u64) !void {
    const value: f64 = @bitCast(input);
    try std.testing.expectEqual(expected, @as(u64, @bitCast(scalar.sqrtd(value))));
}

test "Square root returns exact perfect squares" {
    // f16
    try expectSqrtHBits(0x0000, 0x0000); // sqrt(0) = 0
    try expectSqrtHBits(0x3c00, 0x3c00); // sqrt(1) = 1
    try expectSqrtHBits(0x4400, 0x4000); // sqrt(4) = 2
    try expectSqrtHBits(0x4880, 0x4200); // sqrt(9) = 3
    try expectSqrtHBits(0x4c00, 0x4400); // sqrt(16) = 4

    // f32
    try expectSqrtFBits(0x0000_0000, 0x0000_0000); // sqrt(0) = 0
    try expectSqrtFBits(0x3f80_0000, 0x3f80_0000); // sqrt(1) = 1
    try expectSqrtFBits(0x4080_0000, 0x4000_0000); // sqrt(4) = 2
    try expectSqrtFBits(0x4110_0000, 0x4040_0000); // sqrt(9) = 3
    try expectSqrtFBits(0x41c8_0000, 0x40a0_0000); // sqrt(25) = 5

    // f64
    try expectSqrtDBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000); // sqrt(0) = 0
    try expectSqrtDBits(0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000); // sqrt(1) = 1
    try expectSqrtDBits(0x4010_0000_0000_0000, 0x4000_0000_0000_0000); // sqrt(4) = 2
    try expectSqrtDBits(0x4022_0000_0000_0000, 0x4008_0000_0000_0000); // sqrt(9) = 3
    try expectSqrtDBits(0x4039_0000_0000_0000, 0x4014_0000_0000_0000); // sqrt(25) = 5
}

test "Square root preserves sign of zero" {
    try expectSqrtHBits(0x8000, 0x8000); // -0 -> -0
    try expectSqrtFBits(0x8000_0000, 0x8000_0000); // -0 -> -0
    try expectSqrtDBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000); // -0 -> -0
}

test "Square root handles +Inf" {
    try expectSqrtHBits(0x7c00, 0x7c00); // +Inf -> +Inf
    try expectSqrtFBits(0x7f80_0000, 0x7f80_0000); // +Inf -> +Inf
    try expectSqrtDBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000); // +Inf -> +Inf
}

test "Square root handles negative finite and -Inf as canonical +qNaN" {
    try expectSqrtHBits(0xbc00, 0x7e00); // -1.0 -> canonical qNaN
    try expectSqrtHBits(0xfc00, 0x7e00); // -Inf -> canonical qNaN
    try expectSqrtHBits(0x8200, 0x7e00); // -subnormal -> canonical qNaN
    try expectSqrtFBits(0x8000_0001, 0x7fc0_0000); // -subnormal
    try expectSqrtFBits(0xbf80_0000, 0x7fc0_0000); // -1.0 -> canonical qNaN
    try expectSqrtFBits(0xff80_0000, 0x7fc0_0000); // -Inf -> canonical qNaN
    try expectSqrtDBits(0xbff0_0000_0000_0000, 0x7ff8_0000_0000_0000); // -1.0 -> canonical qNaN
    try expectSqrtDBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000); // -Inf -> canonical qNaN
    try expectSqrtDBits(0x8000_0000_0000_0001, 0x7ff8_0000_0000_0000); // -subnormal
}

test "Square root quiets sNaN and preserves qNaN payload" {
    // f16 signaling NaN (quiet bit clear): 0x7d00 -> 0x7f00 (quieted)
    try expectSqrtHBits(0x7d00, 0x7f00);
    // f16 quiet NaN: 0x7e55 -> 0x7e55 (preserved)
    try expectSqrtHBits(0x7e55, 0x7e55);
    try expectSqrtHBits(0xfd01, 0xff01); // negative signaling NaN
    try expectSqrtHBits(0xfe55, 0xfe55); // negative quiet NaN
    // f32 signaling NaN: 0x7f800001 -> 0x7fc00001 (quieted)
    try expectSqrtFBits(0x7f80_0001, 0x7fc0_0001);
    // f32 quiet NaN: 0x7fc0_0123 -> 0x7fc0_0123 (preserved)
    try expectSqrtFBits(0x7fc0_0123, 0x7fc0_0123);
    try expectSqrtFBits(0xff80_0011, 0xffc0_0011); // negative signaling NaN
    try expectSqrtFBits(0xffc0_0123, 0xffc0_0123); // negative quiet NaN
    // f64 signaling NaN: 0x7ff0_0000_0000_0001 -> 0x7ff8_0000_0000_0001 (quieted)
    try expectSqrtDBits(0x7ff0_0000_0000_0001, 0x7ff8_0000_0000_0001);
    // f64 quiet NaN: 0x7ff8_0000_0000_0123 -> 0x7ff8_0000_0000_0123 (preserved)
    try expectSqrtDBits(0x7ff8_0000_0000_0123, 0x7ff8_0000_0000_0123);
    try expectSqrtDBits(0xfff0_0000_0000_0011, 0xfff8_0000_0000_0011); // negative signaling NaN
    try expectSqrtDBits(0xfff8_0000_0000_0123, 0xfff8_0000_0000_0123); // negative quiet NaN
}

test "Square root handles subnormal inputs" {
    // f16: sqrt(2^-24) = 2^-12 (exact), input is smallest subnormal
    try expectSqrtHBits(0x0001, 0x0c00);
    // f16: sqrt(2^-16) = 2^-8 (exact), input 0x0100 = subnormal 2^-16
    try expectSqrtHBits(0x0100, 0x1c00);

    // f32: sqrt(2^-149) = 2^-74.5 ≈ 3.74e-23, nearest f32 = 0x1a3504f3
    try expectSqrtFBits(0x0000_0001, 0x1a35_04f3);

    // f64: sqrt(2^-1074) = 2^-537 (exact power), nearest f64 = 0x1e60000000000000
    try expectSqrtDBits(0x0000_0000_0000_0001, 0x1e60_0000_0000_0000);
    // Minimum positive normal inputs are exact squares in every format.
    try expectSqrtHBits(0x0400, 0x2000);
    try expectSqrtFBits(0x0080_0000, 0x2000_0000);
    try expectSqrtDBits(0x0010_0000_0000_0000, 0x2000_0000_0000_0000);
}

test "Square root computes correctly rounded result for sqrt(2)" {
    // f16: sqrt(2) rounds to 0x3da8 (1.4140625)
    try expectSqrtHBits(0x4000, 0x3da8);

    // f32 and f64: compare against @sqrt which is our hardware source
    const sqrt2_f32: f32 = @sqrt(2.0);
    try std.testing.expectEqual(
        @as(u32, @bitCast(sqrt2_f32)),
        @as(u32, @bitCast(scalar.sqrtf(@as(f32, 2.0)))),
    );
    const sqrt2_f64: f64 = @sqrt(2.0);
    try std.testing.expectEqual(
        @as(u64, @bitCast(sqrt2_f64)),
        @as(u64, @bitCast(scalar.sqrtd(@as(f64, 2.0)))),
    );
}
