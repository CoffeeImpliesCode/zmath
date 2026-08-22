const std = @import("std");
const scalar = @import("scalar.zig");

test "ULP distance across signed minimum subnormals is two" {
    const negative_minimum_f: f32 = @bitCast(@as(u32, 0x8000_0001));
    const positive_minimum_f: f32 = @bitCast(@as(u32, 0x0000_0001));
    const negative_minimum_d: f64 = @bitCast(@as(u64, 0x8000_0000_0000_0001));
    const positive_minimum_d: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0001));

    try std.testing.expectEqual(
        @as(u32, 2),
        scalar.ulpDistancef(negative_minimum_f, positive_minimum_f),
    );
    try std.testing.expectEqual(
        @as(u64, 2),
        scalar.ulpDistanced(negative_minimum_d, positive_minimum_d),
    );
}

test "ULP distance across infinities spans the finite encoding range" {
    const negative_infinity_f: f32 = @bitCast(@as(u32, 0xff80_0000));
    const positive_infinity_f: f32 = @bitCast(@as(u32, 0x7f80_0000));
    const negative_infinity_d: f64 = @bitCast(@as(u64, 0xfff0_0000_0000_0000));
    const positive_infinity_d: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0000));

    try std.testing.expectEqual(
        @as(u32, 0xff00_0000),
        scalar.ulpDistancef(negative_infinity_f, positive_infinity_f),
    );
    try std.testing.expectEqual(
        @as(u64, 0xffe0_0000_0000_0000),
        scalar.ulpDistanced(negative_infinity_d, positive_infinity_d),
    );
}

test "ULP distance from NaN is the maximum unsigned sentinel" {
    const nan_f: f32 = @bitCast(@as(u32, 0x7fc0_0001));
    const nan_d: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_0001));

    try std.testing.expectEqual(
        @as(u32, 0xffff_ffff),
        scalar.ulpDistancef(nan_f, 0.0),
    );
    try std.testing.expectEqual(
        @as(u64, 0xffff_ffff_ffff_ffff),
        scalar.ulpDistanced(nan_d, 0.0),
    );
}
