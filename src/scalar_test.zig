const std = @import("std");
const scalar = @import("scalar.zig");
const ieee = @import("ieee.zig");

test "absF returns the magnitude of a negative finite value" {
    try std.testing.expectEqual(@as(f32, 3.5), scalar.absF(-3.5));
}

test "absF maps negative zero to positive zero" {
    const negative_zero: f32 = @bitCast(@as(u32, 0x8000_0000));
    const result_bits: u32 = @bitCast(scalar.absF(negative_zero));

    try std.testing.expectEqual(@as(u32, 0x0000_0000), result_bits);
}

test "absD returns the magnitude of a negative finite value" {
    try std.testing.expectEqual(@as(f64, 3.5), scalar.absD(-3.5));
}

test "absD maps negative zero to positive zero" {
    const negative_zero: f64 = @bitCast(@as(u64, 0x8000_0000_0000_0000));
    const result_bits: u64 = @bitCast(scalar.absD(negative_zero));

    try std.testing.expectEqual(@as(u64, 0x0000_0000_0000_0000), result_bits);
}

test "ulpDistanceF returns zero for equal values" {
    try std.testing.expectEqual(@as(u32, 0), scalar.ulpDistanceF(1.5, 1.5));
}

test "ulpDistanceF returns one for adjacent finite values" {
    const adjacent = ieee.nextAfterF(0.0, 1.0);

    try std.testing.expectEqual(@as(u32, 1), scalar.ulpDistanceF(0.0, adjacent));
}

test "ulpDistanceF treats opposite signed zeros as equal" {
    const negative_zero: f32 = @bitCast(@as(u32, 0x8000_0000));
    const positive_zero: f32 = @bitCast(@as(u32, 0x0000_0000));

    try std.testing.expectEqual(@as(u32, 0), scalar.ulpDistanceF(negative_zero, positive_zero));
}

test "ulpDistanceD returns zero for equal values" {
    try std.testing.expectEqual(@as(u64, 0), scalar.ulpDistanceD(1.5, 1.5));
}

test "ulpDistanceD returns one for adjacent finite values" {
    const adjacent = ieee.nextAfterD(0.0, 1.0);

    try std.testing.expectEqual(@as(u64, 1), scalar.ulpDistanceD(0.0, adjacent));
}

test "ulpDistanceD treats opposite signed zeros as equal" {
    const negative_zero: f64 = @bitCast(@as(u64, 0x8000_0000_0000_0000));
    const positive_zero: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0000));

    try std.testing.expectEqual(@as(u64, 0), scalar.ulpDistanceD(negative_zero, positive_zero));
}
