const std = @import("std");
const scalar = @import("scalar.zig");
const ieee = @import("ieee.zig");

test "absf returns the magnitude of a negative finite value" {
    try std.testing.expectEqual(@as(f32, 3.5), scalar.absf(-3.5));
}

test "absf maps negative zero to positive zero" {
    const negative_zero: f32 = @bitCast(@as(u32, 0x8000_0000));
    const result_bits: u32 = @bitCast(scalar.absf(negative_zero));

    try std.testing.expectEqual(@as(u32, 0x0000_0000), result_bits);
}

test "absd returns the magnitude of a negative finite value" {
    try std.testing.expectEqual(@as(f64, 3.5), scalar.absd(-3.5));
}

test "absd maps negative zero to positive zero" {
    const negative_zero: f64 = @bitCast(@as(u64, 0x8000_0000_0000_0000));
    const result_bits: u64 = @bitCast(scalar.absd(negative_zero));

    try std.testing.expectEqual(@as(u64, 0x0000_0000_0000_0000), result_bits);
}

test "ulpDistancef returns zero for equal values" {
    try std.testing.expectEqual(@as(u32, 0), scalar.ulpDistancef(1.5, 1.5));
}

test "ulpDistancef returns one for adjacent finite values" {
    const adjacent = ieee.nextAfterf(0.0, 1.0);

    try std.testing.expectEqual(@as(u32, 1), scalar.ulpDistancef(0.0, adjacent));
}

test "ulpDistancef treats opposite signed zeros as equal" {
    const negative_zero: f32 = @bitCast(@as(u32, 0x8000_0000));
    const positive_zero: f32 = @bitCast(@as(u32, 0x0000_0000));

    try std.testing.expectEqual(@as(u32, 0), scalar.ulpDistancef(negative_zero, positive_zero));
}

test "ulpDistanced returns zero for equal values" {
    try std.testing.expectEqual(@as(u64, 0), scalar.ulpDistanced(1.5, 1.5));
}

test "ulpDistanced returns one for adjacent finite values" {
    const adjacent = ieee.nextAfterd(0.0, 1.0);

    try std.testing.expectEqual(@as(u64, 1), scalar.ulpDistanced(0.0, adjacent));
}

test "ulpDistanced treats opposite signed zeros as equal" {
    const negative_zero: f64 = @bitCast(@as(u64, 0x8000_0000_0000_0000));
    const positive_zero: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0000));

    try std.testing.expectEqual(@as(u64, 0), scalar.ulpDistanced(negative_zero, positive_zero));
}
