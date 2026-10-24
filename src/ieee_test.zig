const std = @import("std");
const ieee = @import("ieee.zig");

test "RoundMode defines the four public rounding modes" {
    const modes = [_]ieee.RoundMode{
        .nearest_even,
        .toward_zero,
        .toward_positive,
        .toward_negative,
    };

    try std.testing.expectEqual(modes.len, @typeInfo(ieee.RoundMode).@"enum".field_names.len);
}

test "classifyf identifies every public IEEE class" {
    const cases = [_]struct {
        bits: u32,
        expected: ieee.Class,
    }{
        .{ .bits = 0x0000_0000, .expected = .zero },
        .{ .bits = 0x8000_0000, .expected = .zero },
        .{ .bits = 0x0000_0001, .expected = .subnormal },
        .{ .bits = 0x3f80_0000, .expected = .normal },
        .{ .bits = 0x7f80_0000, .expected = .infinity },
        .{ .bits = 0x7fc0_0000, .expected = .quiet_nan },
        .{ .bits = 0x7f80_0001, .expected = .signaling_nan },
    };

    for (cases) |case| {
        const value: f32 = @bitCast(case.bits);
        try std.testing.expectEqual(case.expected, ieee.classifyf(value));
    }
}

test "classifyd identifies every public IEEE class" {
    const cases = [_]struct {
        bits: u64,
        expected: ieee.Class,
    }{
        .{ .bits = 0x0000_0000_0000_0000, .expected = .zero },
        .{ .bits = 0x8000_0000_0000_0000, .expected = .zero },
        .{ .bits = 0x0000_0000_0000_0001, .expected = .subnormal },
        .{ .bits = 0x3ff0_0000_0000_0000, .expected = .normal },
        .{ .bits = 0x7ff0_0000_0000_0000, .expected = .infinity },
        .{ .bits = 0x7ff8_0000_0000_0000, .expected = .quiet_nan },
        .{ .bits = 0x7ff0_0000_0000_0001, .expected = .signaling_nan },
    };

    for (cases) |case| {
        const value: f64 = @bitCast(case.bits);
        try std.testing.expectEqual(case.expected, ieee.classifyd(value));
    }
}

test "copySignf copies only the sign bit" {
    const positive: f32 = @bitCast(@as(u32, 0x3fc0_0000));
    const negative: f32 = @bitCast(@as(u32, 0xbfc0_0000));
    const positive_sign: f32 = @bitCast(@as(u32, 0x0000_0000));
    const negative_sign: f32 = @bitCast(@as(u32, 0x8000_0000));

    try std.testing.expectEqual(
        @as(u32, 0xbfc0_0000),
        @as(u32, @bitCast(ieee.copySignf(positive, negative_sign))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3fc0_0000),
        @as(u32, @bitCast(ieee.copySignf(negative, positive_sign))),
    );
}

test "copySignd copies only the sign bit" {
    const positive: f64 = @bitCast(@as(u64, 0x3ff8_0000_0000_0000));
    const negative: f64 = @bitCast(@as(u64, 0xbff8_0000_0000_0000));
    const positive_sign: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0000));
    const negative_sign: f64 = @bitCast(@as(u64, 0x8000_0000_0000_0000));

    try std.testing.expectEqual(
        @as(u64, 0xbff8_0000_0000_0000),
        @as(u64, @bitCast(ieee.copySignd(positive, negative_sign))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3ff8_0000_0000_0000),
        @as(u64, @bitCast(ieee.copySignd(negative, positive_sign))),
    );
}

test "nextAfterf advances signed zero toward positive and negative one" {
    const positive_zero: f32 = @bitCast(@as(u32, 0x0000_0000));
    const negative_zero: f32 = @bitCast(@as(u32, 0x8000_0000));
    const positive_one: f32 = @bitCast(@as(u32, 0x3f80_0000));
    const negative_one: f32 = @bitCast(@as(u32, 0xbf80_0000));

    try std.testing.expectEqual(
        @as(u32, 0x0000_0001),
        @as(u32, @bitCast(ieee.nextAfterf(positive_zero, positive_one))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x0000_0001),
        @as(u32, @bitCast(ieee.nextAfterf(negative_zero, positive_one))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x8000_0001),
        @as(u32, @bitCast(ieee.nextAfterf(positive_zero, negative_one))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x8000_0001),
        @as(u32, @bitCast(ieee.nextAfterf(negative_zero, negative_one))),
    );
}

test "nextAfterd advances signed zero toward positive and negative one" {
    const positive_zero: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0000));
    const negative_zero: f64 = @bitCast(@as(u64, 0x8000_0000_0000_0000));
    const positive_one: f64 = @bitCast(@as(u64, 0x3ff0_0000_0000_0000));
    const negative_one: f64 = @bitCast(@as(u64, 0xbff0_0000_0000_0000));

    try std.testing.expectEqual(
        @as(u64, 0x0000_0000_0000_0001),
        @as(u64, @bitCast(ieee.nextAfterd(positive_zero, positive_one))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x0000_0000_0000_0001),
        @as(u64, @bitCast(ieee.nextAfterd(negative_zero, positive_one))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x8000_0000_0000_0001),
        @as(u64, @bitCast(ieee.nextAfterd(positive_zero, negative_one))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x8000_0000_0000_0001),
        @as(u64, @bitCast(ieee.nextAfterd(negative_zero, negative_one))),
    );
}

test "totalOrderf orders negative zero before positive zero" {
    const positive_zero: f32 = @bitCast(@as(u32, 0x0000_0000));
    const negative_zero: f32 = @bitCast(@as(u32, 0x8000_0000));

    try std.testing.expect(ieee.totalOrderf(negative_zero, positive_zero));
    try std.testing.expect(!ieee.totalOrderf(positive_zero, negative_zero));
}

test "totalOrderd orders negative zero before positive zero" {
    const positive_zero: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0000));
    const negative_zero: f64 = @bitCast(@as(u64, 0x8000_0000_0000_0000));

    try std.testing.expect(ieee.totalOrderd(negative_zero, positive_zero));
    try std.testing.expect(!ieee.totalOrderd(positive_zero, negative_zero));
}
