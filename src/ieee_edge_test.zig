const std = @import("std");
const ieee = @import("ieee.zig");

test "nextAfterf quiets a signaling NaN without changing its sign or payload" {
    const signaling_nan: f32 = @bitCast(@as(u32, 0xff81_2345));
    const toward_zero: f32 = @bitCast(@as(u32, 0x0000_0000));

    try std.testing.expectEqual(
        @as(u32, 0xffc1_2345),
        @as(u32, @bitCast(ieee.nextAfterf(signaling_nan, toward_zero))),
    );
}

test "nextAfterd quiets a signaling NaN without changing its sign or payload" {
    const signaling_nan: f64 = @bitCast(@as(u64, 0xfff0_1234_5678_9abc));
    const toward_zero: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0000));

    try std.testing.expectEqual(
        @as(u64, 0xfff8_1234_5678_9abc),
        @as(u64, @bitCast(ieee.nextAfterd(signaling_nan, toward_zero))),
    );
}

test "nextAfterf steps infinities toward zero" {
    const positive_infinity: f32 = @bitCast(@as(u32, 0x7f80_0000));
    const negative_infinity: f32 = @bitCast(@as(u32, 0xff80_0000));
    const zero: f32 = @bitCast(@as(u32, 0x0000_0000));

    try std.testing.expectEqual(
        @as(u32, 0x7f7f_ffff),
        @as(u32, @bitCast(ieee.nextAfterf(positive_infinity, zero))),
    );
    try std.testing.expectEqual(
        @as(u32, 0xff7f_ffff),
        @as(u32, @bitCast(ieee.nextAfterf(negative_infinity, zero))),
    );
}

test "nextAfterd steps infinities toward zero" {
    const positive_infinity: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0000));
    const negative_infinity: f64 = @bitCast(@as(u64, 0xfff0_0000_0000_0000));
    const zero: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0000));

    try std.testing.expectEqual(
        @as(u64, 0x7fef_ffff_ffff_ffff),
        @as(u64, @bitCast(ieee.nextAfterd(positive_infinity, zero))),
    );
    try std.testing.expectEqual(
        @as(u64, 0xffef_ffff_ffff_ffff),
        @as(u64, @bitCast(ieee.nextAfterd(negative_infinity, zero))),
    );
}

test "nextAfterf returns the destination bit pattern for equal signed zeros" {
    const positive_zero: f32 = @bitCast(@as(u32, 0x0000_0000));
    const negative_zero: f32 = @bitCast(@as(u32, 0x8000_0000));

    try std.testing.expectEqual(
        @as(u32, 0x8000_0000),
        @as(u32, @bitCast(ieee.nextAfterf(positive_zero, negative_zero))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x0000_0000),
        @as(u32, @bitCast(ieee.nextAfterf(negative_zero, positive_zero))),
    );
}

test "nextAfterd returns the destination bit pattern for equal signed zeros" {
    const positive_zero: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0000));
    const negative_zero: f64 = @bitCast(@as(u64, 0x8000_0000_0000_0000));

    try std.testing.expectEqual(
        @as(u64, 0x8000_0000_0000_0000),
        @as(u64, @bitCast(ieee.nextAfterd(positive_zero, negative_zero))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x0000_0000_0000_0000),
        @as(u64, @bitCast(ieee.nextAfterd(negative_zero, positive_zero))),
    );
}
