const std = @import("std");
const scalar = @import("scalar.zig");

test "fmaf NaN handling" {
    const quiet_nan_payload: f32 = @bitCast(@as(u32, 0x7fc00000 | 0x00012345));
    const quiet_nan_neg: f32 = @bitCast(@as(u32, 0xffc00000 | 0x00012345));
    {
        const result = scalar.fmaf(quiet_nan_payload, 1.0, 2.0);
        try std.testing.expectEqual(@as(u32, 0x7fc00000 | 0x00012345), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(1.0, quiet_nan_payload, 2.0);
        try std.testing.expectEqual(@as(u32, 0x7fc00000 | 0x00012345), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(1.0, 2.0, quiet_nan_payload);
        try std.testing.expectEqual(@as(u32, 0x7fc00000 | 0x00012345), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(quiet_nan_neg, 1.0, 2.0);
        try std.testing.expectEqual(@as(u32, 0xffc00000 | 0x00012345), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(1.0, quiet_nan_neg, 2.0);
        try std.testing.expectEqual(@as(u32, 0xffc00000 | 0x00012345), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(1.0, 2.0, quiet_nan_neg);
        try std.testing.expectEqual(@as(u32, 0xffc00000 | 0x00012345), @as(u32, @bitCast(result)));
    }
}

test "fmaf inf * zero -> nan" {
    const inf: f32 = @bitCast(@as(u32, 0x7f800000));
    const neg_inf: f32 = @bitCast(@as(u32, 0xff800000));
    const zero: f32 = 0.0;
    const neg_zero: f32 = @bitCast(@as(u32, 0x80000000));
    {
        const result = scalar.fmaf(inf, zero, 1.0);
        try std.testing.expectEqual(@as(u32, 0x7fc00000), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(neg_inf, zero, 1.0);
        try std.testing.expectEqual(@as(u32, 0x7fc00000), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(inf, neg_zero, 1.0);
        try std.testing.expectEqual(@as(u32, 0x7fc00000), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(neg_inf, neg_zero, 1.0);
        try std.testing.expectEqual(@as(u32, 0x7fc00000), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(zero, inf, 1.0);
        try std.testing.expectEqual(@as(u32, 0x7fc00000), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(neg_zero, inf, 1.0);
        try std.testing.expectEqual(@as(u32, 0x7fc00000), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(zero, neg_inf, 1.0);
        try std.testing.expectEqual(@as(u32, 0x7fc00000), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(neg_zero, neg_inf, 1.0);
        try std.testing.expectEqual(@as(u32, 0x7fc00000), @as(u32, @bitCast(result)));
    }
}

test "fmaf product inf + inf opposite sign -> nan" {
    const inf: f32 = @bitCast(@as(u32, 0x7f800000));
    const neg_inf: f32 = @bitCast(@as(u32, 0xff800000));
    {
        const result = scalar.fmaf(inf, 1.0, neg_inf);
        try std.testing.expectEqual(@as(u32, 0x7fc00000), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(neg_inf, 1.0, inf);
        try std.testing.expectEqual(@as(u32, 0x7fc00000), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(inf, -1.0, inf);
        try std.testing.expectEqual(@as(u32, 0x7fc00000), @as(u32, @bitCast(result)));
    }
    {
        const result = scalar.fmaf(neg_inf, -1.0, neg_inf);
        try std.testing.expectEqual(@as(u32, 0x7fc00000), @as(u32, @bitCast(result)));
    }
}

test "fmaf normal values" {
    const a: f32 = 1.0;
    const b: f32 = 2.0;
    const c: f32 = 3.0;
    try std.testing.expectEqual(@mulAdd(f32, a, b, c), scalar.fmaf(a, b, c));
    try std.testing.expectEqual(@mulAdd(f32, -1.0, 2.0, 3.0), scalar.fmaf(-1.0, 2.0, 3.0));
    try std.testing.expectEqual(@mulAdd(f32, 1.0, -2.0, 3.0), scalar.fmaf(1.0, -2.0, 3.0));
    try std.testing.expectEqual(@mulAdd(f32, 1.0, 2.0, -3.0), scalar.fmaf(1.0, 2.0, -3.0));
    try std.testing.expectEqual(@mulAdd(f32, 0.0, 2.0, 3.0), scalar.fmaf(0.0, 2.0, 3.0));
    try std.testing.expectEqual(@mulAdd(f32, 1.0, 0.0, 3.0), scalar.fmaf(1.0, 0.0, 3.0));
    try std.testing.expectEqual(@mulAdd(f32, 1.0, 2.0, 0.0), scalar.fmaf(1.0, 2.0, 0.0));
}
