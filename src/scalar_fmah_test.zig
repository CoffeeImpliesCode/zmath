const std = @import("std");
const scalar = @import("scalar.zig");

test "fmah NaN handling" {
    const quiet_nan_payload: f16 = @bitCast(@as(u16, 0x7e00 | 0x0123));
    const quiet_nan_neg: f16 = @bitCast(@as(u16, 0xfe00 | 0x0123));
    {
        const result = scalar.fmah(quiet_nan_payload, 1.0, 2.0);
        try std.testing.expectEqual(@as(u16, 0x7e00 | 0x0123), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(1.0, quiet_nan_payload, 2.0);
        try std.testing.expectEqual(@as(u16, 0x7e00 | 0x0123), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(1.0, 2.0, quiet_nan_payload);
        try std.testing.expectEqual(@as(u16, 0x7e00 | 0x0123), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(quiet_nan_neg, 1.0, 2.0);
        try std.testing.expectEqual(@as(u16, 0xfe00 | 0x0123), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(1.0, quiet_nan_neg, 2.0);
        try std.testing.expectEqual(@as(u16, 0xfe00 | 0x0123), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(1.0, 2.0, quiet_nan_neg);
        try std.testing.expectEqual(@as(u16, 0xfe00 | 0x0123), @as(u16, @bitCast(result)));
    }
}

test "fmah inf * zero -> nan" {
    const inf: f16 = @bitCast(@as(u16, 0x7c00));
    const neg_inf: f16 = @bitCast(@as(u16, 0xfc00));
    const zero: f16 = 0.0;
    const neg_zero: f16 = @bitCast(@as(u16, 0x8000));
    {
        const result = scalar.fmah(inf, zero, 1.0);
        try std.testing.expectEqual(@as(u16, 0x7e00), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(neg_inf, zero, 1.0);
        try std.testing.expectEqual(@as(u16, 0x7e00), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(inf, neg_zero, 1.0);
        try std.testing.expectEqual(@as(u16, 0x7e00), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(neg_inf, neg_zero, 1.0);
        try std.testing.expectEqual(@as(u16, 0x7e00), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(zero, inf, 1.0);
        try std.testing.expectEqual(@as(u16, 0x7e00), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(neg_zero, inf, 1.0);
        try std.testing.expectEqual(@as(u16, 0x7e00), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(zero, neg_inf, 1.0);
        try std.testing.expectEqual(@as(u16, 0x7e00), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(neg_zero, neg_inf, 1.0);
        try std.testing.expectEqual(@as(u16, 0x7e00), @as(u16, @bitCast(result)));
    }
}

test "fmah product inf + inf opposite sign -> nan" {
    const inf: f16 = @bitCast(@as(u16, 0x7c00));
    const neg_inf: f16 = @bitCast(@as(u16, 0xfc00));
    {
        const result = scalar.fmah(inf, 1.0, neg_inf);
        try std.testing.expectEqual(@as(u16, 0x7e00), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(neg_inf, 1.0, inf);
        try std.testing.expectEqual(@as(u16, 0x7e00), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(inf, -1.0, inf);
        try std.testing.expectEqual(@as(u16, 0x7e00), @as(u16, @bitCast(result)));
    }
    {
        const result = scalar.fmah(neg_inf, -1.0, neg_inf);
        try std.testing.expectEqual(@as(u16, 0x7e00), @as(u16, @bitCast(result)));
    }
}

test "fmah normal values" {
    const a: f16 = 1.0;
    const b: f16 = 2.0;
    const c: f16 = 3.0;
    try std.testing.expectEqual(@mulAdd(f16, a, b, c), scalar.fmah(a, b, c));
    try std.testing.expectEqual(@mulAdd(f16, -1.0, 2.0, 3.0), scalar.fmah(-1.0, 2.0, 3.0));
    try std.testing.expectEqual(@mulAdd(f16, 1.0, -2.0, 3.0), scalar.fmah(1.0, -2.0, 3.0));
    try std.testing.expectEqual(@mulAdd(f16, 1.0, 2.0, -3.0), scalar.fmah(1.0, 2.0, -3.0));
    try std.testing.expectEqual(@mulAdd(f16, 0.0, 2.0, 3.0), scalar.fmah(0.0, 2.0, 3.0));
    try std.testing.expectEqual(@mulAdd(f16, 1.0, 0.0, 3.0), scalar.fmah(1.0, 0.0, 3.0));
    try std.testing.expectEqual(@mulAdd(f16, 1.0, 2.0, 0.0), scalar.fmah(1.0, 2.0, 0.0));
}
