const std = @import("std");
const scalar = @import("scalar.zig");

test "fmad NaN handling" {
    const quiet_nan_payload: f64 = @bitCast(@as(u64, 0x7ff8000000000000 | 0x0000000000001234));
    const quiet_nan_neg: f64 = @bitCast(@as(u64, 0xfff8000000000000 | 0x0000000000001234));
    {
        const result = scalar.fmad(quiet_nan_payload, 1.0, 2.0);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000 | 0x0000000000001234), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(1.0, quiet_nan_payload, 2.0);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000 | 0x0000000000001234), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(1.0, 2.0, quiet_nan_payload);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000 | 0x0000000000001234), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(quiet_nan_neg, 1.0, 2.0);
        try std.testing.expectEqual(@as(u64, 0xfff8000000000000 | 0x0000000000001234), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(1.0, quiet_nan_neg, 2.0);
        try std.testing.expectEqual(@as(u64, 0xfff8000000000000 | 0x0000000000001234), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(1.0, 2.0, quiet_nan_neg);
        try std.testing.expectEqual(@as(u64, 0xfff8000000000000 | 0x0000000000001234), @as(u64, @bitCast(result)));
    }
}

test "fmad inf * zero -> nan" {
    const inf: f64 = @bitCast(@as(u64, 0x7ff0000000000000));
    const neg_inf: f64 = @bitCast(@as(u64, 0xfff0000000000000));
    const zero: f64 = 0.0;
    const neg_zero: f64 = @bitCast(@as(u64, 0x8000000000000000));
    {
        const result = scalar.fmad(inf, zero, 1.0);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(neg_inf, zero, 1.0);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(inf, neg_zero, 1.0);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(neg_inf, neg_zero, 1.0);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(zero, inf, 1.0);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(neg_zero, inf, 1.0);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(zero, neg_inf, 1.0);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(neg_zero, neg_inf, 1.0);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000), @as(u64, @bitCast(result)));
    }
}

test "fmad product inf + inf opposite sign -> nan" {
    const inf: f64 = @bitCast(@as(u64, 0x7ff0000000000000));
    const neg_inf: f64 = @bitCast(@as(u64, 0xfff0000000000000));
    {
        const result = scalar.fmad(inf, 1.0, neg_inf);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(neg_inf, 1.0, inf);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(inf, -1.0, inf);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000), @as(u64, @bitCast(result)));
    }
    {
        const result = scalar.fmad(neg_inf, -1.0, neg_inf);
        try std.testing.expectEqual(@as(u64, 0x7ff8000000000000), @as(u64, @bitCast(result)));
    }
}

test "fmad normal values" {
    const a: f64 = 1.0;
    const b: f64 = 2.0;
    const c: f64 = 3.0;
    try std.testing.expectEqual(@mulAdd(f64, a, b, c), scalar.fmad(a, b, c));
    try std.testing.expectEqual(@mulAdd(f64, -1.0, 2.0, 3.0), scalar.fmad(-1.0, 2.0, 3.0));
    try std.testing.expectEqual(@mulAdd(f64, 1.0, -2.0, 3.0), scalar.fmad(1.0, -2.0, 3.0));
    try std.testing.expectEqual(@mulAdd(f64, 1.0, 2.0, -3.0), scalar.fmad(1.0, 2.0, -3.0));
    try std.testing.expectEqual(@mulAdd(f64, 0.0, 2.0, 3.0), scalar.fmad(0.0, 2.0, 3.0));
    try std.testing.expectEqual(@mulAdd(f64, 1.0, 0.0, 3.0), scalar.fmad(1.0, 0.0, 3.0));
    try std.testing.expectEqual(@mulAdd(f64, 1.0, 2.0, 0.0), scalar.fmad(1.0, 2.0, 0.0));
}
