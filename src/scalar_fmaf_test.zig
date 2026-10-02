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

test "fmaf keeps tiny perturbations on both sides of a normal midpoint" {
    const cases = [_]struct { a: u32, b: u32, c: u32, expected: u32 }{
        .{ .a = 0x3f80_0800, .b = 0x3f80_0800, .c = 0x0d80_0000, .expected = 0x3f80_1001 },
        .{ .a = 0x3f80_0800, .b = 0x3f80_1800, .c = 0x8d80_0000, .expected = 0x3f80_2001 },
        .{ .a = 0xbf80_0800, .b = 0x3f80_0800, .c = 0x8d80_0000, .expected = 0xbf80_1001 },
        .{ .a = 0xbf80_0800, .b = 0x3f80_1800, .c = 0x0d80_0000, .expected = 0xbf80_2001 },
    };
    for (cases) |entry| {
        var inputs = [_]f32{ @bitCast(entry.a), @bitCast(entry.b), @bitCast(entry.c) };
        std.mem.doNotOptimizeAway(&inputs);
        const actual: u32 = @bitCast(scalar.fmaf(inputs[0], inputs[1], inputs[2]));
        try std.testing.expectEqual(entry.expected, actual);
    }
}

test "fmaf rounds only after underflow and overflow cancellation" {
    const cases = [_]struct { a: u32, b: u32, c: u32, expected: u32 }{
        .{ .a = 0x0000_0001, .b = 0x3f00_0000, .c = 0x0000_0001, .expected = 0x0000_0002 },
        .{ .a = 0x8000_0001, .b = 0x3f00_0000, .c = 0x8000_0001, .expected = 0x8000_0002 },
        .{ .a = 0x7f7f_ffff, .b = 0x4000_0000, .c = 0xff7f_ffff, .expected = 0x7f7f_ffff },
        .{ .a = 0xff7f_ffff, .b = 0x4000_0000, .c = 0x7f7f_ffff, .expected = 0xff7f_ffff },
    };
    for (cases) |entry| {
        var inputs = [_]f32{ @bitCast(entry.a), @bitCast(entry.b), @bitCast(entry.c) };
        std.mem.doNotOptimizeAway(&inputs);
        const actual: u32 = @bitCast(scalar.fmaf(inputs[0], inputs[1], inputs[2]));
        try std.testing.expectEqual(entry.expected, actual);
    }
}

fn expectFmafBits(a: u32, b: u32, c: u32, expected: u32) !void {
    var inputs = [_]f32{ @bitCast(a), @bitCast(b), @bitCast(c) };
    std.mem.doNotOptimizeAway(&inputs);
    const actual: u32 = @bitCast(scalar.fmaf(inputs[0], inputs[1], inputs[2]));
    try std.testing.expectEqual(expected, actual);
}

test "fmaf retains the low product bits through cancellation" {
    // max * (1 - 2^-24) = 2^128 - 2^105 + 2^80. Subtracting max leaves
    // -(2^104 - 2^80), not the two-rounding answer -2^104.
    try expectFmafBits(0x7f7f_ffff, 0x3f7f_ffff, 0xff7f_ffff, 0xf37f_ffff);
    try expectFmafBits(0xff7f_ffff, 0x3f7f_ffff, 0x7f7f_ffff, 0x737f_ffff);
}

test "fmaf keeps a subnormal sum whose product lands on a tie" {
    // The exact products are (2^23 - 1) * 2^-150, a midpoint of the subnormal
    // grid that only the fused evaluation resolves against the addend.
    try expectFmafBits(0x3f00_0000, 0x007f_ffff, 0x0000_0001, 0x0040_0000);
    try expectFmafBits(0xbf00_0000, 0x807f_ffff, 0x8000_0001, 0x003f_fffe);
}

test "fmaf rounds once at the subnormal and normal boundary" {
    // 0.5 * (2^-125 - 2^-149) is the midpoint between the largest subnormal
    // and the smallest normal, so it rounds up to the even one.
    try expectFmafBits(0x3f00_0000, 0x00ff_ffff, 0x0000_0000, 0x0080_0000);
    try expectFmafBits(0xbf00_0000, 0x00ff_ffff, 0x0000_0000, 0x8080_0000);
    // One subnormal step below the midpoint is exact.
    try expectFmafBits(0x3f00_0000, 0x00ff_fffe, 0x0000_0000, 0x007f_ffff);
    // One subnormal step below the midpoint plus one more step rounds back down.
    try expectFmafBits(0x3f00_0000, 0x00ff_ffff, 0x8000_0001, 0x007f_fffe);
}

test "fmaf rounds a midpoint from the near side as well" {
    // Same exact midpoints as the normal-range test, with the perturbation on
    // the other side of the tie.
    try expectFmafBits(0x3f80_0800, 0x3f80_0800, 0x2d80_0000, 0x3f80_1001);
    try expectFmafBits(0x3f80_0800, 0x3f80_0800, 0x8d80_0000, 0x3f80_1000);
    try expectFmafBits(0xbf80_0800, 0x3f80_0800, 0x0d80_0000, 0xbf80_1000);
    try expectFmafBits(0x3f80_0800, 0x3f80_1800, 0x0d80_0000, 0x3f80_2002);
}

test "fmaf propagates a non-invalid infinity" {
    try expectFmafBits(0x7f80_0000, 0x3f80_0000, 0x3f80_0000, 0x7f80_0000);
    try expectFmafBits(0x7f80_0000, 0x3f80_0000, 0xbf80_0000, 0x7f80_0000);
    try expectFmafBits(0xff80_0000, 0x3f80_0000, 0x0000_0000, 0xff80_0000);
    try expectFmafBits(0x3f80_0000, 0x3f80_0000, 0x7f80_0000, 0x7f80_0000);
    try expectFmafBits(0x3f80_0000, 0x3f80_0000, 0xff80_0000, 0xff80_0000);
}

test "fmaf zero product keeps the addend zero sign" {
    try expectFmafBits(0x3f80_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000);
    try expectFmafBits(0xbf80_0000, 0x8000_0000, 0x0000_0000, 0x0000_0000);
    try expectFmafBits(0xbf80_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000);
    try expectFmafBits(0xbf80_0000, 0x0000_0000, 0x8000_0000, 0x8000_0000);
    try expectFmafBits(0x3f80_0000, 0x8000_0000, 0x8000_0000, 0x8000_0000);
    try expectFmafBits(0xbf80_0000, 0x0000_0000, 0xbf80_0000, 0xbf80_0000);
    // A zero product keeps a nonzero addend untouched.
    // Two zeros of opposite signs are +0 under round-to-nearest-even.
    try expectFmafBits(0x3f80_0000, 0x0000_0000, 0x8000_0000, 0x0000_0000);
    try expectFmafBits(0x3f80_0000, 0x8000_0000, 0x0000_0000, 0x0000_0000);
}

test "fmaf exact cancellation of nonzero terms is positive zero" {
    try expectFmafBits(0x3f00_0000, 0x3f80_0000, 0xbf00_0000, 0x0000_0000);
    try expectFmafBits(0x3f80_0000, 0x4000_0000, 0xc000_0000, 0x0000_0000);
    try expectFmafBits(0xbf00_0000, 0x3f80_0000, 0x3f00_0000, 0x0000_0000);
}
