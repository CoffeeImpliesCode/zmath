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

// The expectations below are exact IEEE 754 results of `a*b + c` with a single
// rounding, not samples of one implementation. They therefore hold for the
// fused instruction of a native target and for the integer dyadic engine of
// every other one. Each case goes through `doNotOptimizeAway` so that a
// constant folder cannot replace the portable path with a folded literal, and
// each one isolates a consumer-visible boundary: a tie that a discarded
// alignment bit decides, an underflow or overflow that only a single rounding
// survives, a cancellation that keeps the low product bits, and the sign of an
// exact zero.
fn expectFmadBits(a: u64, b: u64, c: u64, expected: u64) !void {
    var inputs = [_]f64{ @bitCast(a), @bitCast(b), @bitCast(c) };
    std.mem.doNotOptimizeAway(&inputs);
    const actual: u64 = @bitCast(scalar.fmad(inputs[0], inputs[1], inputs[2]));
    try std.testing.expectEqual(expected, actual);
}

test "fmad keeps a subnormal perturbation on both sides of a normal midpoint" {
    // (1 + 2^-26) * (1 + 2^-27) = 1 + 3 * 2^-27 + 2^-53, an exact midpoint of
    // the binary64 grid. Rounding the product first would drop the 2^-53 half
    // and the addend with it, while the fused sum keeps the half as the
    // rounding decision and the addend as a sticky bit that moves it.
    try expectFmadBits(0x3ff0_0000_0400_0000, 0x3ff0_0000_0200_0000, 0x0000_0000_0000_0001, 0x3ff0_0000_0600_0001);
    try expectFmadBits(0x3ff0_0000_0400_0000, 0x3ff0_0000_0200_0000, 0x8000_0000_0000_0001, 0x3ff0_0000_0600_0000);
    try expectFmadBits(0xbff0_0000_0400_0000, 0x3ff0_0000_0200_0000, 0x8000_0000_0000_0001, 0xbff0_0000_0600_0001);
    try expectFmadBits(0xbff0_0000_0400_0000, 0x3ff0_0000_0200_0000, 0x0000_0000_0000_0001, 0xbff0_0000_0600_0000);
    // Without an addend the same product is a tie, and the even neighbour wins.
    try expectFmadBits(0x3ff0_0000_0400_0000, 0x3ff0_0000_0200_0000, 0x0000_0000_0000_0000, 0x3ff0_0000_0600_0000);
}

test "fmad rounds a midpoint whose truncated neighbour is odd" {
    // (1 + 2^-26) * (1 + 2^-26 + 2^-27) = 1 + 5 * 2^-27 + 2^-52 + 2^-53, which
    // is half an ulp above an odd significand. Ties-to-even rounds away from
    // the odd neighbour here, and an addend on either side moves the tie.
    try expectFmadBits(0x3ff0_0000_0400_0000, 0x3ff0_0000_0600_0000, 0x0000_0000_0000_0000, 0x3ff0_0000_0a00_0002);
    try expectFmadBits(0x3ff0_0000_0400_0000, 0x3ff0_0000_0600_0000, 0x0000_0000_0000_0001, 0x3ff0_0000_0a00_0002);
    try expectFmadBits(0x3ff0_0000_0400_0000, 0x3ff0_0000_0600_0000, 0x8000_0000_0000_0001, 0x3ff0_0000_0a00_0001);
}

test "fmad rounds only once across the subnormal and the overflow boundary" {
    // Half of the smallest subnormal is not representable, so the addend is
    // what decides the subnormal result.
    try expectFmadBits(0x0000_0000_0000_0001, 0x3fe0_0000_0000_0000, 0x0000_0000_0000_0001, 0x0000_0000_0000_0002);
    try expectFmadBits(0x8000_0000_0000_0001, 0x3fe0_0000_0000_0000, 0x8000_0000_0000_0001, 0x8000_0000_0000_0002);
    // The exact product overflows, yet the sum is finite, which only a fused
    // evaluation returns as the largest finite value.
    try expectFmadBits(0x7fef_ffff_ffff_ffff, 0x4000_0000_0000_0000, 0xffef_ffff_ffff_ffff, 0x7fef_ffff_ffff_ffff);
    try expectFmadBits(0xffef_ffff_ffff_ffff, 0x4000_0000_0000_0000, 0x7fef_ffff_ffff_ffff, 0xffef_ffff_ffff_ffff);
    // An addend too small to reach the product is still below the midpoint of
    // the overflow range, so the result stays at the largest finite value.
    try expectFmadBits(0x7fef_ffff_ffff_ffff, 0x3ff0_0000_0000_0000, 0x0000_0000_0000_0001, 0x7fef_ffff_ffff_ffff);
    try expectFmadBits(0x7fef_ffff_ffff_ffff, 0x3ff0_0000_0000_0000, 0x7c8f_ffff_ffff_ffff, 0x7fef_ffff_ffff_ffff);
    // One half of an ulp higher is the midpoint of the overflow range, and the
    // even neighbour there is the infinity.
    try expectFmadBits(0x7fef_ffff_ffff_ffff, 0x3ff0_0000_0000_0000, 0x7c90_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectFmadBits(0x7fef_ffff_ffff_ffff, 0x4000_0000_0000_0000, 0x0000_0000_0000_0000, 0x7ff0_0000_0000_0000);
    // A product that underflows to nothing still leaves the addend untouched.
    try expectFmadBits(0x0000_0000_0000_0001, 0x0000_0000_0000_0001, 0x8000_0000_0000_0001, 0x8000_0000_0000_0001);
}

test "fmad retains the low product bits through cancellation" {
    // max * (1 - 2^-53) - max is exactly -(2^971 - 2^918). Rounding the product
    // first loses the 2^918 term and returns -2^971 instead.
    try expectFmadBits(0x7fef_ffff_ffff_ffff, 0x3fef_ffff_ffff_ffff, 0xffef_ffff_ffff_ffff, 0xfc9f_ffff_ffff_ffff);
    try expectFmadBits(0xffef_ffff_ffff_ffff, 0x3fef_ffff_ffff_ffff, 0x7fef_ffff_ffff_ffff, 0x7c9f_ffff_ffff_ffff);
    // (1 + 2^-52) * (1 - 2^-52) - 1 is exactly -2^-104, far below the operand
    // that produced it.
    try expectFmadBits(0x3ff0_0000_0000_0001, 0x3fef_ffff_ffff_fffe, 0xbff0_0000_0000_0000, 0xb970_0000_0000_0000);
    try expectFmadBits(0xbff0_0000_0000_0001, 0x3fef_ffff_ffff_fffe, 0x3ff0_0000_0000_0000, 0x3970_0000_0000_0000);
}

test "fmad keeps a subnormal sum whose product lands on a tie" {
    // The exact products are (2^52 - 1) * 2^-1075, a midpoint of the subnormal
    // grid that only the fused evaluation resolves against the addend.
    try expectFmadBits(0x3fe0_0000_0000_0000, 0x000f_ffff_ffff_ffff, 0x0000_0000_0000_0001, 0x0008_0000_0000_0000);
    try expectFmadBits(0xbfe0_0000_0000_0000, 0x800f_ffff_ffff_ffff, 0x8000_0000_0000_0001, 0x0007_ffff_ffff_fffe);
}

test "fmad rounds once at the subnormal and normal boundary" {
    // 0.5 * (2^-1021 - 2^-1074) is the midpoint between the largest subnormal
    // and the smallest normal, so it rounds up to the even one.
    try expectFmadBits(0x3fe0_0000_0000_0000, 0x001f_ffff_ffff_ffff, 0x0000_0000_0000_0000, 0x0010_0000_0000_0000);
    try expectFmadBits(0xbfe0_0000_0000_0000, 0x001f_ffff_ffff_ffff, 0x0000_0000_0000_0000, 0x8010_0000_0000_0000);
    // One subnormal step below the midpoint is exact.
    try expectFmadBits(0x3fe0_0000_0000_0000, 0x001f_ffff_ffff_fffe, 0x0000_0000_0000_0000, 0x000f_ffff_ffff_ffff);
    // One subnormal step below the midpoint plus one more step rounds back down.
    try expectFmadBits(0x3fe0_0000_0000_0000, 0x001f_ffff_ffff_ffff, 0x8000_0000_0000_0001, 0x000f_ffff_ffff_fffe);
}

test "fmad propagates a non-invalid infinity" {
    try expectFmadBits(0x7ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectFmadBits(0x7ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0xbff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectFmadBits(0xfff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectFmadBits(0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectFmadBits(0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0xfff0_0000_0000_0000, 0xfff0_0000_0000_0000);
}

test "fmad zero product keeps the addend zero sign" {
    try expectFmadBits(0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectFmadBits(0xbff0_0000_0000_0000, 0x8000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectFmadBits(0xbff0_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectFmadBits(0xbff0_0000_0000_0000, 0x0000_0000_0000_0000, 0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectFmadBits(0x3ff0_0000_0000_0000, 0x8000_0000_0000_0000, 0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectFmadBits(0xbff0_0000_0000_0000, 0x0000_0000_0000_0000, 0xbff0_0000_0000_0000, 0xbff0_0000_0000_0000);
    // A zero product keeps a nonzero addend untouched.
    // Two zeros of opposite signs are +0 under round-to-nearest-even.
    try expectFmadBits(0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000, 0x8000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectFmadBits(0x3ff0_0000_0000_0000, 0x8000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
}

test "fmad exact cancellation of nonzero terms is positive zero" {
    try expectFmadBits(0x3fe0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0xbfe0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectFmadBits(0x3ff0_0000_0000_0000, 0x4000_0000_0000_0000, 0xc000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectFmadBits(0xbfe0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x3fe0_0000_0000_0000, 0x0000_0000_0000_0000);
    // Exact cancellation of the smallest subnormal still settles on +0.
    try expectFmadBits(0x8000_0000_0000_0001, 0x3ff0_0000_0000_0000, 0x0000_0000_0000_0001, 0x0000_0000_0000_0000);
}
