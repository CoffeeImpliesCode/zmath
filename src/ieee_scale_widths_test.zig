//! The mode-parameterised `scaleB` and `roundIntegral` families at the
//! binary32 and binary64 widths. The binary16 kernels are the reference: the
//! same construction runs at all three widths, and the binary16 result is what
//! a wider format has to reproduce wherever the binary16 rounding was exact.
//! Every expected pattern is derived from a width's own field layout, so no
//! test can pin a pattern that belongs to another width.
const std = @import("std");
const ieee = @import("ieee.zig");
const scaling = @import("scalar_scaling.zig");

const modes = [_]ieee.RoundMode{
    .nearest_even,
    .toward_zero,
    .toward_positive,
    .toward_negative,
};

/// Binds one width so a check is written once and run at every width, in the
/// same shape as the format descriptors of `ieee.zig`. The fraction width and
/// the bias are the only numbers stated here: 23 and 127, 52 and 1023. Every
/// field constant below is derived from them.
const Width = struct {
    const f = struct {
        const Float = f32;
        const Bits = u32;
        const fraction_bits = 23;
        const bias = 127;
    };
    const d = struct {
        const Float = f64;
        const Bits = u64;
        const fraction_bits = 52;
        const bias = 1023;
    };
};

/// The sign bit, which is the top bit of the width.
fn sign(comptime W: type) W.Bits {
    return @as(W.Bits, 1) << (@bitSizeOf(W.Float) - 1);
}

/// The hidden bit, the lowest bit of the exponent field. With a zero exponent
/// field it reads as the smallest normal, which is how a rounded carry out of
/// the fraction reaches the normal range.
fn hidden(comptime W: type) W.Bits {
    return @as(W.Bits, 1) << W.fraction_bits;
}

/// `2^exponent`, whose field is the bias plus that exponent.
fn powerOfTwo(comptime W: type, exponent: comptime_int) W.Bits {
    return @as(W.Bits, W.bias + exponent) << W.fraction_bits;
}

/// `2^exponent` scaled by `1 + 2^-fraction_below`, so one bit below the field
/// gives three and two bits below gives two and a half.
fn withFractionBelow(
    comptime W: type,
    exponent: comptime_int,
    fraction_below: comptime_int,
) W.Bits {
    return powerOfTwo(W, exponent) | (@as(W.Bits, 1) << (W.fraction_bits - fraction_below));
}

/// The whole exponent field, which is the encoding of infinity.
fn infinity(comptime W: type) W.Bits {
    return @as(W.Bits, 2 * W.bias + 1) << W.fraction_bits;
}

/// One field below infinity, which is the largest finite magnitude.
fn maximumFinite(comptime W: type) W.Bits {
    return infinity(W) - 1;
}

/// The smallest normal, which the hidden bit also spells with a zero field.
fn minimumNormal(comptime W: type) W.Bits {
    return hidden(W);
}

/// The largest subnormal, one fraction field below the smallest normal.
fn maximumSubnormal(comptime W: type) W.Bits {
    return hidden(W) - 1;
}

/// The smallest subnormal, which is one quantum above zero.
fn minimumSubnormal(comptime W: type) W.Bits {
    return 1;
}

/// The quiet bit, the top fraction bit, which is the hidden bit halved.
fn quietBit(comptime W: type) W.Bits {
    return hidden(W) >> 1;
}

/// A NaN payload bit below the quiet bit, so a NaN built from it is signaling.
fn nanPayload(comptime W: type) W.Bits {
    return @as(W.Bits, 1) << (W.fraction_bits - 3);
}

/// A quiet NaN with a payload, which both families must return unchanged.
fn quietNan(comptime W: type) W.Bits {
    return infinity(W) | quietBit(W) | nanPayload(W);
}

/// A signaling NaN with the same payload, which must come back quiet with that
/// payload still set.
fn signalingNan(comptime W: type) W.Bits {
    return infinity(W) | nanPayload(W);
}

fn valueOf(comptime W: type, bits: W.Bits) W.Float {
    return @bitCast(bits);
}

fn bitsOf(comptime W: type, value: W.Float) W.Bits {
    return @bitCast(value);
}

// ---------------------------------------------------------------------------
// scaleB
// ---------------------------------------------------------------------------

/// A finite value scaled past the largest finite value saturates at the format
/// boundary, and the sign each mode gives that boundary is its own.
fn expectScaleOverflow(
    comptime W: type,
    scale: fn (ieee.RoundMode, W.Float, i32) W.Float,
) !void {
    const cases = [_]struct {
        mode: ieee.RoundMode,
        positive_bits: W.Bits,
        negative_bits: W.Bits,
    }{
        .{
            .mode = .nearest_even,
            .positive_bits = infinity(W),
            .negative_bits = infinity(W) | sign(W),
        },
        .{
            .mode = .toward_zero,
            .positive_bits = maximumFinite(W),
            .negative_bits = maximumFinite(W) | sign(W),
        },
        .{
            .mode = .toward_positive,
            .positive_bits = infinity(W),
            .negative_bits = maximumFinite(W) | sign(W),
        },
        .{
            .mode = .toward_negative,
            .positive_bits = maximumFinite(W),
            .negative_bits = infinity(W) | sign(W),
        },
    };

    const positive = valueOf(W, maximumFinite(W));
    const negative = valueOf(W, maximumFinite(W) | sign(W));
    for (cases) |case| {
        const positive_result = bitsOf(W, scale(case.mode, positive, 1));
        const negative_result = bitsOf(W, scale(case.mode, negative, 1));
        try std.testing.expectEqual(case.positive_bits, positive_result);
        try std.testing.expectEqual(case.negative_bits, negative_result);
    }
}

/// The smallest subnormal scaled down is exactly half a quantum, so the tie
/// keeps a zero under nearest-even and toward-zero, and each directed mode
/// reaches one quantum only for the sign it rounds toward.
fn expectScaleUnderflowTie(
    comptime W: type,
    scale: fn (ieee.RoundMode, W.Float, i32) W.Float,
) !void {
    const cases = [_]struct {
        mode: ieee.RoundMode,
        positive_bits: W.Bits,
        negative_bits: W.Bits,
    }{
        .{
            .mode = .nearest_even,
            .positive_bits = 0,
            .negative_bits = sign(W),
        },
        .{
            .mode = .toward_zero,
            .positive_bits = 0,
            .negative_bits = sign(W),
        },
        .{
            .mode = .toward_positive,
            .positive_bits = minimumSubnormal(W),
            .negative_bits = sign(W),
        },
        .{
            .mode = .toward_negative,
            .positive_bits = 0,
            .negative_bits = minimumSubnormal(W) | sign(W),
        },
    };

    const positive = valueOf(W, minimumSubnormal(W));
    const negative = valueOf(W, minimumSubnormal(W) | sign(W));
    for (cases) |case| {
        const positive_result = bitsOf(W, scale(case.mode, positive, -1));
        const negative_result = bitsOf(W, scale(case.mode, negative, -1));
        try std.testing.expectEqual(case.positive_bits, positive_result);
        try std.testing.expectEqual(case.negative_bits, negative_result);
    }
}

/// Gradual underflow and exact promotion hold in every mode, because a scale by
/// a power of two needs no rounding while the significand is unchanged.
fn expectScaleExact(
    comptime W: type,
    scale: fn (ieee.RoundMode, W.Float, i32) W.Float,
) !void {
    const fraction_bits: i32 = W.fraction_bits;
    const tiny_bits = minimumSubnormal(W);
    const half_normal_bits = minimumNormal(W) >> 1;
    const widest_bits = maximumSubnormal(W);
    const doubled_bits = minimumNormal(W) | (maximumSubnormal(W) - 1);
    const tiny = valueOf(W, tiny_bits);
    const negative_tiny = valueOf(W, tiny_bits | sign(W));
    const half_normal = valueOf(W, half_normal_bits);
    const widest = valueOf(W, widest_bits);

    for (modes) |mode| {
        // Every exponent from one to the fraction width is an exact step at
        // one, and the top of that range is the largest integer step the
        // format holds exactly.
        const one = valueOf(W, powerOfTwo(W, 0));
        inline for (1..W.fraction_bits + 1) |step| {
            const expected = powerOfTwo(W, step);
            const actual = bitsOf(W, scale(mode, one, step));
            try std.testing.expectEqual(expected, actual);
        }

        // The smallest subnormal reaches the smallest normal after exactly one
        // fraction width of scaling, half the smallest normal doubles into it,
        // and the largest subnormal doubles into the normal range with no bit
        // lost.
        const promoted = bitsOf(W, scale(mode, tiny, fraction_bits));
        const negative_promoted = bitsOf(W, scale(mode, negative_tiny, fraction_bits));
        const doubled_half = bitsOf(W, scale(mode, half_normal, 1));
        const doubled_widest = bitsOf(W, scale(mode, widest, 1));
        const kept_widest = bitsOf(W, scale(mode, widest, 0));
        try std.testing.expectEqual(minimumNormal(W), promoted);
        try std.testing.expectEqual(minimumNormal(W) | sign(W), negative_promoted);
        try std.testing.expectEqual(minimumNormal(W), doubled_half);
        try std.testing.expectEqual(doubled_bits, doubled_widest);
        try std.testing.expectEqual(widest_bits, kept_widest);
    }
}

/// Half a quantum above a subnormal quantum is a tie, and the even quantum
/// wins under nearest-even while each directed mode takes its own side.
fn expectScaleParityTies(
    comptime W: type,
    scale: fn (ieee.RoundMode, W.Float, i32) W.Float,
) !void {
    // The smallest normal plus one quantum of its own significand is half a
    // subnormal quantum above the half-way subnormal, and plus three is half a
    // quantum above that quantum's successor.
    const half_quantum = minimumNormal(W) >> 1;
    const even_input = minimumNormal(W) | 1;
    const odd_input = minimumNormal(W) | 3;
    // A directed mode moves the magnitude of a negative value the other way, so
    // each mode states both signs.
    const cases = [_]struct {
        mode: ieee.RoundMode,
        even_bits: W.Bits,
        negative_even_bits: W.Bits,
        odd_bits: W.Bits,
        negative_odd_bits: W.Bits,
    }{
        .{
            .mode = .nearest_even,
            .even_bits = half_quantum,
            .negative_even_bits = half_quantum,
            .odd_bits = half_quantum + 2,
            .negative_odd_bits = half_quantum + 2,
        },
        .{
            .mode = .toward_zero,
            .even_bits = half_quantum,
            .negative_even_bits = half_quantum,
            .odd_bits = half_quantum + 1,
            .negative_odd_bits = half_quantum + 1,
        },
        .{
            .mode = .toward_positive,
            .even_bits = half_quantum + 1,
            .negative_even_bits = half_quantum,
            .odd_bits = half_quantum + 2,
            .negative_odd_bits = half_quantum + 1,
        },
        .{
            .mode = .toward_negative,
            .even_bits = half_quantum,
            .negative_even_bits = half_quantum + 1,
            .odd_bits = half_quantum + 1,
            .negative_odd_bits = half_quantum + 2,
        },
    };

    for (cases) |case| {
        const even_positive = bitsOf(W, scale(case.mode, valueOf(W, even_input), -1));
        const odd_positive = bitsOf(W, scale(case.mode, valueOf(W, odd_input), -1));
        try std.testing.expectEqual(case.even_bits, even_positive);
        try std.testing.expectEqual(case.odd_bits, odd_positive);
        const even_value = valueOf(W, even_input | sign(W));
        const odd_value = valueOf(W, odd_input | sign(W));
        const even_negative = bitsOf(W, scale(case.mode, even_value, -1));
        const odd_negative = bitsOf(W, scale(case.mode, odd_value, -1));
        try std.testing.expectEqual(case.negative_even_bits | sign(W), even_negative);
        try std.testing.expectEqual(case.negative_odd_bits | sign(W), odd_negative);
    }
}

/// The largest value one field above the smallest normal is half a subnormal
/// quantum above the largest subnormal, so that tie rounds up into the
/// smallest normal itself under every mode that leaves the subnormal range.
/// A directed mode that rounds a negative value toward zero keeps it in the
/// subnormal range instead, so each mode states both signs.
fn expectScaleCarryToNormal(
    comptime W: type,
    scale: fn (ieee.RoundMode, W.Float, i32) W.Float,
) !void {
    const input = minimumNormal(W) | maximumSubnormal(W);
    const cases = [_]struct {
        mode: ieee.RoundMode,
        positive_bits: W.Bits,
        negative_bits: W.Bits,
    }{
        .{
            .mode = .nearest_even,
            .positive_bits = minimumNormal(W),
            .negative_bits = minimumNormal(W),
        },
        .{
            .mode = .toward_zero,
            .positive_bits = maximumSubnormal(W),
            .negative_bits = maximumSubnormal(W),
        },
        .{
            .mode = .toward_positive,
            .positive_bits = minimumNormal(W),
            .negative_bits = maximumSubnormal(W),
        },
        .{
            .mode = .toward_negative,
            .positive_bits = maximumSubnormal(W),
            .negative_bits = minimumNormal(W),
        },
    };

    for (cases) |case| {
        const positive = bitsOf(W, scale(case.mode, valueOf(W, input), -1));
        const negative = bitsOf(W, scale(case.mode, valueOf(W, input | sign(W)), -1));
        try std.testing.expectEqual(case.positive_bits, positive);
        try std.testing.expectEqual(case.negative_bits | sign(W), negative);
    }
}

/// Zero and infinity pass through at any exponent, a NaN keeps its sign and
/// payload and comes back quiet, and the largest i32 exponent in either
/// direction saturates the result the way the same exponent does at binary16.
fn expectScaleSpecials(
    comptime W: type,
    scale: fn (ieee.RoundMode, W.Float, i32) W.Float,
) !void {
    const cases = [_]struct {
        mode: ieee.RoundMode,
        over_bits: W.Bits,
        negative_over_bits: W.Bits,
        under_bits: W.Bits,
        negative_under_bits: W.Bits,
    }{
        .{
            .mode = .nearest_even,
            .over_bits = infinity(W),
            .negative_over_bits = infinity(W) | sign(W),
            .under_bits = 0,
            .negative_under_bits = sign(W),
        },
        .{
            .mode = .toward_zero,
            .over_bits = maximumFinite(W),
            .negative_over_bits = maximumFinite(W) | sign(W),
            .under_bits = 0,
            .negative_under_bits = sign(W),
        },
        .{
            .mode = .toward_positive,
            .over_bits = infinity(W),
            .negative_over_bits = maximumFinite(W) | sign(W),
            .under_bits = minimumSubnormal(W),
            .negative_under_bits = sign(W),
        },
        .{
            .mode = .toward_negative,
            .over_bits = maximumFinite(W),
            .negative_over_bits = infinity(W) | sign(W),
            .under_bits = 0,
            .negative_under_bits = minimumSubnormal(W) | sign(W),
        },
    };

    const tiny = valueOf(W, minimumSubnormal(W));
    const negative_tiny = valueOf(W, minimumSubnormal(W) | sign(W));
    for (cases) |case| {
        const over = bitsOf(W, scale(case.mode, tiny, std.math.maxInt(i32)));
        const negative_over = bitsOf(W, scale(
            case.mode,
            negative_tiny,
            std.math.maxInt(i32),
        ));
        const under = bitsOf(W, scale(case.mode, tiny, std.math.minInt(i32)));
        const negative_under = bitsOf(W, scale(
            case.mode,
            negative_tiny,
            std.math.minInt(i32),
        ));
        try std.testing.expectEqual(case.over_bits, over);
        try std.testing.expectEqual(case.negative_over_bits, negative_over);
        try std.testing.expectEqual(case.under_bits, under);
        try std.testing.expectEqual(case.negative_under_bits, negative_under);
    }

    for (modes) |mode| {
        const positive_zero = bitsOf(W, scale(mode, valueOf(W, 0), std.math.minInt(i32)));
        const negative_zero = bitsOf(W, scale(
            mode,
            valueOf(W, sign(W)),
            std.math.maxInt(i32),
        ));
        try std.testing.expectEqual(@as(W.Bits, 0), positive_zero);
        try std.testing.expectEqual(sign(W), negative_zero);
        const positive_infinity = bitsOf(W, scale(
            mode,
            valueOf(W, infinity(W)),
            std.math.minInt(i32),
        ));
        const negative_infinity = bitsOf(W, scale(
            mode,
            valueOf(W, infinity(W) | sign(W)),
            std.math.maxInt(i32),
        ));
        try std.testing.expectEqual(infinity(W), positive_infinity);
        try std.testing.expectEqual(infinity(W) | sign(W), negative_infinity);
        const quiet = bitsOf(W, scale(
            mode,
            valueOf(W, quietNan(W)),
            std.math.maxInt(i32),
        ));
        const signaling = bitsOf(W, scale(
            mode,
            valueOf(W, signalingNan(W)),
            std.math.minInt(i32),
        ));
        const negative_signaling = bitsOf(W, scale(
            mode,
            valueOf(W, signalingNan(W) | sign(W)),
            std.math.maxInt(i32),
        ));
        try std.testing.expectEqual(quietNan(W), quiet);
        try std.testing.expectEqual(quietNan(W), signaling);
        try std.testing.expectEqual(quietNan(W) | sign(W), negative_signaling);
    }
}

/// Every exact domain endpoint is a fixed point of a scale by no exponent, and
/// the largest finite value is one too.
fn expectScaleEndpoints(
    comptime W: type,
    scale: fn (ieee.RoundMode, W.Float, i32) W.Float,
) !void {
    const endpoints = [_]W.Bits{
        maximumFinite(W),
        maximumFinite(W) | sign(W),
        minimumNormal(W),
        minimumNormal(W) | sign(W),
        maximumSubnormal(W),
        maximumSubnormal(W) | sign(W),
        minimumSubnormal(W),
        minimumSubnormal(W) | sign(W),
    };

    for (modes) |mode| {
        for (endpoints) |bits| {
            const value = valueOf(W, bits);
            try std.testing.expectEqual(bits, bitsOf(W, scale(mode, value, 0)));
        }
    }
}

/// A binary16 encoding widened to a wider format, which is exact in both
/// directions. The fraction field is 13 bits wider at binary32 and 42 bits
/// wider at binary64, so a normal significand shifts up with it; a binary16
/// subnormal is normalized first, because its quantum is 2^-24 and the
/// quantum of the wider subnormal range is not; and the reserved field stays
/// reserved.
fn widenH(comptime W: type, h_bits: u16) W.Bits {
    const h_sign: u16 = 0x8000;
    const h_exponent_mask: u16 = 0x7c00;
    const h_fraction_mask: u16 = 0x03ff;
    const h_bias: i32 = 15;
    const h_reserved: u16 = 31;
    const h_fraction_width = 10;

    const sign_bits: W.Bits = @as(W.Bits, h_bits & h_sign) << (@bitSizeOf(W.Float) - 16);
    const magnitude: u16 = h_bits & ~h_sign;
    const field: u16 = (magnitude & h_exponent_mask) >> 10;
    const fraction: u16 = magnitude & h_fraction_mask;
    // A subnormal is normalized by the same count the binary16 kernels use:
    // the distance from the top of its fraction to the hidden bit, which is
    // 1..10 for a nonzero subnormal and 0 for a zero.
    const leading: u4 = if (field == 0 and fraction != 0)
        @intCast(@clz(fraction) - @clz(@as(u16, 0x0400)))
    else
        0;
    const significand: W.Bits = @as(W.Bits, fraction) << leading;
    // The fraction field is 13 bits wider at binary32 and 42 bits wider at
    // binary64, so a fraction shifts up by that much. A subnormal has the same
    // shift once it is normalized, but its restored hidden bit has to come
    // off first, because in a normal it is carried by the field instead.
    const fraction_shift = @as(std.math.Log2Int(W.Bits), W.fraction_bits - h_fraction_width);
    const widened_fraction: W.Bits = if (field == 0)
        (significand & @as(W.Bits, h_fraction_mask)) << fraction_shift
    else
        significand << fraction_shift;
    // A normal field is in 1..30, so its unbiased exponent is in -14..15, and
    // a normalized subnormal exponent is in -24..-15. The widened field is
    // then in bias - 24 .. bias + 15, which is 103..142 at binary32 and
    // 999..1038 at binary64: positive and far below the reserved field, so
    // the cast holds.
    const unbiased: i32 = if (field == 0)
        -14 - @as(i32, leading)
    else
        @as(i32, field) - h_bias;
    const widened_field: W.Bits = if (field == h_reserved)
        @as(W.Bits, 2 * W.bias + 1)
    else if (field == 0 and fraction == 0)
        0
    else
        @intCast(unbiased + W.bias);
    return sign_bits | (widened_field << W.fraction_bits) | widened_fraction;
}

/// Scaling a widened binary16 value by a power of two is the same exact
/// operation, so the wider result is the widened binary16 result whenever the
/// binary16 rounding was a no-op. This is the cross-width form of the mode
/// matrix: the binary16 behaviour decides what the binary32 and binary64
/// answers have to be.
///
/// A normal binary16 result does not imply that no-op. The two widths round
/// the same exact value, and they part company exactly when binary16 dropped a
/// bit: a binary16 subnormal quantum is a whole number of binary16 significand
/// bits, so a discarded bit leaves a binary16 result that no wider format
/// reproduces. Scaling the binary16 result back by the same exponent returns
/// the input exactly in the no-op case and misses it in the rounding case, and
/// that is the test this skips on.
fn expectScaleAgreesWithH(
    comptime W: type,
    scale: fn (ieee.RoundMode, W.Float, i32) W.Float,
) !void {
    const exponents = [_]i32{ -8, -1, 0, 1, 8 };

    var raw: u32 = 0;
    while (raw < 0x1_0000) : (raw += 1) {
        const h_bits: u16 = @intCast(raw);
        for (modes) |mode| {
            for (exponents) |exponent| {
                const h_result = ieee.scaleBh(mode, ieee.fromBitsh(h_bits), exponent);
                if (ieee.classifyh(h_result) != .normal) continue;
                const undone = ieee.scaleBh(mode, h_result, -exponent);
                if (ieee.toBitsh(undone) != h_bits) continue;
                const widened = valueOf(W, widenH(W, h_bits));
                const expected = widenH(W, ieee.toBitsh(h_result));
                const actual = bitsOf(W, scale(mode, widened, exponent));
                try std.testing.expectEqual(expected, actual);
            }
        }
    }
}

/// Round-to-nearest is the one mode with a second, independently written
/// implementation in this tree: `scalbn` carries the same contract under a
/// different exponent convention and a different normalisation. Every sampled
/// encoding has to agree with it bit for bit, including the overflow and
/// underflow edges where that kernel returns infinity and zero.
fn expectNearestEvenAgreesWithScalbn(
    comptime W: type,
    scale: fn (ieee.RoundMode, W.Float, i32) W.Float,
    reference: fn (W.Float, i32) W.Float,
) !void {
    const exponents = [_]i32{
        -2_000_000_000,
        -1000,
        -W.fraction_bits,
        -1,
        0,
        1,
        W.fraction_bits,
        1000,
        2_000_000_000,
    };

    // A dense sweep over every exponent field including the reserved one, both
    // signs, and the significand endpoints and tie points.
    const fractions = [_]W.Bits{
        0,
        1,
        2,
        hidden(W) >> 1,
        (hidden(W) >> 1) - 1,
        (hidden(W) >> 1) + 1,
        maximumSubnormal(W),
    };
    var field: u64 = 0;
    while (field <= 2 * @as(u64, W.bias) + 1) : (field += 1) {
        const field_bits: W.Bits = @intCast(field);
        for (fractions) |fraction| {
            for ([_]W.Bits{ 0, sign(W) }) |sign_bits| {
                const bits = (field_bits << W.fraction_bits) | fraction | sign_bits;
                const value = valueOf(W, bits);
                for (exponents) |exponent| {
                    const expected = bitsOf(W, reference(value, exponent));
                    const actual = bitsOf(W, scale(.nearest_even, value, exponent));
                    try std.testing.expectEqual(expected, actual);
                }
            }
        }
    }

    // A full-width pseudorandom sweep, so the significand positions the dense
    // sweep skips are covered too.
    var state: u64 = 0x243f_6a88_85a3_08d3;
    var sample: u64 = 0;
    while (sample < 32_768) : (sample += 1) {
        state = state *% 6364136223846793005 +% 1442695040888963407;
        const value = valueOf(W, @as(W.Bits, @truncate(state)));
        for (exponents) |exponent| {
            const expected = bitsOf(W, reference(value, exponent));
            const actual = bitsOf(W, scale(.nearest_even, value, exponent));
            try std.testing.expectEqual(expected, actual);
        }
    }
}

// ---------------------------------------------------------------------------
// roundIntegral
// ---------------------------------------------------------------------------

/// Two and a half is a tie, two is its even neighbour, and each directed mode
/// takes the neighbour on its own side.
fn expectRoundIntegralTies(
    comptime W: type,
    round: fn (ieee.RoundMode, W.Float) W.Float,
) !void {
    const cases = [_]struct {
        mode: ieee.RoundMode,
        positive_bits: W.Bits,
        negative_bits: W.Bits,
    }{
        .{
            .mode = .nearest_even,
            .positive_bits = powerOfTwo(W, 1),
            .negative_bits = powerOfTwo(W, 1) | sign(W),
        },
        .{
            .mode = .toward_zero,
            .positive_bits = powerOfTwo(W, 1),
            .negative_bits = powerOfTwo(W, 1) | sign(W),
        },
        .{
            .mode = .toward_positive,
            .positive_bits = withFractionBelow(W, 1, 1),
            .negative_bits = powerOfTwo(W, 1) | sign(W),
        },
        .{
            .mode = .toward_negative,
            .positive_bits = powerOfTwo(W, 1),
            .negative_bits = withFractionBelow(W, 1, 1) | sign(W),
        },
    };

    const tie = valueOf(W, withFractionBelow(W, 1, 2));
    const negative_tie = valueOf(W, withFractionBelow(W, 1, 2) | sign(W));
    for (cases) |case| {
        const positive = bitsOf(W, round(case.mode, tie));
        const negative = bitsOf(W, round(case.mode, negative_tie));
        try std.testing.expectEqual(case.positive_bits, positive);
        try std.testing.expectEqual(case.negative_bits, negative);
    }
}

/// Three quarters is the smallest value above the half-way point, so it rounds
/// up under the two modes that may and to a signed zero under the two that may
/// not, and the sign of the input is the sign of that zero.
fn expectRoundIntegralBelowOne(
    comptime W: type,
    round: fn (ieee.RoundMode, W.Float) W.Float,
) !void {
    const cases = [_]struct {
        mode: ieee.RoundMode,
        positive_bits: W.Bits,
        negative_bits: W.Bits,
    }{
        .{
            .mode = .nearest_even,
            .positive_bits = powerOfTwo(W, 0),
            .negative_bits = powerOfTwo(W, 0) | sign(W),
        },
        .{
            .mode = .toward_zero,
            .positive_bits = 0,
            .negative_bits = sign(W),
        },
        .{
            .mode = .toward_positive,
            .positive_bits = powerOfTwo(W, 0),
            .negative_bits = sign(W),
        },
        .{
            .mode = .toward_negative,
            .positive_bits = 0,
            .negative_bits = powerOfTwo(W, 0) | sign(W),
        },
    };

    const above_half = valueOf(W, withFractionBelow(W, -1, 1));
    const negative_above_half = valueOf(W, withFractionBelow(W, -1, 1) | sign(W));
    for (cases) |case| {
        const positive = bitsOf(W, round(case.mode, above_half));
        const negative = bitsOf(W, round(case.mode, negative_above_half));
        try std.testing.expectEqual(case.positive_bits, positive);
        try std.testing.expectEqual(case.negative_bits, negative);
    }
}

/// Half is the tie that keeps a signed zero under nearest-even, and the whole
/// subnormal range sits below the half-way point, so the subnormal endpoints
/// and the smallest normal take the same path as three quarters.
fn expectRoundIntegralBelowOneEndpoints(
    comptime W: type,
    round: fn (ieee.RoundMode, W.Float) W.Float,
) !void {
    const cases = [_]struct {
        mode: ieee.RoundMode,
        positive_bits: W.Bits,
        negative_bits: W.Bits,
    }{
        .{
            .mode = .nearest_even,
            .positive_bits = 0,
            .negative_bits = sign(W),
        },
        .{
            .mode = .toward_zero,
            .positive_bits = 0,
            .negative_bits = sign(W),
        },
        .{
            .mode = .toward_positive,
            .positive_bits = powerOfTwo(W, 0),
            .negative_bits = sign(W),
        },
        .{
            .mode = .toward_negative,
            .positive_bits = 0,
            .negative_bits = powerOfTwo(W, 0) | sign(W),
        },
    };

    // The smallest normal, the largest subnormal, and the smallest subnormal,
    // each in both signs.
    const endpoints = [_]W.Bits{
        minimumNormal(W),
        maximumSubnormal(W),
        minimumSubnormal(W),
    };
    const half_bits = powerOfTwo(W, -1);
    for (cases) |case| {
        for (endpoints) |bits| {
            const positive = bitsOf(W, round(case.mode, valueOf(W, bits)));
            const negative = bitsOf(W, round(case.mode, valueOf(W, bits | sign(W))));
            try std.testing.expectEqual(case.positive_bits, positive);
            try std.testing.expectEqual(case.negative_bits, negative);
        }
        const half = bitsOf(W, round(case.mode, valueOf(W, half_bits)));
        const negative_half = bitsOf(W, round(case.mode, valueOf(W, half_bits | sign(W))));
        try std.testing.expectEqual(case.positive_bits, half);
        try std.testing.expectEqual(case.negative_bits, negative_half);
    }
}

/// A value at or above one ulp of the bias exponent is a multiple of one
/// already, and that starts exactly at the fraction width, where both the
/// smallest such integer and a multiple of it have to come back untouched.
fn expectRoundIntegralAlreadyIntegral(
    comptime W: type,
    round: fn (ieee.RoundMode, W.Float) W.Float,
) !void {
    const integral = [_]W.Bits{
        powerOfTwo(W, 0),
        powerOfTwo(W, 1),
        powerOfTwo(W, W.fraction_bits),
        withFractionBelow(W, W.fraction_bits, 1),
        powerOfTwo(W, W.fraction_bits + 1),
        powerOfTwo(W, 100),
        maximumFinite(W),
        maximumFinite(W) | sign(W),
    };

    for (modes) |mode| {
        for (integral) |bits| {
            const value = valueOf(W, bits);
            try std.testing.expectEqual(bits, bitsOf(W, round(mode, value)));
        }
    }
}

/// Zero and infinity pass through, a NaN comes back quiet with its sign and
/// payload, and the negative zero of the input is the negative zero of the
/// result.
fn expectRoundIntegralSpecials(
    comptime W: type,
    round: fn (ieee.RoundMode, W.Float) W.Float,
) !void {
    const special = [_]struct { bits: W.Bits, expected: W.Bits }{
        .{ .bits = 0, .expected = 0 },
        .{ .bits = sign(W), .expected = sign(W) },
        .{ .bits = infinity(W), .expected = infinity(W) },
        .{ .bits = infinity(W) | sign(W), .expected = infinity(W) | sign(W) },
        .{ .bits = quietNan(W), .expected = quietNan(W) },
        .{ .bits = quietNan(W) | sign(W), .expected = quietNan(W) | sign(W) },
        .{ .bits = signalingNan(W), .expected = quietNan(W) },
        .{ .bits = signalingNan(W) | sign(W), .expected = quietNan(W) | sign(W) },
    };

    for (modes) |mode| {
        for (special) |case| {
            const actual = bitsOf(W, round(mode, valueOf(W, case.bits)));
            try std.testing.expectEqual(case.expected, actual);
        }
    }
}

/// Rounding a widened binary16 value to an integer is the same operation, so
/// wherever the binary16 result is finite, the wider result is that result
/// widened.
fn expectRoundIntegralAgreesWithH(
    comptime W: type,
    round: fn (ieee.RoundMode, W.Float) W.Float,
) !void {
    var raw: u32 = 0;
    while (raw < 0x1_0000) : (raw += 1) {
        const h_bits: u16 = @intCast(raw);
        const widened = valueOf(W, widenH(W, h_bits));
        for (modes) |mode| {
            const h_result = ieee.roundIntegralh(mode, ieee.fromBitsh(h_bits));
            const category = ieee.classifyh(h_result);
            if (category == .infinity or
                category == .quiet_nan or
                category == .signaling_nan)
            {
                continue;
            }
            const expected = widenH(W, ieee.toBitsh(h_result));
            const actual = bitsOf(W, round(mode, widened));
            try std.testing.expectEqual(expected, actual);
        }
    }
}

// ---------------------------------------------------------------------------
// The binary16 kernels from before the shared cores
// ---------------------------------------------------------------------------

// Copied from `ieee.zig` as it stood while only the binary16 forms existed,
// with the private mask names spelled out. This is the oracle for the rewrite:
// if the width-generic core moves a single binary16 result, the comparisons at
// the end of this file find it over the whole binary16 encoding space.
const reference_sign_mask: u16 = 0x8000;
const reference_exponent_mask: u16 = 0x7c00;
const reference_fraction_mask: u16 = 0x03ff;
const reference_quiet_mask: u16 = 0x0200;

fn referenceIsNaNh(bits: u16) bool {
    return bits & ~reference_sign_mask > reference_exponent_mask;
}

fn referenceQuietNaNh(bits: u16) f16 {
    return ieee.fromBitsh(bits | reference_quiet_mask);
}

fn referenceOverflowH(mode: ieee.RoundMode, negative: bool) f16 {
    const magnitude: u16 = switch (mode) {
        .nearest_even => reference_exponent_mask,
        .toward_zero => reference_exponent_mask - 1,
        .toward_positive => if (negative) reference_exponent_mask - 1 else reference_exponent_mask,
        .toward_negative => if (negative) reference_exponent_mask else reference_exponent_mask - 1,
    };
    return ieee.fromBitsh((if (negative) reference_sign_mask else 0) | magnitude);
}

fn referenceRoundShiftRightH(
    mode: ieee.RoundMode,
    negative: bool,
    significand: u16,
    shift: i64,
) u16 {
    if (shift > 11) {
        return switch (mode) {
            .nearest_even, .toward_zero => 0,
            .toward_positive => if (negative) 0 else 1,
            .toward_negative => if (negative) 1 else 0,
        };
    }

    const shift_amount: u4 = @intCast(shift);
    const quotient = significand >> shift_amount;
    const discarded_mask = (@as(u16, 1) << shift_amount) - 1;
    const discarded = significand & discarded_mask;
    const half_shift: u4 = @intCast(shift - 1);
    const half = @as(u16, 1) << half_shift;
    const increment = switch (mode) {
        .nearest_even => discarded > half or
            (discarded == half and quotient & 1 != 0),
        .toward_zero => false,
        .toward_positive => !negative and discarded != 0,
        .toward_negative => negative and discarded != 0,
    };
    return if (increment) quotient + 1 else quotient;
}

// Copied from `ieee.zig` as it stood while only the binary16 forms existed,
// with the private mask names spelled out and the local sign variable renamed
// so it does not shadow the width helper above. Nothing else is changed, so
// this is the oracle for the rewrite: if the width-generic core moves a single
// binary16 result, the comparisons at the end of this file find it over the
// whole binary16 encoding space.
fn referenceScaleBh(mode: ieee.RoundMode, value: f16, exponent: i32) f16 {
    const bits = ieee.toBitsh(value);
    const sign_bits = bits & reference_sign_mask;
    const magnitude = bits & ~reference_sign_mask;

    if (referenceIsNaNh(bits)) return referenceQuietNaNh(bits);
    if (magnitude == 0 or magnitude == reference_exponent_mask) return value;

    var significand = (magnitude & reference_fraction_mask) | 0x0400;
    var exponent_field: i64 = @intCast((magnitude & reference_exponent_mask) >> 10);
    if (exponent_field == 0) {
        significand = magnitude & reference_fraction_mask;
        exponent_field = 1;
        const normalization_shift: u4 =
            @intCast(@clz(significand) - @clz(@as(u16, 0x0400)));
        significand <<= normalization_shift;
        exponent_field -= normalization_shift;
    }

    const scaled_exponent = exponent_field + @as(i64, exponent);
    if (scaled_exponent > 30) {
        return referenceOverflowH(mode, sign_bits != 0);
    }
    if (scaled_exponent >= 1) {
        const result_exponent: u16 = @intCast(scaled_exponent);
        const result_bits = sign_bits |
            (result_exponent << 10) |
            (significand & reference_fraction_mask);
        return ieee.fromBitsh(result_bits);
    }

    const shift = 1 - scaled_exponent;
    const rounded_significand = referenceRoundShiftRightH(
        mode,
        sign_bits != 0,
        significand,
        shift,
    );
    return ieee.fromBitsh(sign_bits | rounded_significand);
}

fn referenceRoundIntegralh(mode: ieee.RoundMode, value: f16) f16 {
    const bits = ieee.toBitsh(value);
    const sign_bits = bits & reference_sign_mask;
    const magnitude = bits & ~reference_sign_mask;

    if (referenceIsNaNh(bits)) return referenceQuietNaNh(bits);
    if (magnitude == 0 or magnitude == reference_exponent_mask) return value;

    const exponent_field: i32 = @intCast((magnitude & reference_exponent_mask) >> 10);
    const unbiased_exponent = exponent_field - 15;
    if (unbiased_exponent >= 10) return value;

    if (unbiased_exponent < 0) {
        const rounded_magnitude: u16 = switch (mode) {
            .nearest_even => if (unbiased_exponent == -1 and
                magnitude & reference_fraction_mask != 0)
                0x3c00
            else
                0,
            .toward_zero => 0,
            .toward_positive => if (sign_bits == 0) 0x3c00 else 0,
            .toward_negative => if (sign_bits != 0) 0x3c00 else 0,
        };
        return ieee.fromBitsh(sign_bits | rounded_magnitude);
    }

    const fractional_bits: u4 = @intCast(10 - unbiased_exponent);
    const discarded_mask = (@as(u16, 1) << fractional_bits) - 1;
    const discarded = magnitude & discarded_mask;
    if (discarded == 0) return value;

    const truncated = magnitude & ~discarded_mask;
    const unit = discarded_mask + 1;
    const half = unit >> 1;
    const increment = switch (mode) {
        .nearest_even => discarded > half or
            (discarded == half and (truncated >> fractional_bits) & 1 != 0),
        .toward_zero => false,
        .toward_positive => sign_bits == 0,
        .toward_negative => sign_bits != 0,
    };
    const rounded_magnitude = truncated + if (increment) unit else 0;
    return ieee.fromBitsh(sign_bits | rounded_magnitude);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "the derived binary32 constants match the special-value table" {
    try std.testing.expectEqual(@as(u32, 0x0000_0000), bitsOf(Width.f, 0.0));
    try std.testing.expectEqual(@as(u32, 0x8000_0000), bitsOf(Width.f, -0.0));
    try std.testing.expectEqual(@as(u32, 0x3f80_0000), powerOfTwo(Width.f, 0));
    try std.testing.expectEqual(@as(u32, 0x7f80_0000), infinity(Width.f));
    try std.testing.expectEqual(@as(u32, 0xff80_0000), infinity(Width.f) | sign(Width.f));
    try std.testing.expectEqual(@as(u32, 0x7fc0_0000), infinity(Width.f) | hidden(Width.f) >> 1);
    try std.testing.expectEqual(@as(u32, 0x7f7f_ffff), maximumFinite(Width.f));
    try std.testing.expectEqual(@as(u32, 0x0080_0000), minimumNormal(Width.f));
    try std.testing.expectEqual(@as(u32, 0x007f_ffff), maximumSubnormal(Width.f));
    try std.testing.expectEqual(@as(u32, 0x0000_0001), minimumSubnormal(Width.f));
    try std.testing.expectEqual(@as(f32, 1.0), valueOf(Width.f, powerOfTwo(Width.f, 0)));
    try std.testing.expectEqual(@as(f32, 0x1.fffffep127), valueOf(Width.f, maximumFinite(Width.f)));
}

test "the derived binary64 constants match the special-value table" {
    try std.testing.expectEqual(@as(u64, 0x0000_0000_0000_0000), bitsOf(Width.d, 0.0));
    try std.testing.expectEqual(@as(u64, 0x8000_0000_0000_0000), bitsOf(Width.d, -0.0));
    try std.testing.expectEqual(@as(u64, 0x3ff0_0000_0000_0000), powerOfTwo(Width.d, 0));
    try std.testing.expectEqual(@as(u64, 0x7ff0_0000_0000_0000), infinity(Width.d));
    try std.testing.expectEqual(@as(u64, 0xfff0_0000_0000_0000), infinity(Width.d) | sign(Width.d));
    try std.testing.expectEqual(
        @as(u64, 0x7ff8_0000_0000_0000),
        infinity(Width.d) | hidden(Width.d) >> 1,
    );
    try std.testing.expectEqual(@as(u64, 0x7fef_ffff_ffff_ffff), maximumFinite(Width.d));
    try std.testing.expectEqual(@as(u64, 0x0010_0000_0000_0000), minimumNormal(Width.d));
    try std.testing.expectEqual(@as(u64, 0x000f_ffff_ffff_ffff), maximumSubnormal(Width.d));
    try std.testing.expectEqual(@as(u64, 1), minimumSubnormal(Width.d));
    try std.testing.expectEqual(@as(f64, 1.0), valueOf(Width.d, powerOfTwo(Width.d, 0)));
}

test "scaleBf saturates a positive and a negative overflow in every mode" {
    try expectScaleOverflow(Width.f, ieee.scaleBf);
}

test "scaleBd saturates a positive and a negative overflow in every mode" {
    try expectScaleOverflow(Width.d, ieee.scaleBd);
}

test "scaleBf rounds half a subnormal quantum down in every mode" {
    try expectScaleUnderflowTie(Width.f, ieee.scaleBf);
}

test "scaleBd rounds half a subnormal quantum down in every mode" {
    try expectScaleUnderflowTie(Width.d, ieee.scaleBd);
}

test "scaleBf keeps subnormal and normal steps exact in every mode" {
    try expectScaleExact(Width.f, ieee.scaleBf);
}

test "scaleBd keeps subnormal and normal steps exact in every mode" {
    try expectScaleExact(Width.d, ieee.scaleBd);
}

test "scaleBf rounds subnormal ties by parity in every mode" {
    try expectScaleParityTies(Width.f, ieee.scaleBf);
}

test "scaleBd rounds subnormal ties by parity in every mode" {
    try expectScaleParityTies(Width.d, ieee.scaleBd);
}

test "scaleBf carries a subnormal tie into the minimum normal" {
    try expectScaleCarryToNormal(Width.f, ieee.scaleBf);
}

test "scaleBd carries a subnormal tie into the minimum normal" {
    try expectScaleCarryToNormal(Width.d, ieee.scaleBd);
}

test "scaleBf keeps special values and saturates extreme exponents" {
    try expectScaleSpecials(Width.f, ieee.scaleBf);
}

test "scaleBd keeps special values and saturates extreme exponents" {
    try expectScaleSpecials(Width.d, ieee.scaleBd);
}

test "scaleBf returns every domain endpoint unchanged" {
    try expectScaleEndpoints(Width.f, ieee.scaleBf);
}

test "scaleBd returns every domain endpoint unchanged" {
    try expectScaleEndpoints(Width.d, ieee.scaleBd);
}

test "scaleBf agrees with the binary16 result wherever the binary16 scale is exact" {
    try expectScaleAgreesWithH(Width.f, ieee.scaleBf);
}

test "scaleBd agrees with the binary16 result wherever the binary16 scale is exact" {
    try expectScaleAgreesWithH(Width.d, ieee.scaleBd);
}

test "scaleBf round-to-nearest agrees with the accepted binary32 scalbn" {
    try expectNearestEvenAgreesWithScalbn(Width.f, ieee.scaleBf, scaling.scalbnf);
}

test "scaleBd round-to-nearest agrees with the accepted binary64 scalbn" {
    try expectNearestEvenAgreesWithScalbn(Width.d, ieee.scaleBd, scaling.scalbnd);
}

test "roundIntegralf rounds positive and negative ties in every mode" {
    try expectRoundIntegralTies(Width.f, ieee.roundIntegralf);
}

test "roundIntegrald rounds positive and negative ties in every mode" {
    try expectRoundIntegralTies(Width.d, ieee.roundIntegrald);
}

test "roundIntegralf keeps a signed zero below the half-way point" {
    try expectRoundIntegralBelowOne(Width.f, ieee.roundIntegralf);
    try expectRoundIntegralBelowOneEndpoints(Width.f, ieee.roundIntegralf);
}

test "roundIntegrald keeps a signed zero below the half-way point" {
    try expectRoundIntegralBelowOne(Width.d, ieee.roundIntegrald);
    try expectRoundIntegralBelowOneEndpoints(Width.d, ieee.roundIntegrald);
}

test "roundIntegralf leaves every integral domain endpoint unchanged" {
    try expectRoundIntegralAlreadyIntegral(Width.f, ieee.roundIntegralf);
    try expectRoundIntegralSpecials(Width.f, ieee.roundIntegralf);
}

test "roundIntegrald leaves every integral domain endpoint unchanged" {
    try expectRoundIntegralAlreadyIntegral(Width.d, ieee.roundIntegrald);
    try expectRoundIntegralSpecials(Width.d, ieee.roundIntegrald);
}

test "roundIntegralf agrees with the binary16 result wherever binary16 stays finite" {
    try expectRoundIntegralAgreesWithH(Width.f, ieee.roundIntegralf);
}

test "roundIntegrald agrees with the binary16 result wherever binary16 stays finite" {
    try expectRoundIntegralAgreesWithH(Width.d, ieee.roundIntegrald);
}

test "scaleBh keeps every binary16 result after the shared-core rewrite" {
    const exponents = [_]i32{ -20, -11, -1, 0, 1, 11, 20 };

    var raw: u32 = 0;
    while (raw < 0x1_0000) : (raw += 1) {
        const bits: u16 = @intCast(raw);
        const value = ieee.fromBitsh(bits);
        for (modes) |mode| {
            for (exponents) |exponent| {
                const expected = ieee.toBitsh(referenceScaleBh(mode, value, exponent));
                const actual = ieee.toBitsh(ieee.scaleBh(mode, value, exponent));
                try std.testing.expectEqual(expected, actual);
            }
        }
    }
}

test "roundIntegralh keeps every binary16 result after the shared-core rewrite" {
    var raw: u32 = 0;
    while (raw < 0x1_0000) : (raw += 1) {
        const bits: u16 = @intCast(raw);
        const value = ieee.fromBitsh(bits);
        for (modes) |mode| {
            const expected = ieee.toBitsh(referenceRoundIntegralh(mode, value));
            const actual = ieee.toBitsh(ieee.roundIntegralh(mode, value));
            try std.testing.expectEqual(expected, actual);
        }
    }
}
