//! Tests for the `zmath.complex` surface.
//!
//! The expected values here are identities and ranges wherever a bit pattern
//! would have to be guessed. Every pinned constant carries the derivation that
//! produced it, and every magnitude that the naive reduction would overflow is
//! paired with a proof that the naive reduction really does overflow, so a
//! regression cannot make the test pass for the wrong reason.

const std = @import("std");
const testing = std.testing;
const complex = @import("complex.zig");
const scalar = @import("scalar.zig");

const Complexf = complex.Complexf;
const Complexd = complex.Complexd;

const inf_f: f32 = @bitCast(@as(u32, 0x7f80_0000));
const inf_d: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0000));
const neg_inf_f: f32 = @bitCast(@as(u32, 0xff80_0000));
const neg_inf_d: f64 = @bitCast(@as(u64, 0xfff0_0000_0000_0000));

/// The largest finite binary64, whose reciprocal is a subnormal.
const maximum_d: f64 = @bitCast(@as(u64, 0x7fef_ffff_ffff_ffff));

/// The positive canonical quiet NaN of each format, per the special-value
/// table of the specification.
const nan_f: f32 = @bitCast(@as(u32, 0x7fc0_0000));
const nan_d: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_0000));

/// A negative canonical quiet NaN, whose sign and payload must survive.
const neg_nan_f: f32 = @bitCast(@as(u32, 0xffc0_0000));
const neg_nan_d: f64 = @bitCast(@as(u64, 0xfff8_0000_0000_0000));

/// A signaling NaN with a diagnostic payload, which must be quieted in place.
const signaling_nan_f: f32 = @bitCast(@as(u32, 0x7f80_0001));
const signaling_nan_d: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0001));

const pi_f: f32 = @floatCast(std.math.pi);
const pi_d: f64 = std.math.pi;

fn expectBitsF(value: f32, bits: u32) !void {
    try testing.expectEqual(bits, @as(u32, @bitCast(value)));
}

fn expectBitsD(value: f64, bits: u64) !void {
    try testing.expectEqual(bits, @as(u64, @bitCast(value)));
}

/// Asserts that `value` is the positive canonical quiet NaN of `f32`.
fn expectCanonicalF(value: f32) !void {
    try expectBitsF(value, 0x7fc0_0000);
}

/// Asserts that the binary64 squared sum of `x` and `y` overflows, which is
/// the condition that forces the division and magnitude reductions off their
/// direct form. The operands are bound to typed locals because a comptime
/// float would overflow the compiler instead of the format.
fn expectSquaredSumOverflowsD(x: f64, y: f64) !void {
    const sum: f64 = x * x + y * y;
    try testing.expect(std.math.isInf(sum));
}

/// The binary32 counterpart of `expectSquaredSumOverflowsD`.
fn expectSquaredSumOverflowsF(x: f32, y: f32) !void {
    const sum: f32 = x * x + y * y;
    try testing.expect(std.math.isInf(sum));
}

/// Asserts that the binary64 square of `x` overflows to infinity.
fn expectSquareOverflowsD(x: f64) !void {
    const square: f64 = x * x;
    try testing.expect(std.math.isInf(square));
}

/// Asserts that the binary64 square of `x` underflows to zero.
fn expectSquareUnderflowsD(x: f64) !void {
    const square: f64 = x * x;
    try testing.expectEqual(@as(f64, 0), square);
}

/// Asserts that `value` is a binary32 infinity, with a typed local so that the
/// check sees a float.
fn expectOverflowsToInfinityF(value: f32) !void {
    try testing.expect(std.math.isInf(value));
}

/// Asserts that `value` is the positive canonical quiet NaN of `f64`.
fn expectCanonicalD(value: f64) !void {
    try expectBitsD(value, 0x7ff8_0000_0000_0000);
}

/// Asserts that two complex values agree to a relative tolerance in each
/// component. A component that is exactly zero admits any magnitude within the
/// tolerance, because a relative test against zero has no scale.
fn expectComplexApprox(expected: Complexd, actual: Complexd, tolerance: f64) !void {
    if (expected.re == 0) {
        try testing.expectApproxEqAbs(@as(f64, 0), actual.re, tolerance);
    } else {
        try testing.expectApproxEqRel(expected.re, actual.re, tolerance);
    }
    if (expected.im == 0) {
        try testing.expectApproxEqAbs(@as(f64, 0), actual.im, tolerance);
    } else {
        try testing.expectApproxEqRel(expected.im, actual.im, tolerance);
    }
}

/// Asserts that a value is finite, so an infinity or a NaN fails the test.
fn expectFiniteF(value: f32) !void {
    try testing.expect(std.math.isFinite(value));
}

fn expectFiniteD(value: f64) !void {
    try testing.expect(std.math.isFinite(value));
}

test "accessors expose the stored components" {
    const zf: Complexf = .{ .re = 3.5, .im = -4.25 };
    const zd: Complexd = .{ .re = 1e200, .im = -2.5e-200 };
    try testing.expectEqual(zf.re, complex.realf(zf));
    try testing.expectEqual(zf.im, complex.imagf(zf));
    try testing.expectEqual(zd.re, complex.reald(zd));
    try testing.expectEqual(zd.im, complex.imagd(zd));
    // The accessors are the identity, so even a NaN payload must round-trip.
    try expectBitsF(complex.realf(.{ .re = nan_f, .im = 1 }), 0x7fc0_0000);
    try expectBitsD(complex.imagd(.{ .re = 0, .im = signaling_nan_d }), 0x7ff0_0000_0000_0001);
}

test "conj negates only the imaginary part and is an exact involution" {
    try testing.expectEqual(Complexf{ .re = 3, .im = -4 }, complex.conjf(.{ .re = 3, .im = 4 }));
    try testing.expectEqual(Complexd{ .re = 3, .im = -4 }, complex.conjd(.{ .re = 3, .im = 4 }));
    // Negation is exact, so the imaginary part must be the exact negative.
    try expectBitsF(complex.conjf(.{ .re = 1, .im = 0 }).im, 0x8000_0000);
    try expectBitsF(complex.conjf(.{ .re = 1, .im = -0.0 }).im, 0x0000_0000);

    const samples = [_]Complexf{
        .{ .re = 0, .im = 0 },
        .{ .re = -0.0, .im = -0.0 },
        .{ .re = 1e-40, .im = -1e-40 },
        .{ .re = inf_f, .im = -inf_f },
        .{ .re = neg_nan_f, .im = signaling_nan_f },
    };
    for (samples) |z| {
        const twice = complex.conjf(complex.conjf(z));
        try expectBitsF(twice.re, @bitCast(z.re));
        try expectBitsF(twice.im, @bitCast(z.im));
    }
}

test "proj drops the imaginary part and absorbs an infinity" {
    const projected = complex.projf(.{ .re = 3, .im = -4 });
    try testing.expectEqual(@as(f32, 3), projected.re);
    try expectBitsF(projected.im, 0x0000_0000);
    // The real part is stored, so a negative zero keeps its sign.
    try expectBitsD(complex.projd(.{ .re = -0.0, .im = 5 }).re, 0x8000_0000_0000_0000);
    try expectBitsD(complex.projd(.{ .re = -0.0, .im = 5 }).im, 0x0000_0000_0000_0000);

    // An infinite component projects to +inf + i0, where the naive expression
    // inf - i * inf is a NaN.
    for ([_]Complexf{
        .{ .re = inf_f, .im = 1 },
        .{ .re = 1, .im = inf_f },
        .{ .re = inf_f, .im = inf_f },
        .{ .re = inf_f, .im = -inf_f },
        .{ .re = neg_inf_f, .im = 0 },
    }) |z| {
        const p = complex.projf(z);
        try expectBitsF(p.re, 0x7f80_0000);
        try expectBitsF(p.im, 0x0000_0000);
    }
    const p = complex.projd(.{ .re = neg_inf_d, .im = -1 });
    try expectBitsD(p.re, 0x7ff0_0000_0000_0000);
    try expectBitsD(p.im, 0x0000_0000_0000_0000);

    // A NaN has no projection onto a real axis, so it propagates to both parts.
    try expectCanonicalF(complex.projf(.{ .re = nan_f, .im = 1 }).re);
    try expectCanonicalF(complex.projf(.{ .re = 1, .im = nan_f }).im);
    // The negative canonical NaN keeps its sign, so the result is the negative
    // canonical quiet NaN and not the positive one.
    try expectBitsD(complex.projd(.{ .re = 1, .im = neg_nan_d }).im, 0xfff8_0000_0000_0000);
    // A NaN keeps its own sign and payload.
    try expectBitsF(complex.projf(.{ .re = 1, .im = neg_nan_f }).im, 0xffc0_0000);
    try expectBitsF(complex.projf(.{ .re = 1, .im = signaling_nan_f }).im, 0x7fc0_0001);
}

test "add and sub are componentwise" {
    try testing.expectEqual(Complexf{ .re = 4, .im = 6 }, complex.addf(
        .{ .re = 1, .im = 2 },
        .{ .re = 3, .im = 4 },
    ));
    try testing.expectEqual(Complexd{ .re = -1, .im = 0 }, complex.subd(
        .{ .re = 2, .im = 3 },
        .{ .re = 3, .im = 3 },
    ));
    // The IEEE identity must survive the sign of a zero exactly.
    try expectBitsF(complex.addf(.{ .re = -0.0, .im = -0.0 }, .{ .re = -0.0, .im = -0.0 }).re, 0x8000_0000);
    try expectBitsF(complex.subf(.{ .re = 0.0, .im = 0.0 }, .{ .re = 0.0, .im = 0.0 }).im, 0x0000_0000);

    // A NaN contaminates only its own component, keeps its sign and payload,
    // and is quieted.
    const partial = complex.addf(.{ .re = nan_f, .im = 2 }, .{ .re = 3, .im = 4 });
    try expectCanonicalF(partial.re);
    try testing.expectEqual(@as(f32, 6), partial.im);
    // The same invalid rule covers subtraction, where the operands are equal
    // infinities.
    try expectCanonicalF(complex.subf(.{ .re = inf_f, .im = 0 }, .{ .re = inf_f, .im = 0 }).re);
    try expectCanonicalD(complex.subd(.{ .re = 0, .im = neg_inf_d }, .{ .re = 0, .im = neg_inf_d }).im);
    try expectBitsF(complex.addf(.{ .re = neg_nan_f, .im = 0 }, .{ .re = 1, .im = 0 }).re, 0xffc0_0000);
    try expectBitsF(complex.subf(.{ .re = signaling_nan_f, .im = 0 }, .{ .re = 1, .im = 0 }).re, 0x7fc0_0001);
    // With a NaN in each operand of one component, the first one wins.
    try expectBitsF(complex.addf(
        .{ .re = neg_nan_f, .im = 0 },
        .{ .re = signaling_nan_f, .im = 0 },
    ).re, 0xffc0_0000);
    try testing.expectEqual(@as(f32, 4), complex.addf(
        .{ .re = 1, .im = 2 },
        .{ .re = 3, .im = nan_f },
    ).re);

    // Infinities add as IEEE infinities do, and an opposing pair is invalid.
    try expectBitsF(complex.addf(.{ .re = inf_f, .im = 0 }, .{ .re = 1, .im = 0 }).re, 0x7f80_0000);
    try expectCanonicalF(complex.addf(.{ .re = inf_f, .im = 0 }, .{ .re = neg_inf_f, .im = 0 }).re);
    try expectCanonicalD(complex.addd(.{ .re = inf_d, .im = 0 }, .{ .re = neg_inf_d, .im = 0 }).re);
    // An infinity stays in its own component, so the other component adds
    // zero to zero and is a positive zero.
    try expectBitsF(complex.addf(.{ .re = inf_f, .im = 0 }, .{ .re = 0, .im = 0 }).im, 0x0000_0000);
    // Underflow rounds to a signed zero by IEEE addition.
    try expectBitsF(complex.addf(
        .{ .re = 1e-45, .im = 0 },
        .{ .re = -1e-45, .im = 0 },
    ).re, 0x0000_0000);
}

test "mul is the naive four-multiply form" {
    // (3 + 4i)(5 + 6i) = 15 - 24 + (18 + 20)i.
    try testing.expectEqual(Complexf{ .re = -9, .im = 38 }, complex.mulf(
        .{ .re = 3, .im = 4 },
        .{ .re = 5, .im = 6 },
    ));
    // (3 + 4i)(5 - 6i) = 15 + 24 + (-18 + 20)i.
    try testing.expectEqual(Complexd{ .re = 39, .im = 2 }, complex.muld(
        .{ .re = 3, .im = 4 },
        .{ .re = 5, .im = -6 },
    ));
    // A real product is the real part with a zero imaginary part.
    try testing.expectEqual(@as(f32, 6), complex.mulf(.{ .re = 2, .im = 3 }, .{ .re = 3, .im = 0 }).re);
    try testing.expectEqual(@as(f32, 9), complex.mulf(.{ .re = 2, .im = 3 }, .{ .re = 3, .im = 0 }).im);

    // Any NaN operand poisons both components, because every product mixes the
    // two operands of one component with those of the other.
    const poisoned = complex.mulf(.{ .re = 0, .im = nan_f }, .{ .re = 3, .im = 0 });
    try expectCanonicalF(poisoned.re);
    try expectCanonicalF(poisoned.im);
    try expectBitsF(complex.mulf(
        .{ .re = neg_nan_f, .im = 0 },
        .{ .re = 1, .im = 0 },
    ).re, 0xffc0_0000);
    try expectBitsF(complex.mulf(
        .{ .re = 1, .im = 0 },
        .{ .re = 1, .im = signaling_nan_f },
    ).im, 0x7fc0_0001);

    // infinity * 0 is invalid, and the product overflows to infinity.
    try expectCanonicalF(complex.mulf(.{ .re = inf_f, .im = 0 }, .{ .re = 0, .im = 0 }).re);
    try expectBitsF(complex.mulf(.{ .re = 3e38, .im = 0 }, .{ .re = 3, .im = 0 }).re, 0x7f80_0000);
    // Underflow gives a signed zero.
    try expectBitsF(complex.mulf(
        .{ .re = 1e-30, .im = 0 },
        .{ .re = 1e-30, .im = 0 },
    ).re, 0x0000_0000);
}

test "mul by the conjugate of a real value is the squared magnitude" {
    for ([_]f32{ 2, -3, 0.5, 1e-30, 1e30 }) |r| {
        const product = complex.mulf(
            .{ .re = r, .im = 0 },
            complex.conjf(.{ .re = r, .im = 0 }),
        );
        try testing.expectEqual(r * r, product.re);
        // r*(-0) + 0*r is +0 for every nonzero r, and for r = +0 too.
        try expectBitsF(product.im, 0x0000_0000);
    }
    for ([_]f64{ 2, -3, 0.5, 1e-300, 1e300 }) |r| {
        const product = complex.muld(
            .{ .re = r, .im = 0 },
            complex.conjd(.{ .re = r, .im = 0 }),
        );
        try testing.expectEqual(r * r, product.re);
        try expectBitsD(product.im, 0x0000_0000_0000_0000);
    }
}

test "div survives operands the naive denominator overflows" {
    // The naive denominator b.re^2 + b.im^2 is infinite here, so the direct
    // form would answer inf / inf = NaN.
    try expectSquaredSumOverflowsD(1e200, 1e200);
    const diagonal = complex.divd(
        .{ .re = 1e200, .im = 1e200 },
        .{ .re = 1e200, .im = 1e200 },
    );
    try testing.expectApproxEqAbs(@as(f64, 1), diagonal.re, 1e-15);
    try testing.expectApproxEqAbs(@as(f64, 0), diagonal.im, 1e-15);
    try expectFiniteD(diagonal.re);
    try expectFiniteD(diagonal.im);

    // (1 - i) / (1 + i) is exactly -i, so the quotient is a pure imaginary
    // negative one.
    const ratio = complex.divd(
        .{ .re = 1e200, .im = -1e200 },
        .{ .re = 1e200, .im = 1e200 },
    );
    try testing.expectApproxEqAbs(@as(f64, 0), ratio.re, 1e-15);
    try testing.expectApproxEqAbs(@as(f64, -1), ratio.im, 1e-15);

    // The same reduction at a magnitude the naive denominator underflows to
    // zero, where the direct form would answer a nonzero / 0.
    try expectSquareUnderflowsD(1e-200);
    const tiny = complex.divd(
        .{ .re = 1e-200, .im = 0 },
        .{ .re = 1e-200, .im = 0 },
    );
    try testing.expectApproxEqAbs(@as(f64, 1), tiny.re, 1e-15);
    try expectBitsD(tiny.im, 0x0000_0000_0000_0000);

    // A tiny numerator against a large divisor is a genuine underflow.
    const underflow = complex.divd(
        .{ .re = 1e-200, .im = 0 },
        .{ .re = 1e200, .im = 0 },
    );
    try expectBitsD(underflow.re, 0x0000_0000_0000_0000);

    // The binary32 pair, where 1e38 * 1e38 is also infinite.
    try expectSquaredSumOverflowsF(1e38, 1e38);
    const narrow = complex.divf(
        .{ .re = 1e38, .im = 1e38 },
        .{ .re = 1e38, .im = 1e38 },
    );
    try testing.expectApproxEqAbs(@as(f32, 1), narrow.re, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 0), narrow.im, 1e-6);
    try expectFiniteF(narrow.re);
    try expectFiniteF(narrow.im);

    // A moderate divisor still takes the direct path and is exact.
    const exact = complex.divd(
        .{ .re = 3, .im = 4 },
        .{ .re = 1, .im = 2 },
    );
    try testing.expectEqual(@as(f64, 2.2), exact.re);
    try testing.expectEqual(@as(f64, -0.4), exact.im);
}

test "div by zero, by an infinity, and zero over zero" {
    // A finite nonzero numerator over the origin has infinite magnitude and an
    // undefined direction, so the real part is the signed infinity of the
    // naive numerator and the imaginary part is invalid.
    // 1 * +0 + 2 * +0 is +0, so the real part is +infinity.
    const over_zero = complex.divd(.{ .re = 1, .im = 2 }, .{ .re = 0, .im = 0 });
    try expectBitsD(over_zero.re, 0x7ff0_0000_0000_0000);
    try expectCanonicalD(over_zero.im);
    // 1 * -0 + 2 * -0 is -0, so a negative zero divisor reverses the sign.
    const over_negative_zero = complex.divd(.{ .re = 1, .im = 2 }, .{ .re = -0.0, .im = -0.0 });
    try expectBitsD(over_negative_zero.re, 0xfff0_0000_0000_0000);
    try expectCanonicalD(over_negative_zero.im);
    // The imaginary part is the invalid operation infinity minus infinity.
    try expectCanonicalF(complex.divf(.{ .re = -1, .im = 1 }, .{ .re = 0, .im = 0 }).im);

    // Zero over zero is invalid in both components.
    const indeterminate = complex.divd(.{ .re = 0, .im = -0.0 }, .{ .re = 0, .im = 0 });
    try expectCanonicalD(indeterminate.re);
    try expectCanonicalD(indeterminate.im);
    try expectCanonicalF(complex.divf(.{ .re = 0, .im = 0 }, .{ .re = -0.0, .im = -0.0 }).re);

    // A finite numerator over an infinite divisor has zero magnitude. The
    // reciprocal of b is 1 / b.re when b.re is infinite, and -i / b.im when
    // b.im is infinite, so each component is the sign of one numerator
    // component against that one.
    const over_positive_real_infinity = complex.divd(
        .{ .re = -1, .im = 2 },
        .{ .re = inf_d, .im = 0 },
    );
    // -1 / +inf is -0 and 2 / +inf is +0.
    try expectBitsD(over_positive_real_infinity.re, 0x8000_0000_0000_0000);
    try expectBitsD(over_positive_real_infinity.im, 0x0000_0000_0000_0000);
    const over_negative_real_infinity = complex.divd(
        .{ .re = -1, .im = 2 },
        .{ .re = neg_inf_d, .im = 0 },
    );
    // -1 / -inf is +0 and 2 / -inf is -0.
    try expectBitsD(over_negative_real_infinity.re, 0x0000_0000_0000_0000);
    try expectBitsD(over_negative_real_infinity.im, 0x8000_0000_0000_0000);
    // (1 + 2i) / (i inf) is (2 - i) / inf, which is +0 - i0.
    const over_imaginary_infinity = complex.divd(
        .{ .re = 1, .im = 2 },
        .{ .re = 0, .im = inf_d },
    );
    try expectBitsD(over_imaginary_infinity.re, 0x0000_0000_0000_0000);
    try expectBitsD(over_imaginary_infinity.im, 0x8000_0000_0000_0000);

    // An infinite numerator against an infinite divisor, and a divisor with two
    // infinite components, have no direction at all.
    const infinite_over_infinite = complex.divd(
        .{ .re = inf_d, .im = 1 },
        .{ .re = inf_d, .im = 0 },
    );
    try expectCanonicalD(infinite_over_infinite.re);
    try expectCanonicalD(infinite_over_infinite.im);
    try expectCanonicalF(complex.divf(
        .{ .re = 1, .im = 0 },
        .{ .re = inf_f, .im = inf_f },
    ).re);
    // An infinite numerator against a zero divisor is also invalid.
    try expectCanonicalD(complex.divd(.{ .re = inf_d, .im = 0 }, .{ .re = 0, .im = 0 }).re);

    // A NaN operand contaminates both components and keeps its own sign.
    try expectBitsD(complex.divd(
        .{ .re = 0, .im = neg_nan_d },
        .{ .re = 1, .im = 0 },
    ).re, 0xfff8_0000_0000_0000);
    try expectCanonicalF(complex.divf(.{ .re = 1, .im = 0 }, .{ .re = nan_f, .im = 0 }).re);
}

test "div keeps a subnormal quotient that a flushing divide would drop" {
    // The exact quotient of 1 and the maximum finite binary64 rounds to the
    // subnormal 2^-1024 (0x0004000000000000): the exact value 2^-971 /
    // (2^53 - 1) is 2^50 ulps plus one eighth of an ulp, so it rounds down to
    // 2^50 ulps. A target whose binary64 divide flushes the subnormal
    // range answers a signed zero for it, so these rows are the
    // consumer-visible difference.
    // The reduced denominator is one there, so the imaginary part is the
    // positive zero of the exact numerator.
    var real_numerator: Complexd = .{ .re = 1, .im = 0 };
    var real_divisor: Complexd = .{ .re = maximum_d, .im = 0 };
    std.mem.doNotOptimizeAway(&real_numerator);
    std.mem.doNotOptimizeAway(&real_divisor);
    const over_real = complex.divd(real_numerator, real_divisor);
    try expectBitsD(over_real.re, 0x0004_0000_0000_0000);
    try expectBitsD(over_real.im, 0x0000_0000_0000_0000);

    // An imaginary-only divisor turns the same reciprocal a quarter turn
    // first, so the subnormal reaches the real part through the other
    // numerator term.
    var imaginary_numerator: Complexd = .{ .re = 0, .im = 1 };
    var imaginary_divisor: Complexd = .{ .re = 0, .im = maximum_d };
    std.mem.doNotOptimizeAway(&imaginary_numerator);
    std.mem.doNotOptimizeAway(&imaginary_divisor);
    const over_imaginary = complex.divd(imaginary_numerator, imaginary_divisor);
    try expectBitsD(over_imaginary.re, 0x0004_0000_0000_0000);
    try expectBitsD(over_imaginary.im, 0x0000_0000_0000_0000);

    // A divisor with two maximal components leaves a reduced denominator of
    // two, so the quotient halves to the subnormal 2^-1025 in both parts at
    // once and every sign quadrant of it is reachable.
    const quadrants = [_]struct { numerator: Complexd, divisor: Complexd, re: u64, im: u64 }{
        .{
            .numerator = .{ .re = 1, .im = 0 },
            .divisor = .{ .re = maximum_d, .im = maximum_d },
            .re = 0x0002_0000_0000_0000,
            .im = 0x8002_0000_0000_0000,
        },
        .{
            .numerator = .{ .re = -1, .im = 0 },
            .divisor = .{ .re = maximum_d, .im = maximum_d },
            .re = 0x8002_0000_0000_0000,
            .im = 0x0002_0000_0000_0000,
        },
        .{
            .numerator = .{ .re = 1, .im = 0 },
            .divisor = .{ .re = -maximum_d, .im = -maximum_d },
            .re = 0x8002_0000_0000_0000,
            .im = 0x0002_0000_0000_0000,
        },
        .{
            .numerator = .{ .re = -1, .im = 0 },
            .divisor = .{ .re = -maximum_d, .im = -maximum_d },
            .re = 0x0002_0000_0000_0000,
            .im = 0x8002_0000_0000_0000,
        },
    };
    var quadrant_fixtures = quadrants;
    std.mem.doNotOptimizeAway(&quadrant_fixtures);
    for (quadrant_fixtures) |row| {
        const quotient = complex.divd(row.numerator, row.divisor);
        try expectBitsD(quotient.re, row.re);
        try expectBitsD(quotient.im, row.im);
    }

    // One order of magnitude further inside the range, where 1 / 1e308 is the
    // subnormal 0x000730d67819e8d2 rather than a zero.
    var moderate_divisor: Complexd = .{ .re = 1e308, .im = 0 };
    std.mem.doNotOptimizeAway(&moderate_divisor);
    const moderate = complex.divd(real_numerator, moderate_divisor);
    try expectBitsD(moderate.re, 0x0007_30d6_7819_e8d2);
    try expectBitsD(moderate.im, 0x0000_0000_0000_0000);
}

test "abs is overflow-safe and dominates both components" {
    // The samples whose squares are normal numbers, so the textbook magnitude
    // is a valid expectation for them.
    const well_scaled = [_]Complexd{
        .{ .re = 3, .im = 4 },
        .{ .re = -3, .im = 4 },
        .{ .re = 1, .im = 1 },
        .{ .re = -0.0, .im = -0.0 },
        .{ .re = 5, .im = 0 },
        .{ .re = 1e-150, .im = 1e-150 },
    };
    for (well_scaled) |z| {
        const magnitude = complex.absd(z);
        // hypot is never below either component, and round-to-nearest of a
        // value at least as large as a representable one stays at least as large.
        try testing.expect(magnitude >= @abs(z.re));
        try testing.expect(magnitude >= @abs(z.im));
        try testing.expectApproxEqRel(
            @sqrt(z.re * z.re + z.im * z.im),
            magnitude,
            1e-15,
        );
    }
    // The same dominance holds where the textbook magnitude is unusable.
    for ([_]Complexd{
        .{ .re = 1e200, .im = 1e200 },
        .{ .re = 1e-200, .im = 1e-200 },
        .{ .re = 1e300, .im = 1e-300 },
    }) |z| {
        const magnitude = complex.absd(z);
        try testing.expect(magnitude >= @abs(z.re));
        try testing.expect(magnitude >= @abs(z.im));
        try expectFiniteD(magnitude);
    }
    try expectBitsD(complex.absd(.{ .re = 0, .im = -0.0 }), 0x0000_0000_0000_0000);

    // 1e200 * 1e200 is infinite, so the squared form cannot be the reduction.
    try expectSquareOverflowsD(1e200);
    const large = complex.absd(.{ .re = 1e200, .im = 1e200 });
    try testing.expectApproxEqRel(@as(f64, 1.4142135623730951e200), large, 1e-15);
    try expectFiniteD(large);
    const tiny = complex.absd(.{ .re = 1e-200, .im = 1e-200 });
    try testing.expectApproxEqRel(@as(f64, 1.4142135623730951e-200), tiny, 1e-15);
    try testing.expect(tiny != 0);

    // The binary32 pair, where 3e38 is past the reduction's safe window.
    try expectSquaredSumOverflowsF(3e38, 3e38);
    const narrow = complex.absf(.{ .re = 3e38, .im = 3e38 });
    // 3e38 * sqrt(2) is past the largest binary32, so the magnitude overflows.
    try expectBitsF(narrow, 0x7f80_0000);
    const narrow_finite = complex.absf(.{ .re = 1e38, .im = 1e38 });
    try expectFiniteF(narrow_finite);
    try testing.expectApproxEqRel(@as(f32, 1.4142135e38), narrow_finite, 1e-6);

    // An infinite component is an infinite magnitude, and it wins over a NaN
    // one because the value is not merely unknown, it is unbounded.
    try expectBitsD(complex.absd(.{ .re = inf_d, .im = 0 }), 0x7ff0_0000_0000_0000);
    try expectBitsD(complex.absd(.{ .re = 0, .im = neg_inf_d }), 0x7ff0_0000_0000_0000);
    try expectBitsD(complex.absd(.{ .re = inf_d, .im = nan_d }), 0x7ff0_0000_0000_0000);
    // A NaN component has no magnitude at all, so the result is canonical.
    try expectCanonicalD(complex.absd(.{ .re = nan_d, .im = 0 }));
    try expectCanonicalD(complex.absd(.{ .re = 1, .im = signaling_nan_d }));
    try expectCanonicalF(complex.absf(.{ .re = neg_nan_f, .im = 0 }));
}

test "arg covers the quadrants and the signed zeros" {
    try testing.expectApproxEqAbs(pi_f / 4, complex.argf(.{ .re = 1, .im = 1 }), 1e-7);
    try testing.expectApproxEqAbs(3 * pi_f / 4, complex.argf(.{ .re = -1, .im = 1 }), 1e-7);
    try testing.expectApproxEqAbs(-3 * pi_f / 4, complex.argf(.{ .re = -1, .im = -1 }), 1e-7);
    try testing.expectApproxEqAbs(-pi_f / 4, complex.argf(.{ .re = 1, .im = -1 }), 1e-7);
    try testing.expectApproxEqAbs(pi_d / 2, complex.argd(.{ .re = 0, .im = 1 }), 1e-15);
    try testing.expectApproxEqAbs(-pi_d / 2, complex.argd(.{ .re = 0, .im = -1 }), 1e-15);
    try testing.expectApproxEqAbs(pi_d, complex.argd(.{ .re = -1, .im = 1e-300 }), 1e-15);
    try testing.expectApproxEqAbs(pi_d, complex.argd(.{ .re = -1e-300, .im = 0 }), 1e-15);

    // The four signed-zero corners, which the accepted atan2 kernel fixes.
    try expectBitsF(complex.argf(.{ .re = 0.0, .im = 0.0 }), 0x0000_0000);
    try expectBitsF(complex.argf(.{ .re = 0.0, .im = -0.0 }), 0x8000_0000);
    try testing.expectApproxEqAbs(pi_f, complex.argf(.{ .re = -0.0, .im = 0.0 }), 1e-7);
    try testing.expectApproxEqAbs(-pi_f, complex.argf(.{ .re = -0.0, .im = -0.0 }), 1e-7);
    try expectBitsD(complex.argd(.{ .re = 0.0, .im = 0.0 }), 0x0000_0000_0000_0000);
    try expectBitsD(complex.argd(.{ .re = 0.0, .im = -0.0 }), 0x8000_0000_0000_0000);
    try testing.expectApproxEqAbs(pi_d, complex.argd(.{ .re = -0.0, .im = 0.0 }), 1e-15);
    try testing.expectApproxEqAbs(-pi_d, complex.argd(.{ .re = -0.0, .im = -0.0 }), 1e-15);

    // An infinite component has a quadrant but no argument error.
    try testing.expectApproxEqAbs(pi_f / 2, complex.argf(.{ .re = 1, .im = inf_f }), 1e-7);
    try testing.expectApproxEqAbs(-pi_f / 4, complex.argf(.{ .re = inf_f, .im = -inf_f }), 1e-7);
    try testing.expectApproxEqAbs(3 * pi_f / 4, complex.argf(.{ .re = neg_inf_f, .im = inf_f }), 1e-7);
    // A NaN is quieted in place.
    try expectBitsF(complex.argf(.{ .re = neg_nan_f, .im = 0 }), 0xffc0_0000);
    try expectBitsD(complex.argd(.{ .re = 0, .im = signaling_nan_d }), 0x7ff8_0000_0000_0001);
}

test "exp follows the C decomposition" {
    // exp(x + i0) is e^x, so the imaginary part is the input zero.
    try expectBitsF(complex.expf(.{ .re = 0, .im = 0 }).re, 0x3f80_0000);
    try expectBitsF(complex.expf(.{ .re = 0, .im = -0.0 }).im, 0x8000_0000);
    try testing.expectEqual(scalar.expf(1.0), complex.expf(.{ .re = 1, .im = 0 }).re);
    try expectBitsD(complex.expd(.{ .re = 1, .im = -0.0 }).im, 0x8000_0000_0000_0000);

    // |exp(z)| is e^z.re exactly, which ties both components to the real part.
    for ([_]Complexd{ .{ .re = 0.5, .im = 0.5 }, .{ .re = -2, .im = 30 } }) |z| {
        const w = complex.expd(z);
        try testing.expectApproxEqRel(scalar.expd(z.re), complex.absd(w), 1e-15);
    }
    // A purely imaginary argument is a unit magnitude rotation.
    const rotation = complex.expd(.{ .re = 0, .im = 0.7 });
    try testing.expectApproxEqRel(@as(f64, 1), complex.absd(rotation), 1e-15);
    try testing.expectApproxEqAbs(scalar.cosd(0.7), rotation.re, 1e-15);
    try testing.expectApproxEqAbs(scalar.sind(0.7), rotation.im, 1e-15);

    // Real overflow and underflow, and the infinities.
    try expectBitsF(complex.expf(.{ .re = 1000, .im = 0 }).re, 0x7f80_0000);
    try expectBitsF(complex.expf(.{ .re = -1000, .im = 0 }).re, 0x0000_0000);
    try expectBitsD(complex.expd(.{ .re = inf_d, .im = 0 }).re, 0x7ff0_0000_0000_0000);
    try expectBitsD(complex.expd(.{ .re = -inf_d, .im = 0 }).re, 0x0000_0000_0000_0000);
    // +inf with a nonzero argument is an overflow times an undefined rotation.
    const positive_overflow = complex.expd(.{ .re = inf_d, .im = 1 });
    try expectBitsD(positive_overflow.re, 0x7ff0_0000_0000_0000);
    try expectCanonicalD(positive_overflow.im);
    const negative_overflow = complex.expd(.{ .re = neg_inf_d, .im = -1 });
    try expectBitsD(negative_overflow.re, 0x0000_0000_0000_0000);
    try expectCanonicalD(negative_overflow.im);
    // An infinite argument is not on any rotation.
    try expectCanonicalD(complex.expd(.{ .re = 1, .im = inf_d }).re);
    try expectCanonicalF(complex.expf(.{ .re = 1, .im = inf_f }).im);

    // A NaN contaminates both components, because the scale and both rotation
    // factors feed both products.
    try expectCanonicalD(complex.expd(.{ .re = nan_d, .im = 1 }).re);
    try expectCanonicalD(complex.expd(.{ .re = nan_d, .im = 1 }).im);
    try expectBitsD(complex.expd(.{ .re = neg_nan_d, .im = 1 }).re, 0xfff8_0000_0000_0000);
    try expectBitsF(complex.expf(.{ .re = 1, .im = signaling_nan_f }).im, 0x7fc0_0001);
}

test "log is the principal argument decomposition" {
    // log(+1) is a positive zero, and the negative real axis gives +-pi.
    try expectBitsD(complex.logd(.{ .re = 1, .im = 0 }).re, 0x0000_0000_0000_0000);
    try expectBitsD(complex.logd(.{ .re = 1, .im = 0 }).im, 0x0000_0000_0000_0000);
    try testing.expectApproxEqAbs(pi_d, complex.logd(.{ .re = -1, .im = 0 }).im, 1e-15);
    try testing.expectApproxEqAbs(-pi_d, complex.logd(.{ .re = -1, .im = -0.0 }).im, 1e-15);
    try expectBitsD(complex.logd(.{ .re = -1, .im = -0.0 }).re, 0x0000_0000_0000_0000);

    // The origin has an infinite logarithm and a zero angle, not a signed pi.
    const origin = complex.logd(.{ .re = 0, .im = 0 });
    try expectBitsD(origin.re, 0xfff0_0000_0000_0000);
    try expectBitsD(origin.im, 0x0000_0000_0000_0000);
    try expectBitsD(complex.logd(.{ .re = -0.0, .im = -0.0 }).im, 0x8000_0000_0000_0000);
    try expectBitsF(complex.logf(.{ .re = 0, .im = -0.0 }).re, 0xff80_0000);

    // exp(log(z)) is z on the principal branch, and log(z) is ln|z| + i arg(z).
    for ([_]Complexd{ .{ .re = 3, .im = 4 }, .{ .re = -3, .im = 4 }, .{ .re = 1e-200, .im = 1e-200 } }) |z| {
        const w = complex.expd(complex.logd(z));
        try testing.expectApproxEqRel(z.re, w.re, 1e-14);
        try testing.expectApproxEqRel(z.im, w.im, 1e-14);
    }
    const sized = complex.logd(.{ .re = 3, .im = 4 });
    try testing.expectApproxEqRel(scalar.logd(scalar.hypotd(3, 4)), sized.re, 1e-15);
    try testing.expectApproxEqAbs(complex.argd(.{ .re = 3, .im = 4 }), sized.im, 1e-15);

    // log10 divides only the imaginary part by ln(10), so the real part is the
    // base-10 kernel and log10(1000) is exactly three.
    try testing.expectApproxEqAbs(@as(f64, 3), complex.log10d(.{ .re = 1000, .im = 0 }).re, 1e-14);
    const tenth = complex.log10d(.{ .re = 3, .im = 4 });
    try testing.expectApproxEqRel(sized.re / scalar.logd(10), tenth.re, 1e-15);
    try testing.expectApproxEqRel(sized.im / scalar.logd(10), tenth.im, 1e-15);
    try testing.expectApproxEqRel(@as(f64, 2), complex.log10d(.{ .re = 100, .im = 0 }).re, 1e-15);
    // The origin and the infinities keep the log10 form of the log answers.
    try expectBitsD(complex.log10d(.{ .re = 0, .im = 0 }).re, 0xfff0_0000_0000_0000);
    try expectBitsD(complex.log10d(.{ .re = inf_d, .im = 0 }).re, 0x7ff0_0000_0000_0000);
    try testing.expectApproxEqAbs(pi_d / 2 / scalar.logd(10), complex.log10d(.{ .re = 0, .im = inf_d }).im, 1e-14);
    try expectCanonicalD(complex.log10d(.{ .re = inf_d, .im = inf_d }).im);

    // The infinities and the NaN.
    try expectBitsD(complex.logd(.{ .re = inf_d, .im = 1 }).re, 0x7ff0_0000_0000_0000);
    try expectBitsD(complex.logd(.{ .re = inf_d, .im = 1 }).im, 0x0000_0000_0000_0000);
    try expectBitsD(complex.logd(.{ .re = inf_d, .im = -1 }).im, 0x8000_0000_0000_0000);
    try testing.expectApproxEqAbs(pi_d / 2, complex.logd(.{ .re = 1, .im = inf_d }).im, 1e-15);
    try testing.expectApproxEqAbs(-pi_d / 2, complex.logd(.{ .re = -1, .im = neg_inf_d }).im, 1e-15);
    try expectBitsD(complex.logd(.{ .re = inf_d, .im = inf_d }).re, 0x7ff0_0000_0000_0000);
    try expectCanonicalD(complex.logd(.{ .re = inf_d, .im = inf_d }).im);
    try expectCanonicalD(complex.logd(.{ .re = nan_d, .im = 1 }).re);
    try expectBitsD(complex.logd(.{ .re = 1, .im = neg_nan_d }).im, 0xfff8_0000_0000_0000);
    try expectBitsF(complex.log10f(.{ .re = 1, .im = signaling_nan_f }).im, 0x7fc0_0001);
}

test "sqrt is the principal root" {
    // The negative real axis gives a positive imaginary part whose sign the
    // input imaginary zero selects, and a positive zero real part.
    const negative = complex.sqrtf(.{ .re = -4, .im = 0 });
    try expectBitsF(negative.re, 0x0000_0000);
    try testing.expectEqual(@as(f32, 2), negative.im);
    // The imaginary part is -2, not a signed zero, and the real part is +0.
    try testing.expectEqual(@as(f32, -2), complex.sqrtf(.{ .re = -4, .im = -0.0 }).im);
    try expectBitsF(complex.sqrtf(.{ .re = -4, .im = -0.0 }).re, 0x0000_0000);
    try testing.expectApproxEqAbs(@as(f64, 2), complex.sqrtd(.{ .re = -4, .im = 0 }).im, 1e-15);
    // (2 + i)^2 is 3 + 4i exactly, and (1 + 2i)^2 is -3 + 4i exactly.
    const exact = complex.sqrtd(.{ .re = 3, .im = 4 });
    try testing.expectEqual(@as(f64, 2), exact.re);
    try testing.expectEqual(@as(f64, 1), exact.im);
    const exact_negative = complex.sqrtd(.{ .re = -3, .im = 4 });
    try testing.expectEqual(@as(f64, 1), exact_negative.re);
    try testing.expectEqual(@as(f64, 2), exact_negative.im);
    // The root of 3 - 4i is the conjugate of the root of 3 + 4i.
    try testing.expectApproxEqAbs(@as(f64, 2), complex.sqrtd(.{ .re = 3, .im = -4 }).re, 1e-15);
    try testing.expectApproxEqAbs(@as(f64, -1), complex.sqrtd(.{ .re = 3, .im = -4 }).im, 1e-15);

    // The root squared is the input, in both quadrants of the real part.
    for ([_]Complexd{
        .{ .re = 3, .im = 4 },
        .{ .re = -3, .im = 4 },
        .{ .re = -0.5, .im = -0.5 },
        .{ .re = 1e-200, .im = 1e-200 },
    }) |z| {
        const w = complex.muld(complex.sqrtd(z), complex.sqrtd(z));
        try testing.expectApproxEqRel(z.re, w.re, 1e-14);
        try testing.expectApproxEqRel(z.im, w.im, 1e-14);
    }

    // The subnormal and the near-overflow reductions, where halving the
    // magnitude would lose the answer.
    const subnormal = complex.sqrtf(.{ .re = 1e-40, .im = 0 });
    try testing.expectApproxEqRel(@as(f32, 1e-20), subnormal.re, 1e-5);
    try expectBitsF(subnormal.im, 0x0000_0000);
    const huge = complex.sqrtf(.{ .re = 3e38, .im = 0 });
    try expectFiniteF(huge.re);
    try testing.expectApproxEqRel(@as(f32, 1.7320508e19), huge.re, 1e-6);
    // The real part of the root squares to (|z| + a) / 2 and the imaginary
    // part to (|z| - a) / 2, so for z = 1e200 + i 1e200 with |z| = sqrt(2)
    // times 1e200 the two parts are 1.09868e100 and 4.55090e99. They are not
    // equal, because the root of 1 + i is not on the 45 degree line.
    const wide = complex.sqrtd(.{ .re = 1e200, .im = 1e200 });
    try expectFiniteD(wide.re);
    try expectFiniteD(wide.im);
    try testing.expectApproxEqRel(
        @sqrt((@sqrt(2.0) + 1.0) / 2.0) * 1e100,
        wide.re,
        1e-15,
    );
    try testing.expectApproxEqRel(
        @sqrt((@sqrt(2.0) - 1.0) / 2.0) * 1e100,
        wide.im,
        1e-15,
    );
    // Squaring the scaled root twice is the identity: |sqrt(z)|^2 = |z|.
    try testing.expectApproxEqRel(@sqrt(2.0) * 1e200, complex.absd(wide) * complex.absd(wide), 1e-15);

    // The origin keeps the signs of the input, and the infinities grow the
    // parts the input makes unbounded.
    const origin = complex.sqrtf(.{ .re = -0.0, .im = -0.0 });
    try expectBitsF(origin.re, 0x8000_0000);
    try expectBitsF(origin.im, 0x8000_0000);
    const positive_infinity = complex.sqrtd(.{ .re = inf_d, .im = 0 });
    try expectBitsD(positive_infinity.re, 0x7ff0_0000_0000_0000);
    try expectBitsD(positive_infinity.im, 0x0000_0000_0000_0000);
    try expectBitsD(complex.sqrtd(.{ .re = inf_d, .im = 1 }).im, 0x7ff0_0000_0000_0000);
    try expectBitsD(complex.sqrtd(.{ .re = inf_d, .im = -1 }).im, 0xfff0_0000_0000_0000);
    try expectBitsD(complex.sqrtd(.{ .re = 1, .im = inf_d }).im, 0x7ff0_0000_0000_0000);
    try expectBitsD(complex.sqrtd(.{ .re = neg_inf_d, .im = neg_inf_d }).re, 0x7ff0_0000_0000_0000);
    try expectBitsD(complex.sqrtd(.{ .re = neg_inf_d, .im = neg_inf_d }).im, 0x7ff0_0000_0000_0000);
    try expectCanonicalD(complex.sqrtd(.{ .re = nan_d, .im = 1 }).re);
    try expectBitsD(complex.sqrtd(.{ .re = 1, .im = neg_nan_d }).im, 0xfff8_0000_0000_0000);
    try expectBitsF(complex.sqrtf(.{ .re = 1, .im = signaling_nan_f }).im, 0x7fc0_0001);
}

test "sqrt keeps a subnormal component that a flushing divide would drop" {
    // The root of 2 +- i times three units in the last place of the binary64
    // subnormal range divides that imaginary part by exactly 2 * sqrt(2), and
    // 3 / 2.828... is 1.06 units, which rounds to one unit of the subnormal
    // range. A target whose binary64 divide flushes that range answers a
    // signed zero instead. The normal component is sqrt(0.5 * 2 + 0.5 * 2),
    // the root of two, and the input imaginary sign carries it.
    const three_units: f64 = @bitCast(@as(u64, 3));
    const rows = [_]struct { z: Complexd, re: u64, im: u64 }{
        .{
            .z = .{ .re = 2, .im = three_units },
            .re = 0x3ff6_a09e_667f_3bcd,
            .im = 0x0000_0000_0000_0001,
        },
        .{
            .z = .{ .re = 2, .im = -three_units },
            .re = 0x3ff6_a09e_667f_3bcd,
            .im = 0x8000_0000_0000_0001,
        },
        .{
            .z = .{ .re = -2, .im = three_units },
            .re = 0x0000_0000_0000_0001,
            .im = 0x3ff6_a09e_667f_3bcd,
        },
        .{
            .z = .{ .re = -2, .im = -three_units },
            .re = 0x0000_0000_0000_0001,
            .im = 0xbff6_a09e_667f_3bcd,
        },
    };
    var fixtures = rows;
    std.mem.doNotOptimizeAway(&fixtures);
    for (fixtures) |row| {
        const root = complex.sqrtd(row.z);
        try expectBitsD(root.re, row.re);
        try expectBitsD(root.im, row.im);
    }
}

test "pow reproduces the C special cases" {
    // pow(0, 0) is one, and the zero exponent keeps its own sign.
    const origin = complex.powf(.{ .re = 0, .im = 0 }, .{ .re = 0, .im = 0 });
    try testing.expectEqual(@as(f32, 1), origin.re);
    try expectBitsF(origin.im, 0x0000_0000);
    // The sign of the zero that the result takes is the sign of the exponent's
    // imaginary zero.
    const negative_zero_exponent = complex.powf(
        .{ .re = 2, .im = 3 },
        .{ .re = 0.0, .im = -0.0 },
    );
    try testing.expectEqual(@as(f32, 1), negative_zero_exponent.re);
    try expectBitsF(negative_zero_exponent.im, 0x8000_0000);
    try testing.expectEqual(@as(f64, 1), complex.powd(
        .{ .re = inf_d, .im = inf_d },
        .{ .re = 0, .im = 0 },
    ).re);

    // A zero base has no logarithm: a positive real exponent gives zero and
    // every other exponent has infinite magnitude with an undefined direction.
    const zero_to_half = complex.powf(.{ .re = 0, .im = 0 }, .{ .re = 0.5, .im = 0 });
    try expectBitsF(zero_to_half.re, 0x0000_0000);
    try expectBitsF(zero_to_half.im, 0x0000_0000);
    const zero_to_negative = complex.powd(.{ .re = -0.0, .im = 0 }, .{ .re = -1, .im = 0 });
    try expectBitsD(zero_to_negative.re, 0x7ff0_0000_0000_0000);
    try expectCanonicalD(zero_to_negative.im);
    try expectCanonicalF(complex.powf(.{ .re = 0, .im = 0 }, .{ .re = 1, .im = 1 }).im);

    // A negative real base with an exact integer exponent is real, so the
    // imaginary part is an exact zero rather than a reduction error.
    const cube = complex.powf(.{ .re = -8, .im = 0 }, .{ .re = 3, .im = 0 });
    try testing.expectEqual(@as(f32, -512), cube.re);
    try expectBitsF(cube.im, 0x0000_0000);
    const fourth = complex.powd(.{ .re = -8, .im = 0 }, .{ .re = 4, .im = 0 });
    try testing.expectEqual(@as(f64, 4096), fourth.re);
    try expectBitsD(fourth.im, 0x0000_0000_0000_0000);
    // A negative base with a negative integer exponent keeps the parity of the
    // exponent, so (-2)^-3 is -1/8.
    try testing.expectApproxEqRel(
        @as(f64, -0.125),
        complex.powd(.{ .re = -2, .im = 0 }, .{ .re = -3, .im = 0 }).re,
        1e-15,
    );

    // A negative real base with a non-integral exponent takes the principal
    // argument at +pi, so the half power of -4 is 2i with a vanishing real
    // part. That real part is the argument reduction error of the cosine at
    // pi/2, which is 6.12e-17 and far below the tolerance.
    const half = complex.powd(.{ .re = -4, .im = 0 }, .{ .re = 0.5, .im = 0 });
    try testing.expectApproxEqAbs(@as(f64, 0), half.re, 1e-15);
    try testing.expectApproxEqAbs(@as(f64, 2), half.im, 1e-15);
    // The square of the half power is the base, so the branch was taken at
    // arg = +pi and not at arg = -pi.
    const half_squared = complex.muld(half, half);
    try testing.expectApproxEqRel(@as(f64, -4), half_squared.re, 1e-15);
    try testing.expectApproxEqAbs(@as(f64, 0), half_squared.im, 1e-15);
    const narrow_half = complex.powf(.{ .re = -4, .im = 0 }, .{ .re = 0.5, .im = 0 });
    try testing.expectApproxEqAbs(@as(f32, 0), narrow_half.re, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 2), narrow_half.im, 1e-6);
    // exp(y log x) is the definition, so a squaring round trip is exact.
    const squared = complex.powd(.{ .re = 2, .im = 3 }, .{ .re = 2, .im = 0 });
    const product = complex.muld(.{ .re = 2, .im = 3 }, .{ .re = 2, .im = 3 });
    try testing.expectApproxEqRel(product.re, squared.re, 1e-15);
    try testing.expectApproxEqRel(product.im, squared.im, 1e-15);
    // The exponent one is the identity on the principal branch.
    try testing.expectApproxEqRel(@as(f64, 2), complex.powd(.{ .re = 2, .im = 3 }, .{ .re = 1, .im = 0 }).re, 1e-15);
    try testing.expectApproxEqRel(@as(f64, 3), complex.powd(.{ .re = 2, .im = 3 }, .{ .re = 1, .im = 0 }).im, 1e-15);

    // Overflow and underflow come from the base, and the infinities from the
    // same decomposition.
    try expectBitsF(complex.powf(.{ .re = 1e30, .im = 0 }, .{ .re = 2, .im = 0 }).re, 0x7f80_0000);
    try expectBitsF(complex.powf(.{ .re = 1e-30, .im = 0 }, .{ .re = 2, .im = 0 }).re, 0x0000_0000);
    try expectCanonicalD(complex.powd(.{ .re = 1, .im = 0 }, .{ .re = inf_d, .im = 0 }).re);
    // A negative NaN exponent keeps its own sign.
    try expectBitsD(complex.powd(.{ .re = 1, .im = 0 }, .{ .re = neg_nan_d, .im = 0 }).re, 0xfff8_0000_0000_0000);
    try expectBitsF(complex.powf(
        .{ .re = 1, .im = signaling_nan_f },
        .{ .re = 1, .im = 0 },
    ).re, 0x7fc0_0001);
}

test "sin, cos and tan use the exact decompositions" {
    // A real argument leaves the imaginary part at the input zero.
    for ([_]f64{ 0, 1, -1, 0.7, 3 }) |x| {
        try testing.expectEqual(scalar.sind(x), complex.sind(.{ .re = x, .im = 0 }).re);
        try expectBitsD(complex.sind(.{ .re = x, .im = -0.0 }).im, 0x8000_0000_0000_0000);
        try testing.expectEqual(scalar.cosd(x), complex.cosd(.{ .re = x, .im = 0 }).re);
        try testing.expectEqual(scalar.tand(x), complex.tand(.{ .re = x, .im = 0 }).re);
        try expectBitsD(complex.tand(.{ .re = x, .im = -0.0 }).im, 0x8000_0000_0000_0000);
    }
    try expectBitsF(complex.sinf(.{ .re = 0, .im = -0.0 }).re, 0x0000_0000);
    try testing.expectApproxEqAbs(@as(f32, 1), complex.sinf(.{ .re = pi_f / 2, .im = 0 }).re, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 1), complex.cosf(.{ .re = 0, .im = 0 }).re, 1e-7);
    try expectBitsF(complex.tanf(.{ .re = 0, .im = 0 }).re, 0x0000_0000);

    // A pure imaginary argument: sin(ib) = i sinh(b) and cos(ib) = cosh(b).
    try expectBitsD(complex.sind(.{ .re = 0, .im = 1 }).re, 0x0000_0000_0000_0000);
    try testing.expectApproxEqAbs(scalar.sinhd(1), complex.sind(.{ .re = 0, .im = 1 }).im, 1e-15);
    try testing.expectApproxEqAbs(scalar.coshd(1), complex.cosd(.{ .re = 0, .im = 1 }).re, 1e-15);
    try expectBitsD(complex.cosd(.{ .re = 0, .im = 1 }).im, 0x8000_0000_0000_0000);
    // The tangent quotient is sinh(2b) / (cos 2a + cosh 2b) = sinh(2b) /
    // (1 + cosh 2b) = tanh(b), the identity behind tan(ib) = i tanh(b).
    try testing.expectApproxEqRel(
        scalar.tanhd(1),
        complex.tand(.{ .re = 0, .im = 1 }).im,
        1e-14,
    );
    try expectBitsD(complex.tand(.{ .re = 0, .im = 1 }).re, 0x0000_0000_0000_0000);

    // cosh(2b) overflows long before the tangent quotient does, so the binary32
    // pair at b = 100 exercises the overflow reduction.
    try expectOverflowsToInfinityF(scalar.coshf(200));
    const overflowed = complex.tanf(.{ .re = 0, .im = 100 });
    try expectFiniteF(overflowed.re);
    try expectFiniteF(overflowed.im);
    try expectBitsF(overflowed.re, 0x0000_0000);
    try testing.expectApproxEqAbs(@as(f32, 1), overflowed.im, 1e-6);
    // sinh and cosh overflow drives the corresponding parts to infinity.
    try expectBitsF(complex.sinf(.{ .re = 0, .im = 100 }).im, 0x7f80_0000);
    try expectBitsF(complex.cosf(.{ .re = 0, .im = 100 }).re, 0x7f80_0000);
    // sin(0) is an exact zero, so the infinite cosh leaves the real part zero
    // rather than the IEEE NaN that inf * 0 would give.
    try expectBitsF(complex.sinf(.{ .re = 0, .im = -100 }).re, 0x0000_0000);

    // An infinite argument has no defined rotation, and an infinite imaginary
    // argument leaves the real part invalid while the imaginary part grows.
    try expectCanonicalD(complex.sind(.{ .re = inf_d, .im = 1 }).re);
    try expectCanonicalD(complex.cosd(.{ .re = inf_d, .im = 1 }).im);
    const infinite_imaginary = complex.sind(.{ .re = 1, .im = inf_d });
    try expectCanonicalD(infinite_imaginary.re);
    try expectBitsD(infinite_imaginary.im, 0x7ff0_0000_0000_0000);
    const negative_infinite_imaginary = complex.cosd(.{ .re = 1, .im = inf_d });
    try expectBitsD(negative_infinite_imaginary.im, 0xfff0_0000_0000_0000);
    try expectBitsD(complex.tand(.{ .re = 1, .im = neg_inf_d }).im, 0xfff0_0000_0000_0000);

    // A NaN contaminates both parts of every family.
    for ([_]Complexd{ .{ .re = nan_d, .im = 1 }, .{ .re = 1, .im = neg_nan_d } }) |z| {
        try testing.expect(std.math.isNan(complex.sind(z).re));
        try testing.expect(std.math.isNan(complex.sind(z).im));
        try testing.expect(std.math.isNan(complex.cosd(z).re));
        try testing.expect(std.math.isNan(complex.tand(z).im));
    }
    try expectBitsD(complex.sind(.{ .re = 1, .im = signaling_nan_d }).re, 0x7ff8_0000_0000_0001);
}

test "sinh, cosh and tanh use the exact decompositions" {
    for ([_]f64{ 0, 1, -1, 0.7, 20 }) |x| {
        try testing.expectEqual(scalar.sinhd(x), complex.sinhd(.{ .re = x, .im = 0 }).re);
        try testing.expectEqual(scalar.coshd(x), complex.coshd(.{ .re = x, .im = 0 }).re);
        try testing.expectEqual(scalar.tanhd(x), complex.tanhd(.{ .re = x, .im = 0 }).re);
        try expectBitsD(complex.sinhd(.{ .re = x, .im = -0.0 }).im, 0x8000_0000_0000_0000);
    }
    // A pure imaginary argument: sinh(ib) = i sin(b) and cosh(ib) = cos(b).
    try expectBitsD(complex.sinhd(.{ .re = 0, .im = 1 }).re, 0x0000_0000_0000_0000);
    try testing.expectApproxEqAbs(scalar.sind(1), complex.sinhd(.{ .re = 0, .im = 1 }).im, 1e-15);
    try testing.expectApproxEqAbs(scalar.cosd(1), complex.coshd(.{ .re = 0, .im = 1 }).re, 1e-15);
    // The tangent quotient is sin(2b) / cos(2b) = tan(b).
    try testing.expectApproxEqRel(
        scalar.tand(1),
        complex.tanhd(.{ .re = 0, .im = 1 }).im,
        1e-14,
    );
    try expectBitsD(complex.tanhd(.{ .re = 0, .im = 1 }).re, 0x0000_0000_0000_0000);
    // The tangent of a large real argument is a real limit.
    try expectOverflowsToInfinityF(scalar.coshf(100));
    const limit = complex.tanhf(.{ .re = 50, .im = 1 });
    try expectFiniteF(limit.re);
    try testing.expectApproxEqAbs(@as(f32, 1), limit.re, 1e-7);
    try testing.expectApproxEqAbs(@as(f32, 0), limit.im, 1e-7);
    // cosh(50)^2 overflows in binary32, so the limit is exactly one.
    const exactly_one = complex.tanhf(.{ .re = 50, .im = 0 });
    try expectBitsF(exactly_one.re, 0x3f80_0000);
    // sinh(ib) = i sin(b), which stays bounded, while a large real part grows.
    try testing.expectApproxEqAbs(scalar.sinf(100), complex.sinhf(.{ .re = 0, .im = 100 }).im, 1e-5);
    try expectBitsF(complex.sinhf(.{ .re = 100, .im = 1 }).re, 0x7f80_0000);
    try expectBitsF(complex.coshf(.{ .re = 100, .im = 1 }).re, 0x7f80_0000);
    try expectCanonicalF(complex.coshf(.{ .re = 0, .im = inf_f }).re);
    // An infinite real argument has a real tangent limit.
    try expectBitsF(complex.tanhf(.{ .re = inf_f, .im = 1 }).re, 0x3f80_0000);
    try expectBitsF(complex.tanhf(.{ .re = neg_inf_f, .im = 1 }).re, 0xbf80_0000);
    try expectBitsF(complex.tanhf(.{ .re = inf_f, .im = 1 }).im, 0x0000_0000);
    // The pole of tanh at i pi/2, where the denominator cosh(2a) + cos(2b)
    // vanishes. The real numerator is a zero and the imaginary one is not, so
    // the reported value is the canonical quiet NaN and a signed infinity
    // rather than a zero over zero.
    const pole = complex.tanhd(.{ .re = 0, .im = 1.5707963267948966 });
    try expectCanonicalD(pole.re);
    try expectBitsD(pole.im, 0x7ff0_0000_0000_0000);
    try expectBitsD(complex.tanhd(.{ .re = 0, .im = -1.5707963267948966 }).im, 0xfff0_0000_0000_0000);
    // An infinity in either component leaves the decomposition undefined.
    try expectCanonicalD(complex.sinhd(.{ .re = inf_d, .im = 1 }).re);
    try expectCanonicalD(complex.coshd(.{ .re = 1, .im = inf_d }).im);
    try expectCanonicalD(complex.tanhd(.{ .re = 1, .im = inf_d }).re);

    for ([_]Complexd{ .{ .re = nan_d, .im = 1 }, .{ .re = 1, .im = neg_nan_d } }) |z| {
        try testing.expect(std.math.isNan(complex.sinhd(z).im));
        try testing.expect(std.math.isNan(complex.coshd(z).re));
        try testing.expect(std.math.isNan(complex.tanhd(z).re));
    }
    try expectBitsF(complex.coshf(.{ .re = signaling_nan_f, .im = 0 }).re, 0x7fc0_0001);
}

test "asin, acos and atan use the C99 principal branches" {
    // asin(-2) is -pi/2 + i acosh(2), with the real part at the lower end of
    // the principal range and the imaginary part positive.
    const below = complex.asind(.{ .re = -2, .im = 0 });
    try testing.expectApproxEqAbs(-pi_d / 2, below.re, 1e-15);
    try testing.expectApproxEqRel(scalar.acoshd(2), below.im, 1e-15);
    const above = complex.asind(.{ .re = 2, .im = 0 });
    try testing.expectApproxEqAbs(pi_d / 2, above.re, 1e-15);
    try testing.expectApproxEqRel(-scalar.acoshd(2), above.im, 1e-15);
    // A real argument inside the unit interval has no imaginary part.
    try expectBitsD(complex.asind(.{ .re = 0.5, .im = 0 }).im, 0x8000_0000_0000_0000);
    try testing.expectApproxEqAbs(scalar.asind(0.5), complex.asind(.{ .re = 0.5, .im = 0 }).re, 1e-15);
    try expectBitsD(complex.asind(.{ .re = 0, .im = 0 }).re, 0x0000_0000_0000_0000);
    try expectBitsD(complex.asind(.{ .re = 0, .im = 0 }).im, 0x8000_0000_0000_0000);

    // atan(-2) is the real atan, because -2 is inside the principal domain of
    // the argument, whose branch cut is the imaginary axis.
    const atan_negative = complex.atand(.{ .re = -2, .im = 0 });
    try expectBitsD(atan_negative.im, 0x0000_0000_0000_0000);
    try testing.expectApproxEqRel(scalar.atand(2), complex.atand(.{ .re = 2, .im = 0 }).re, 1e-15);
    try expectBitsD(complex.atand(.{ .re = 0, .im = 0 }).re, 0x0000_0000_0000_0000);

    // acos is pi/2 - asin, so acos(0) is pi/2 and the real part stays in
    // [0, pi] on both sides of the unit interval.
    try testing.expectApproxEqAbs(pi_d / 2, complex.acosd(.{ .re = 0, .im = 0 }).re, 1e-15);
    try expectBitsD(complex.acosd(.{ .re = 0, .im = 0 }).im, 0x0000_0000_0000_0000);
    // acos(-2) = pi - i acosh(2), because acos is pi/2 - asin and asin(-2) is
    // -pi/2 + i acosh(2). The real acos of -2 is a NaN and cannot be the
    // reference.
    try testing.expectApproxEqAbs(pi_d, complex.acosd(.{ .re = -2, .im = 0 }).re, 1e-15);
    try testing.expectApproxEqRel(-scalar.acoshd(2), complex.acosd(.{ .re = -2, .im = 0 }).im, 1e-15);
    // acos(-1) is pi exactly.
    try testing.expectApproxEqAbs(pi_d, complex.acosd(.{ .re = -1, .im = 0 }).re, 1e-15);
    try expectBitsD(complex.acosd(.{ .re = -1, .im = 0 }).im, 0x0000_0000_0000_0000);
    // acos(1) is a zero real part, so the real part is pi/2 - pi/2.
    try expectBitsD(complex.acosd(.{ .re = 1, .im = 0 }).re, 0x0000_0000_0000_0000);

    // The inverses that do round trip: sin inverts asin, cos inverts acos, and
    // tan inverts atan. cos(asin(z)) is not z, because the principal range of
    // asin is narrower than that of acos; it is the principal root of 1 - z^2,
    // which is the branch the C rule pins and the stronger check here.
    for ([_]Complexd{ .{ .re = 0.3, .im = 0.4 }, .{ .re = -0.6, .im = 0.8 }, .{ .re = 0, .im = -0.25 } }) |z| {
        try expectComplexApprox(z, complex.sind(complex.asind(z)), 1e-12);
        try expectComplexApprox(z, complex.cosd(complex.acosd(z)), 1e-12);
        try expectComplexApprox(z, complex.tand(complex.atand(z)), 1e-12);
        const square = complex.muld(z, z);
        const one_minus_square: Complexd = .{ .re = 1 - square.re, .im = -square.im };
        try expectComplexApprox(complex.sqrtd(one_minus_square), complex.cosd(complex.asind(z)), 1e-12);
        // A real argument inside the unit interval is the classic round trip.
        try testing.expectApproxEqRel(@as(f64, 0.5), complex.sind(complex.asind(.{ .re = 0.5, .im = 0 })).re, 1e-15);
        try testing.expectApproxEqRel(@as(f64, 0.5), complex.cosd(complex.acosd(.{ .re = 0.5, .im = 0 })).re, 1e-15);
        try testing.expectApproxEqRel(@as(f64, 0.5), complex.tand(complex.atand(.{ .re = 0.5, .im = 0 })).re, 1e-15);
    }

    // The infinities and the NaN.
    for ([_]Complexd{ .{ .re = inf_d, .im = 1 }, .{ .re = 1, .im = neg_inf_d } }) |z| {
        try expectCanonicalD(complex.asind(z).re);
        try expectCanonicalD(complex.acosd(z).im);
        try expectCanonicalD(complex.atand(z).re);
    }
    try expectBitsD(complex.asind(.{ .re = neg_nan_d, .im = 0 }).re, 0xfff8_0000_0000_0000);
    try expectBitsD(complex.atand(.{ .re = signaling_nan_d, .im = 0 }).im, 0x7ff8_0000_0000_0001);
}

test "asinh, acosh and atanh use the C99 principal branches" {
    // acosh(0) is i pi/2, whose real part is an exact zero, and acosh(-1) is
    // i pi because the two roots multiply to zero there.
    const zero = complex.acoshd(.{ .re = 0, .im = 0 });
    try expectBitsD(zero.re, 0x0000_0000_0000_0000);
    try testing.expectApproxEqAbs(pi_d / 2, zero.im, 1e-15);
    const minus_one = complex.acoshd(.{ .re = -1, .im = 0 });
    try expectBitsD(minus_one.re, 0x0000_0000_0000_0000);
    try testing.expectApproxEqAbs(pi_d, minus_one.im, 1e-15);
    // acosh(0.5) is i acos(0.5) = i pi/3, on the cut below +1.
    const below_one = complex.acoshd(.{ .re = 0.5, .im = 0 });
    try expectBitsD(below_one.re, 0x0000_0000_0000_0000);
    try testing.expectApproxEqAbs(pi_d / 3, below_one.im, 1e-15);
    // acosh(1) is the branch point, where the argument is a zero real part.
    try expectBitsD(complex.acoshd(.{ .re = 1, .im = 0 }).re, 0x0000_0000_0000_0000);
    try expectBitsD(complex.acoshd(.{ .re = 1, .im = 0 }).im, 0x0000_0000_0000_0000);
    // Above the cut the result is real and positive.
    try testing.expectApproxEqRel(scalar.acoshd(2), complex.acoshd(.{ .re = 2, .im = 0 }).re, 1e-15);
    try expectBitsD(complex.acoshd(.{ .re = 2, .im = 0 }).im, 0x0000_0000_0000_0000);

    // atanh(2) is log(3) / 2 - i pi/2, where the second logarithm has its
    // argument on the negative real axis.
    const beyond = complex.atanhd(.{ .re = 2, .im = 0 });
    try testing.expectApproxEqRel(scalar.logd(3) / 2, beyond.re, 1e-15);
    try testing.expectApproxEqAbs(-pi_d / 2, beyond.im, 1e-15);
    // atanh(-2) is its conjugate, because -2 is on the same cut.
    const before = complex.atanhd(.{ .re = -2, .im = 0 });
    try testing.expectApproxEqRel(-scalar.logd(3) / 2, before.re, 1e-15);
    try testing.expectApproxEqAbs(pi_d / 2, before.im, 1e-15);
    // Inside the unit interval the result is real.
    try testing.expectApproxEqRel(scalar.atanhd(0.5), complex.atanhd(.{ .re = 0.5, .im = 0 }).re, 1e-15);
    try expectBitsD(complex.atanhd(.{ .re = 0.5, .im = 0 }).im, 0x0000_0000_0000_0000);
    try expectBitsD(complex.atanhd(.{ .re = 0, .im = 0 }).re, 0x0000_0000_0000_0000);
    // The poles of atanh are the branch points where a logarithm diverges.
    try expectBitsD(complex.atanhd(.{ .re = 1, .im = 0 }).re, 0x7ff0_0000_0000_0000);
    try expectBitsD(complex.atanhd(.{ .re = 1, .im = 0 }).im, 0x0000_0000_0000_0000);
    try expectBitsD(complex.atanhd(.{ .re = -1, .im = 0 }).re, 0xfff0_0000_0000_0000);

    // asinh is odd, and the imaginary part is in [-pi/2, pi/2].
    try expectBitsD(complex.asinhd(.{ .re = 0, .im = 0 }).re, 0x0000_0000_0000_0000);
    try testing.expectApproxEqRel(scalar.asinhd(2), complex.asinhd(.{ .re = 2, .im = 0 }).re, 1e-15);
    // asinh is odd, so the negative real part is a real negative number whose
    // imaginary part is the positive zero of the root product.
    try testing.expectApproxEqRel(-scalar.asinhd(2), complex.asinhd(.{ .re = -2, .im = 0 }).re, 1e-15);
    try expectBitsD(complex.asinhd(.{ .re = -2, .im = 0 }).im, 0x0000_0000_0000_0000);
    // asinh(i) = -i asin(i * i) = -i asin(-1) = i pi/2, so the imaginary part
    // is pi/2 and the real part is the exact zero that acosh(1) is.
    const imaginary = complex.asinhd(.{ .re = 0, .im = 1 });
    try testing.expectApproxEqAbs(pi_d / 2, imaginary.im, 1e-15);
    // asinh(i) = -i asin(i * i) = -i asin(-1) = i pi/2, so the real part is
    // the exact zero that acosh(1) is.
    try expectBitsD(imaginary.re, 0x0000_0000_0000_0000);
    // A large real argument is real, because the root dominates the argument.
    const large = complex.asinhd(.{ .re = 1e200, .im = 0 });
    try expectFiniteD(large.re);
    try expectBitsD(large.im, 0x0000_0000_0000_0000);
    try testing.expectApproxEqRel(scalar.logd(2) + scalar.logd(1e200), large.re, 1e-15);

    // The inverses round trip.
    for ([_]Complexd{ .{ .re = 0.3, .im = 0.4 }, .{ .re = -0.6, .im = 0.8 } }) |z| {
        const round = complex.sinhd(complex.asinhd(z));
        try testing.expectApproxEqRel(z.re, round.re, 1e-12);
        try testing.expectApproxEqRel(z.im, round.im, 1e-12);
        const root = complex.coshd(complex.acoshd(z));
        try testing.expectApproxEqRel(z.re, root.re, 1e-12);
        try testing.expectApproxEqRel(z.im, root.im, 1e-12);
        const ratio = complex.tanhd(complex.atanhd(z));
        try testing.expectApproxEqRel(z.re, ratio.re, 1e-12);
        try testing.expectApproxEqRel(z.im, ratio.im, 1e-12);
    }
    // cosh(acosh(2)) is 2, and sinh(asinh(2)) is 2, because those are the
    // defining identities of the inverses.
    const root_two = complex.coshd(complex.acoshd(.{ .re = 2, .im = 0 }));
    try testing.expectApproxEqRel(@as(f64, 2), root_two.re, 1e-14);
    try expectBitsD(root_two.im, 0x0000_0000_0000_0000);
    const sine_two = complex.sinhd(complex.asinhd(.{ .re = 2, .im = 0 }));
    try testing.expectApproxEqRel(@as(f64, 2), sine_two.re, 1e-14);
    try expectBitsD(sine_two.im, 0x0000_0000_0000_0000);

    // The infinities and the NaN.
    for ([_]Complexd{ .{ .re = inf_d, .im = 1 }, .{ .re = 1, .im = neg_inf_d } }) |z| {
        try expectCanonicalD(complex.asinhd(z).re);
        try expectCanonicalD(complex.acoshd(z).im);
        try expectCanonicalD(complex.atanhd(z).re);
    }
    try expectBitsD(complex.acoshd(.{ .re = neg_nan_d, .im = 0 }).re, 0xfff8_0000_0000_0000);
    try expectBitsF(complex.atanhf(.{ .re = 0, .im = signaling_nan_f }).im, 0x7fc0_0001);
}

test "a NaN operand reaches every family as the canonical quiet NaN" {
    // The positive canonical quiet NaN is the input, so the affected component
    // must come back as exactly that bit pattern.
    try expectCanonicalF(complex.addf(.{ .re = nan_f, .im = 0 }, .{ .re = 0, .im = 0 }).re);
    try expectCanonicalF(complex.subf(.{ .re = 0, .im = nan_f }, .{ .re = 0, .im = 0 }).im);
    try expectCanonicalF(complex.mulf(.{ .re = nan_f, .im = 0 }, .{ .re = 0, .im = 0 }).im);
    try expectCanonicalF(complex.divf(.{ .re = 0, .im = 0 }, .{ .re = nan_f, .im = 0 }).re);
    try expectCanonicalF(complex.expf(.{ .re = 0, .im = nan_f }).re);
    try expectCanonicalF(complex.logf(.{ .re = nan_f, .im = 0 }).im);
    try expectCanonicalF(complex.log10f(.{ .re = 0, .im = nan_f }).re);
    try expectCanonicalF(complex.sqrtf(.{ .re = nan_f, .im = 0 }).im);
    try expectCanonicalF(complex.powf(.{ .re = 0, .im = 0 }, .{ .re = nan_f, .im = 0 }).im);
    try expectCanonicalF(complex.sinf(.{ .re = nan_f, .im = 0 }).im);
    try expectCanonicalF(complex.cosf(.{ .re = 0, .im = nan_f }).re);
    try expectCanonicalF(complex.tanf(.{ .re = nan_f, .im = 0 }).re);
    try expectCanonicalF(complex.asinf(.{ .re = 0, .im = nan_f }).im);
    try expectCanonicalF(complex.acosf(.{ .re = nan_f, .im = 0 }).re);
    try expectCanonicalF(complex.atanf(.{ .re = 0, .im = nan_f }).im);
    try expectCanonicalF(complex.sinhf(.{ .re = nan_f, .im = 0 }).im);
    try expectCanonicalF(complex.coshf(.{ .re = 0, .im = nan_f }).re);
    try expectCanonicalF(complex.tanhf(.{ .re = nan_f, .im = 0 }).im);
    try expectCanonicalF(complex.asinhf(.{ .re = 0, .im = nan_f }).re);
    try expectCanonicalF(complex.acoshf(.{ .re = nan_f, .im = 0 }).im);
    try expectCanonicalF(complex.atanhf(.{ .re = 0, .im = nan_f }).re);

    // The binary64 pair, and the two mixed families with two operands.
    try expectCanonicalD(complex.addd(.{ .re = 0, .im = nan_d }, .{ .re = 0, .im = 0 }).im);
    try expectCanonicalD(complex.muld(.{ .re = 0, .im = 0 }, .{ .re = nan_d, .im = 0 }).re);
    try expectCanonicalD(complex.divd(.{ .re = nan_d, .im = 0 }, .{ .re = 1, .im = 0 }).im);
    try expectCanonicalD(complex.powd(.{ .re = 2, .im = 0 }, .{ .re = 0, .im = nan_d }).re);
    try expectCanonicalD(complex.sind(.{ .re = 0, .im = nan_d }).re);
    try expectCanonicalD(complex.tand(.{ .re = 0, .im = nan_d }).im);
    try expectCanonicalD(complex.acosd(.{ .re = 0, .im = nan_d }).re);
    try expectCanonicalD(complex.coshd(.{ .re = 0, .im = nan_d }).im);
    try expectCanonicalD(complex.tanhd(.{ .re = 0, .im = nan_d }).re);
    try expectCanonicalD(complex.asinhd(.{ .re = 0, .im = nan_d }).re);
    try expectCanonicalD(complex.acoshd(.{ .re = 0, .im = nan_d }).im);
    try expectCanonicalD(complex.atanhd(.{ .re = 0, .im = nan_d }).re);

    // A signaling NaN becomes quiet in place, without losing its payload, and
    // no result may leave a signaling NaN behind.
    const quieted_f = complex.expf(.{ .re = signaling_nan_f, .im = 0 });
    try expectBitsF(quieted_f.re, 0x7fc0_0001);
    try expectBitsF(quieted_f.im, 0x7fc0_0001);
    const quieted_d = complex.logd(.{ .re = 0, .im = signaling_nan_d });
    try expectBitsD(quieted_d.re, 0x7ff8_0000_0000_0001);
    try expectBitsD(quieted_d.im, 0x7ff8_0000_0000_0001);
    // A quiet NaN is already the answer, and keeps its sign.
    const quiet_d = complex.sqrtd(.{ .re = neg_nan_d, .im = 0 });
    try expectBitsD(quiet_d.re, 0xfff8_0000_0000_0000);
    try expectBitsD(quiet_d.im, 0xfff8_0000_0000_0000);
}
