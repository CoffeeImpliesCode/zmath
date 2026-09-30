//! Tests for the `zmath.matrix` surface.
//!
//! The expected values here are identities and hand derivations, never pinned
//! bit patterns that were not derived. Two families appear, and each carries the
//! reason it is checked the way it is.
//!
//! EXACT comparisons. The operands are small integers, so every product, every
//! partial sum, and every cofactor term is an integer of magnitude far below
//! 2^24 in binary32 and 2^53 in binary64. Those are exactly representable, the
//! operations are correctly rounded, and the results are therefore exact
//! integers. On those data a bit comparison is the right comparison, and it is
//! what catches a transposed index or a transposed sign.
//!
//! TOLERANT comparisons. Where the operands are decimals that no binary format
//! holds exactly, a rounding difference between two mathematically equal
//! expressions is not a defect, and the tolerance is stated with the magnitude
//! actually observed. Associativity is the clear case: `(A*B)*C` and `A*(B*C)`
//! group the same products differently, so the group that rounds first differs
//! and bit equality is not something a correct implementation can promise.

const std = @import("std");
const testing = std.testing;
const matrix = @import("matrix.zig");
const vector = @import("vector.zig");
const complex = @import("complex.zig");

const Complexf = complex.Complexf;
const Vec2f = vector.Vec2f;
const Vec3f = vector.Vec3f;
const Vec4f = vector.Vec4f;
const Vec4d = vector.Vec4d;
const Mat = matrix.Mat;
const Mat2f = matrix.Mat2f;
const Mat3f = matrix.Mat3f;
const Mat4f = matrix.Mat4f;
const Mat5f = matrix.Mat5f;
const Mat2d = matrix.Mat2d;
const Mat3d = matrix.Mat3d;
const Mat4d = matrix.Mat4d;
const Mat5d = matrix.Mat5d;

/// `1 / n` in binary64, formed as a comptime expression.
///
/// A hand-typed decimal for `1/3` or `1/7` is a guess at a rounding, and a
/// guess that happens to be wrong would make the near-singularity below weaker
/// or stronger than intended. This is the correctly rounded double by
/// construction.
fn recip(comptime n: f64) f64 {
    return 1.0 / n;
}

/// A `Complexf` from its two components.
fn cf(re: f32, im: f32) Complexf {
    return .{ .re = re, .im = im };
}

/// The positive canonical quiet NaN of binary64, from the special-value table
/// of the specification.
const nan_d: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_0000));

/// The positive infinity of binary64, from the same table.
const inf_d: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0000));

// ---------------------------------------------------------------------------
// Storage and access
// ---------------------------------------------------------------------------

test "storage is column-major while the accessor is row-major" {
    // The literal is ROW-major: the first inner array is row 0, the second is
    // row 1. It is deliberately not symmetric, so a transposed implementation
    // cannot pass by accident.
    const m = Mat2f.init(.{ .{ 1, 2 }, .{ 3, 4 } });

    try testing.expectEqual(@as(f32, 1), m.at(0, 0));
    try testing.expectEqual(@as(f32, 2), m.at(0, 1));
    try testing.expectEqual(@as(f32, 3), m.at(1, 0));
    try testing.expectEqual(@as(f32, 4), m.at(1, 1));

    // The same four values, addressed as `data[column][row]`. Each expectation
    // below is the one that breaks if the two orders are confused: 2 lives at
    // row 0 column 1, so it is `at(0, 1)` and `data[1][0]`, never
    // `at(1, 0)` and never `data[0][1]`.
    try testing.expectEqual(@as(f32, 1), m.data[0][0]);
    try testing.expectEqual(@as(f32, 3), m.data[0][1]);
    try testing.expectEqual(@as(f32, 2), m.data[1][0]);
    try testing.expectEqual(@as(f32, 4), m.data[1][1]);

    // The shape is `[columns][rows]`, so the outer index of a 2 by 3 runs over
    // the three columns and the inner index over the two rows.
    const r = Mat(2, 3, f32).init(.{ .{ 1, 2, 3 }, .{ 4, 5, 6 } });
    try testing.expectEqual(@as(usize, 3), r.data.len);
    try testing.expectEqual(@as(usize, 2), r.data[0].len);
    try testing.expectEqual(@as(usize, 2), @TypeOf(r).rows);
    try testing.expectEqual(@as(usize, 3), @TypeOf(r).columns);
    try testing.expectEqual(@as(f32, 4), r.data[0][1]);
    try testing.expectEqual(@as(f32, 2), r.data[1][0]);
    try testing.expectEqual(@as(f32, 6), r.data[2][1]);

    // A column is contiguous in the storage, which is the property that makes
    // the contraction loop worth writing this way round.
    const column_two = r.data[2];
    try testing.expectEqual(@as(f32, 3), column_two[0]);
    try testing.expectEqual(@as(f32, 6), column_two[1]);
}

test "set writes the row the accessor reads" {
    var m = Mat2f.init(.{ .{ 1, 2 }, .{ 3, 4 } });
    m.set(0, 1, 20);
    try testing.expectEqual(@as(f32, 20), m.at(0, 1));
    try testing.expectEqual(@as(f32, 20), m.data[1][0]);
    // The neighbouring row keeps its own value, which a transposed `set` would
    // have overwritten.
    try testing.expectEqual(@as(f32, 3), m.at(1, 0));
    try testing.expectEqual(@as(f32, 3), m.data[0][1]);
    try testing.expectEqual(@as(f32, 1), m.at(0, 0));
    try testing.expectEqual(@as(f32, 4), m.at(1, 1));
}

// ---------------------------------------------------------------------------
// Identity
// ---------------------------------------------------------------------------

test "identity is a two-sided unit for Mat3f and Mat4d" {
    // Entries and components are in [-2, 2], so every product is at most 4 in
    // magnitude and every partial sum at most 16. All of those are exactly
    // representable, so the comparisons below are bit comparisons and a
    // negative-zero or sign slip cannot hide behind a tolerance.
    const m = Mat3f.init(.{
        .{ 1, -2, 0 },
        .{ 2, 1, 1 },
        .{ 0, 1, -2 },
    });
    const id = Mat3f.identity();

    const right = matrix.mul(3, id, m);
    const left = matrix.mul(3, m, id);
    for (0..3) |row| {
        for (0..3) |column| {
            try testing.expectEqual(m.at(row, column), right.at(row, column));
            try testing.expectEqual(m.at(row, column), left.at(row, column));
        }
    }

    const v = Vec3f.init(1, -2, 3);
    const applied = m.mulVec(v);
    const round_trip = m.mulVec(id.mulVec(v));
    try testing.expectEqual(applied.x, round_trip.x);
    try testing.expectEqual(applied.y, round_trip.y);
    try testing.expectEqual(applied.z, round_trip.z);
    // The identity leaves a vector alone, and it does so bit for bit: each
    // lane is `1 * x` followed by two `0 * y` terms, and a non-zero `x` plus a
    // signed zero is that same non-zero `x`.
    const through_identity = id.mulVec(v);
    try testing.expectEqual(@as(f32, 1), through_identity.x);
    try testing.expectEqual(@as(f32, -2), through_identity.y);
    try testing.expectEqual(@as(f32, 3), through_identity.z);
    // The hand computation of `m * v` for the literal above:
    // row 0: 1*1 + (-2)*(-2) + 0*3 = 1 + 4 + 0 = 5
    // row 1: 2*1 + 1*(-2) + 1*3 = 2 - 2 + 3 = 3
    // row 2: 0*1 + 1*(-2) + (-2)*3 = 0 - 2 - 6 = -8
    try testing.expectEqual(@as(f32, 5), applied.x);
    try testing.expectEqual(@as(f32, 3), applied.y);
    try testing.expectEqual(@as(f32, -8), applied.z);

    // The same identity in binary64 at size four, on a matrix with no zero
    // entry, so the leading `0 * m` terms of a contraction cannot turn a
    // negative zero into a positive one on the way through.
    const m4 = Mat4d.init(.{
        .{ 1, -2, 0, 1 },
        .{ 2, 1, 1, 0 },
        .{ 0, 1, -2, 2 },
        .{ 1, 0, 2, -1 },
    });
    const id4 = Mat4d.identity();
    const right4 = matrix.mul(4, id4, m4);
    const left4 = matrix.mul(4, m4, id4);
    for (0..4) |row| {
        for (0..4) |column| {
            try testing.expectEqual(m4.at(row, column), right4.at(row, column));
            try testing.expectEqual(m4.at(row, column), left4.at(row, column));
        }
    }
    const v4 = Vec4d.init(1, -2, 3, -1);
    const applied4 = m4.mulVec(v4);
    const through4 = m4.mulVec(id4.mulVec(v4));
    try testing.expect(applied4.eq(through4));
    try testing.expect(id4.mulVec(v4).eq(v4));
}

test "identity is available for a non-square shape and means the rectangular one" {
    // Documented behaviour, not an accident: `at(row, column)` is one where the
    // indices agree. The rectangular identity is a unit for every shape that
    // chains with it, so a square-only guard would protect nothing.
    const r = Mat(2, 3, f32).identity();
    try testing.expectEqual(@as(f32, 1), r.at(0, 0));
    try testing.expectEqual(@as(f32, 1), r.at(1, 1));
    try testing.expectEqual(@as(f32, 0), r.at(0, 1));
    try testing.expectEqual(@as(f32, 0), r.at(0, 2));
    try testing.expectEqual(@as(f32, 0), r.at(1, 2));
}

// ---------------------------------------------------------------------------
// Transpose and rectangular products
// ---------------------------------------------------------------------------

test "transpose is an involution and swaps the indices" {
    const m = Mat(2, 3, f32).init(.{ .{ 1, 2, 3 }, .{ 4, 5, 6 } });
    const t = m.transpose();
    try testing.expectEqual(@as(usize, 3), @TypeOf(t).rows);
    try testing.expectEqual(@as(usize, 2), @TypeOf(t).columns);
    for (0..2) |row| {
        for (0..3) |column| {
            try testing.expectEqual(m.at(row, column), t.at(column, row));
        }
    }
    // Derived from the literal above: column 0 of the transpose is the first
    // row of the original, and so on.
    try testing.expectEqual(@as(f32, 1), t.at(0, 0));
    try testing.expectEqual(@as(f32, 4), t.at(0, 1));
    try testing.expectEqual(@as(f32, 2), t.at(1, 0));
    try testing.expectEqual(@as(f32, 5), t.at(1, 1));
    try testing.expectEqual(@as(f32, 3), t.at(2, 0));
    try testing.expectEqual(@as(f32, 6), t.at(2, 1));

    const back = t.transpose();
    for (0..2) |row| {
        for (0..3) |column| {
            try testing.expectEqual(m.at(row, column), back.at(row, column));
        }
    }
}

test "a rectangular product has the contracted shape" {
    // The row-major literal is
    //     [ 1 2 3 ]
    //     [ 4 5 6 ]
    // and its transpose is
    //     [ 1 4 ]
    //     [ 2 5 ]
    //     [ 3 6 ]
    // so the product is the two dot products of each original row with each
    // original column:
    //     [ 1*1+2*2+3*3  1*4+2*5+3*6 ]   [ 14  32 ]
    //     [ 4*1+5*2+6*3  4*4+5*5+6*6 ] = [ 32  77 ]
    // Every product and sum is an integer below 100, so these are exact.
    const m = Mat(2, 3, f32).init(.{ .{ 1, 2, 3 }, .{ 4, 5, 6 } });
    const product = matrix.mul(2, m, m.transpose());
    try testing.expectEqual(@as(usize, 2), @TypeOf(product).rows);
    try testing.expectEqual(@as(usize, 2), @TypeOf(product).columns);
    try testing.expectEqual(@as(f32, 14), product.at(0, 0));
    try testing.expectEqual(@as(f32, 32), product.at(0, 1));
    try testing.expectEqual(@as(f32, 32), product.at(1, 0));
    try testing.expectEqual(@as(f32, 77), product.at(1, 1));
}

// ---------------------------------------------------------------------------
// Entrywise arithmetic
// ---------------------------------------------------------------------------

test "add, sub, neg, and scale against hand-computed matrices" {
    const a = Mat2f.init(.{ .{ 1, 2 }, .{ 3, 4 } });
    const b = Mat2f.init(.{ .{ 10, 20 }, .{ 30, 40 } });

    const sum = a.add(b);
    try testing.expectEqual(@as(f32, 11), sum.at(0, 0));
    try testing.expectEqual(@as(f32, 22), sum.at(0, 1));
    try testing.expectEqual(@as(f32, 33), sum.at(1, 0));
    try testing.expectEqual(@as(f32, 44), sum.at(1, 1));

    const difference = b.sub(a);
    try testing.expectEqual(@as(f32, 9), difference.at(0, 0));
    try testing.expectEqual(@as(f32, 18), difference.at(0, 1));
    try testing.expectEqual(@as(f32, 27), difference.at(1, 0));
    try testing.expectEqual(@as(f32, 36), difference.at(1, 1));

    const negated = a.neg();
    try testing.expectEqual(@as(f32, -1), negated.at(0, 0));
    try testing.expectEqual(@as(f32, -2), negated.at(0, 1));
    try testing.expectEqual(@as(f32, -3), negated.at(1, 0));
    try testing.expectEqual(@as(f32, -4), negated.at(1, 1));

    const doubled = a.scale(2);
    try testing.expectEqual(@as(f32, 2), doubled.at(0, 0));
    try testing.expectEqual(@as(f32, 4), doubled.at(0, 1));
    try testing.expectEqual(@as(f32, 6), doubled.at(1, 0));
    try testing.expectEqual(@as(f32, 8), doubled.at(1, 1));

    // `neg` maps `-0` to `+0` and `+0` to `-0`, so it is an involution on a
    // zero: one application destroys the sign and two applications give it
    // back. That is worth pinning because `x == -x` holds for both zeros, so
    // an equality test alone could not see either behaviour.
    var with_negative_zero = a;
    with_negative_zero.set(0, 0, -0.0);
    const once = with_negative_zero.neg();
    try testing.expectEqual(@as(f32, 0.0), once.at(0, 0));
    try testing.expect(std.math.signbit(once.at(0, 0)) == false);
    const twice = once.neg();
    try testing.expectEqual(@as(f32, 0.0), twice.at(0, 0));
    try testing.expect(std.math.signbit(twice.at(0, 0)) == true);

    // Additive and scalar identities on exact data.
    const self_sum = a.add(a);
    for (0..2) |row| {
        for (0..2) |column| {
            try testing.expectEqual(doubled.at(row, column), self_sum.at(row, column));
        }
    }
    const zero_difference = a.sub(a);
    for (0..2) |row| {
        for (0..2) |column| {
            try testing.expectEqual(@as(f32, 0.0), zero_difference.at(row, column));
        }
    }
}

test "entry arithmetic follows the special-value table" {
    // `spec/core-math-2026-08-22.md:769-801` has two rules that a host unit
    // does not satisfy on its own, and both are observable here.
    //
    // One: a NaN operand propagates with its sign and diagnostic payload and
    // only its quiet bit added. A SIGNALING NaN, 0x7f800001, therefore comes
    // back as the quiet 0x7fc00001 rather than as a payload the host chose.
    const signaling: f32 = @bitCast(@as(u32, 0x7f80_0001));
    var with_signaling = Mat2f.fill(0);
    with_signaling.set(0, 1, signaling);
    const propagated = with_signaling.add(Mat2f.fill(1));
    try testing.expectEqual(@as(u32, 0x7fc0_0001), @as(u32, @bitCast(propagated.at(0, 1))));
    // The neighbouring entry is untouched, so the propagation is per entry and
    // not a whole-matrix effect.
    try testing.expectEqual(@as(f32, 1), propagated.at(0, 0));

    // A quiet NaN is returned unchanged, and the sign of a NaN is preserved.
    const negative_canonical: f32 = @bitCast(@as(u32, 0xffc0_0000));
    var with_negative = Mat2f.fill(0);
    with_negative.set(1, 0, negative_canonical);
    const carried = with_negative.add(Mat2f.fill(1));
    try testing.expectEqual(@as(u32, 0xffc0_0000), @as(u32, @bitCast(carried.at(1, 0))));

    // Two: an invalid operation on NON-NaN operands gives the POSITIVE
    // canonical quiet NaN, 0x7fc00000. A host SSE unit produces the NEGATIVE
    // indefinite, 0xffc00000, for `+infinity + -infinity`, so without the
    // substitution this module would emit a different payload on a different
    // target. That substitution is the whole point of the check.
    const positive_infinity: f32 = @bitCast(@as(u32, 0x7f80_0000));
    const negative_infinity: f32 = @bitCast(@as(u32, 0xff80_0000));
    var left = Mat2f.fill(0);
    var right = Mat2f.fill(0);
    left.set(0, 0, positive_infinity);
    right.set(0, 0, negative_infinity);
    const invalid = left.add(right);
    try testing.expect(std.math.isNan(invalid.at(0, 0)));
    try testing.expectEqual(@as(u32, 0x7fc0_0000), @as(u32, @bitCast(invalid.at(0, 0))));
    // The rest of the matrix is still arithmetic, not a blanket NaN.
    try testing.expectEqual(@as(f32, 0), invalid.at(0, 1));
    try testing.expectEqual(@as(f32, 0), invalid.at(1, 1));
    // And a matrix carrying a NaN is not finite, which is what an inverse
    // rejects.
    try testing.expect(!propagated.isFinite());
    try testing.expect(Mat2f.fill(1).isFinite());
    // An infinity is not finite either, so the two are not conflated.
    try testing.expect(!left.isFinite());
}

test "fill puts one value in every slot" {
    const m = Mat3f.fill(7);
    for (0..3) |row| {
        for (0..3) |column| {
            try testing.expectEqual(@as(f32, 7), m.at(row, column));
        }
    }
}

// ---------------------------------------------------------------------------
// Products
// ---------------------------------------------------------------------------

test "mul is associative within a tight tolerance" {
    // `(A*B)*C` and `A*(B*C)` are the same matrix over the rationals, but the
    // two groupings round their partial sums in a different order, so bit
    // equality is not something a correct implementation can promise. The
    // tolerance is well above the difference actually observed and well below
    // the size of the result: replaying these three matrices in binary64 with
    // the same evaluation order gives a largest entry difference of 1.78e-15
    // against a largest entry of 29.199.
    const a = Mat3d.init(.{
        .{ 1.1, 0.3, 1.3 },
        .{ 3.7, 0.1, 0.2 },
        .{ 1.9, 2.3, 0.2 },
    });
    const b = Mat3d.init(.{
        .{ 1.1, 2.9, 0.1 },
        .{ 2.7, 2.3, 0.4 },
        .{ 0.1, 0.2, 1.3 },
    });
    const c = Mat3d.init(.{
        .{ 1.3, 0.2, 0.4 },
        .{ 0.2, 2.3, 1.3 },
        .{ 0.1, 1.9, 2.9 },
    });

    const left_associated = matrix.mul(3, matrix.mul(3, a, b), c);
    const right_associated = matrix.mul(3, a, matrix.mul(3, b, c));
    for (0..3) |row| {
        for (0..3) |column| {
            try testing.expectApproxEqAbs(
                right_associated.at(row, column),
                left_associated.at(row, column),
                1.0e-12,
            );
        }
    }
}

test "mul distributes over add within a tight tolerance" {
    // Same caveat and same reason as the associativity test: the two sides sum
    // the same products in a different grouping. The observed largest
    // difference on these operands is 1.78e-15 against a largest entry of
    // 16.89.
    const a = Mat3d.init(.{
        .{ 1.1, 0.3, 1.3 },
        .{ 3.7, 0.1, 0.2 },
        .{ 1.9, 2.3, 0.2 },
    });
    const b = Mat3d.init(.{
        .{ 1.1, 2.9, 0.1 },
        .{ 2.7, 2.3, 0.4 },
        .{ 0.1, 0.2, 1.3 },
    });
    const c = Mat3d.init(.{
        .{ 1.3, 0.2, 0.4 },
        .{ 0.2, 2.3, 1.3 },
        .{ 0.1, 1.9, 2.9 },
    });

    const distributed = matrix.mul(3, a, b.add(c));
    const expanded = matrix.mul(3, a, b).add(matrix.mul(3, a, c));
    for (0..3) |row| {
        for (0..3) |column| {
            try testing.expectApproxEqAbs(
                expanded.at(row, column),
                distributed.at(row, column),
                1.0e-12,
            );
        }
    }
}

test "composition is right-to-left" {
    // This is the test the contract's order clause is about.
    //
    //     A = [ 1 2 0 ]     B = [ 0 1 0 ]     v = [ 1 2 3 ]
    //         [ 0 1 3 ]         [ 1 0 2 ]
    //         [ 2 0 1 ]         [ 0 0 1 ]
    //
    // Applying B first: B v = [0*1+1*2+0*3, 1*1+0*2+2*3, 0*1+0*2+1*3]
    //                   = [2, 7, 3],
    // and then A gives A (B v) = [1*2+2*7+0*3, 0*2+1*7+3*3, 2*2+0*7+1*3]
    //                            = [16, 16, 7].
    //
    // Applying A first: A v = [1+4+0, 0+2+9, 2+0+3] = [5, 11, 5],
    // and then B gives B (A v) = [0*5+1*11+0*5, 1*5+0*11+2*5, 0*5+0*11+1*5]
    //                            = [11, 15, 5].
    //
    // The two are 5 apart in the first component, so a product that reversed
    // the composition could not come near either. Every product and sum above
    // is an exact small integer, so the tolerance below is slack that the
    // values do not need; the assertion is stronger than it looks.
    const a = Mat3f.init(.{
        .{ 1, 2, 0 },
        .{ 0, 1, 3 },
        .{ 2, 0, 1 },
    });
    const b = Mat3f.init(.{
        .{ 0, 1, 0 },
        .{ 1, 0, 2 },
        .{ 0, 0, 1 },
    });
    const v = Vec3f.init(1, 2, 3);
    const tolerance: f32 = 1.0e-5;

    const b_then_a = a.mulVec(b.mulVec(v));
    const composed = matrix.mul(3, a, b).mulVec(v);
    try testing.expectApproxEqAbs(@as(f32, 16), b_then_a.x, tolerance);
    try testing.expectApproxEqAbs(@as(f32, 16), b_then_a.y, tolerance);
    try testing.expectApproxEqAbs(@as(f32, 7), b_then_a.z, tolerance);
    try testing.expectApproxEqAbs(@as(f32, 16), composed.x, tolerance);
    try testing.expectApproxEqAbs(@as(f32, 16), composed.y, tolerance);
    try testing.expectApproxEqAbs(@as(f32, 7), composed.z, tolerance);

    const a_then_b = b.mulVec(a.mulVec(v));
    try testing.expectApproxEqAbs(@as(f32, 11), a_then_b.x, tolerance);
    try testing.expectApproxEqAbs(@as(f32, 15), a_then_b.y, tolerance);
    try testing.expectApproxEqAbs(@as(f32, 5), a_then_b.z, tolerance);
    // The other order is a different matrix and gives a different answer. This
    // is the pair that fails if `mul` is reversed, so it is asserted
    // separately rather than left to the three lines above.
    try testing.expectApproxEqAbs(@as(f32, 5), b_then_a.x - a_then_b.x, tolerance);

    // The same derivation at size two, where the two products differ in every
    // entry: A [1 2; 3 4] times B [0 1; 1 0] is [2 1; 4 3], and the reverse
    // product B*A is [3 4; 1 2].
    const a2 = Mat2f.init(.{ .{ 1, 2 }, .{ 3, 4 } });
    const b2 = Mat2f.init(.{ .{ 0, 1 }, .{ 1, 0 } });
    const v2 = Vec2f.init(1, 0);
    const ab = matrix.mul(2, a2, b2);
    const ba = matrix.mul(2, b2, a2);
    try testing.expectEqual(@as(f32, 2), ab.at(0, 0));
    try testing.expectEqual(@as(f32, 1), ab.at(0, 1));
    try testing.expectEqual(@as(f32, 4), ab.at(1, 0));
    try testing.expectEqual(@as(f32, 3), ab.at(1, 1));
    try testing.expectEqual(@as(f32, 3), ba.at(0, 0));
    try testing.expectEqual(@as(f32, 4), ba.at(0, 1));
    try testing.expectEqual(@as(f32, 1), ba.at(1, 0));
    try testing.expectEqual(@as(f32, 2), ba.at(1, 1));
    const ab_v = ab.mulVec(v2);
    const a_after_b = a2.mulVec(b2.mulVec(v2));
    try testing.expectEqual(@as(f32, 2), ab_v.x);
    try testing.expectEqual(@as(f32, 4), ab_v.y);
    try testing.expectEqual(@as(f32, 2), a_after_b.x);
    try testing.expectEqual(@as(f32, 4), a_after_b.y);
}

test "mulVec matches a hand-computed 3x3 product" {
    //     [ 1 2 3 ]   [ 1 ]   [ 1*1+2*1+3*1 ]   [  6 ]
    //     [ 0 1 4 ] * [ 1 ] = [ 0*1+1*1+4*1 ] = [  5 ]
    //     [ 5 6 0 ]   [ 1 ]   [ 5*1+6*1+0*1 ]   [ 11 ]
    const m = Mat3f.init(.{
        .{ 1, 2, 3 },
        .{ 0, 1, 4 },
        .{ 5, 6, 0 },
    });
    const result = m.mulVec(Vec3f.init(1, 1, 1));
    try testing.expectEqual(@as(f32, 6), result.x);
    try testing.expectEqual(@as(f32, 5), result.y);
    try testing.expectEqual(@as(f32, 11), result.z);

    // The zero at row 1 column 0 must contribute nothing, and the four at
    // row 1 column 2 must contribute its full value, so moving the input
    // changes the result in exactly the way the derivation says.
    const weighted = m.mulVec(Vec3f.init(2, 1, 1));
    try testing.expectEqual(@as(f32, 2 + 2 + 3), weighted.x);
    try testing.expectEqual(@as(f32, 0 + 1 + 4), weighted.y);
    try testing.expectEqual(@as(f32, 10 + 6 + 0), weighted.z);
}

test "mulVec uses a plain array where no named vector type matches" {
    // Five rows, so `vector_type` is `[5]f32` and the lane plumbing takes the
    // array branch. Every row is a constant, so the result is four times that
    // constant.
    const m = Mat(5, 4, f32).init(.{
        .{ 1, 1, 1, 1 },
        .{ 2, 2, 2, 2 },
        .{ 3, 3, 3, 3 },
        .{ 4, 4, 4, 4 },
        .{ 5, 5, 5, 5 },
    });
    const result = m.mulVec([5]f32{ 1, 1, 1, 1, 0 });
    for (0..5) |row| {
        try testing.expectEqual(@as(f32, @floatFromInt(row + 1)) * 4, result[row]);
    }
}

// ---------------------------------------------------------------------------
// Determinant
// ---------------------------------------------------------------------------

test "determinantMat3f against hand-computed values" {
    //     [ 1 2 3 ]
    //     [ 4 5 6 ]   = 1*(5*10-6*8) - 2*(4*10-6*7) + 3*(4*8-5*7)
    //     [ 7 8 10]                = 1*2      - 2*(-2)     + 3*(-3)
    //                               = 2 + 4 - 9 = -3
    const m = Mat3f.init(.{
        .{ 1, 2, 3 },
        .{ 4, 5, 6 },
        .{ 7, 8, 10 },
    });
    try testing.expectEqual(@as(f32, -3), matrix.determinantMat3f(m));
    // Every product and partial sum is an integer below 100, so -3 is the
    // exact value and not a rounded one.
    try testing.expectEqual(@as(f32, -3), matrix.determinantMat3d(Mat3d.init(.{
        .{ 1, 2, 3 },
        .{ 4, 5, 6 },
        .{ 7, 8, 10 },
    })));

    // A transpose leaves the determinant alone, which is the classic identity
    // and one that a transposed storage layout would break.
    try testing.expectEqual(@as(f32, -3), matrix.determinantMat3f(m.transpose()));

    // Swapping two rows negates the determinant.
    try testing.expectEqual(@as(f32, 3), matrix.determinantMat3f(Mat3f.init(.{
        .{ 4, 5, 6 },
        .{ 1, 2, 3 },
        .{ 7, 8, 10 },
    })));

    // A positive value, to pin the sign of the second cofactor as well:
    // 0*(0*6-3*5) - 1*(1*6-3*4) + 2*(1*5-0*4) = 0 + 6 + 10 = 16.
    try testing.expectEqual(@as(f32, 16), matrix.determinantMat3f(Mat3f.init(.{
        .{ 0, 1, 2 },
        .{ 1, 0, 3 },
        .{ 4, 5, 6 },
    })));

    try testing.expectEqual(@as(f32, 1), matrix.determinantMat3f(Mat3f.identity()));

    // Two equal rows, so the determinant is mathematically zero:
    // 1*(5*3-6*2) - 2*(4*3-6*1) + 3*(4*2-5*1) = 3 - 12 + 9 = 0.
    // The three parenthesized differences are exact small integers, so the sum
    // is exactly zero rather than nearly zero, and that is what is asserted.
    const singular = Mat3f.init(.{
        .{ 1, 2, 3 },
        .{ 4, 5, 6 },
        .{ 1, 2, 3 },
    });
    try testing.expectEqual(@as(f32, 0), matrix.determinantMat3f(singular));
    // Any two equal rows, by the same argument.
    try testing.expectEqual(@as(f32, 0), matrix.determinantMat3f(Mat3f.init(.{
        .{ 0, 0, 0 },
        .{ 7, 1, 3 },
        .{ 7, 1, 3 },
    })));
}

test "the closed forms and the cofactor recursion agree with hand values" {
    // 1*4 - 2*3 = -2.
    try testing.expectEqual(@as(f32, -2), matrix.determinantMat2f(Mat2f.init(.{
        .{ 1, 2 },
        .{ 3, 4 },
    })));
    try testing.expectEqual(@as(f64, -2), matrix.determinantMat2d(Mat2d.init(.{
        .{ 1, 2 },
        .{ 3, 4 },
    })));

    // A triangular matrix has the product of its diagonal entries, and the
    // cofactor expansion must find that without being told the matrix is
    // triangular. All the off-diagonal entries here are non-zero on purpose, so
    // an implementation that only works on a diagonal cannot pass.
    const t4 = Mat4d.init(.{
        .{ 1, 0, 0, 0 },
        .{ 7, 2, 0, 0 },
        .{ 3, 5, 3, 0 },
        .{ 1, 1, 1, 4 },
    });
    try testing.expectEqual(@as(f64, 24), matrix.determinantMat4d(t4));
    try testing.expectEqual(@as(f32, 24), matrix.determinantMat4f(Mat4f.init(.{
        .{ 1, 0, 0, 0 },
        .{ 7, 2, 0, 0 },
        .{ 3, 5, 3, 0 },
        .{ 1, 1, 1, 4 },
    })));

    const t5 = Mat5d.init(.{
        .{ 1, 0, 0, 0, 0 },
        .{ 7, 2, 0, 0, 0 },
        .{ 3, 5, 3, 0, 0 },
        .{ 1, 1, 1, 4, 0 },
        .{ 2, 0, 1, 1, 5 },
    });
    try testing.expectEqual(@as(f64, 120), matrix.determinantMat5d(t5));
    try testing.expectEqual(@as(f32, 120), matrix.determinantMat5f(Mat5f.init(.{
        .{ 1, 0, 0, 0, 0 },
        .{ 7, 2, 0, 0, 0 },
        .{ 3, 5, 3, 0, 0 },
        .{ 1, 1, 1, 4, 0 },
        .{ 2, 0, 1, 1, 5 },
    })));

    // An identity at every size the contract names, so a broken recursion
    // cannot hide behind a size that is never exercised.
    try testing.expectEqual(@as(f32, 1), matrix.determinantMat2f(Mat2f.identity()));
    try testing.expectEqual(@as(f64, 1), matrix.determinantMat2d(Mat2d.identity()));
    try testing.expectEqual(@as(f32, 1), matrix.determinantMat4f(Mat4f.identity()));
    try testing.expectEqual(@as(f64, 1), matrix.determinantMat4d(Mat4d.identity()));
    try testing.expectEqual(@as(f32, 1), matrix.determinantMat5f(Mat5f.identity()));
    try testing.expectEqual(@as(f64, 1), matrix.determinantMat5d(Mat5d.identity()));

    // 1*2*3 = 6 by the same argument, and the two triangular shapes above with
    // their last diagonal entry raised give 2, 6, and 8 respectively, which is
    // what makes the value sensitive to that entry.
    const triangle = Mat3f.init(.{
        .{ 1, 0, 0 },
        .{ 9, 2, 0 },
        .{ 4, 3, 3 },
    });
    try testing.expectEqual(@as(f32, 6), matrix.determinantMat3f(triangle));
}

test "the determinant of a singular Mat4d is far below its entry scale" {
    // The first and last rows are equal, so the exact determinant is zero. The
    // entries are reciprocals of small primes, so the arithmetic does NOT
    // cancel exactly: replaying the cofactor expansion in binary64 in the
    // documented order gives -4.743e-20.
    const singular = Mat4d.init(.{
        .{ 1, recip(3), recip(7), recip(11) },
        .{ recip(13), recip(17), recip(19), recip(23) },
        .{ recip(29), recip(31), recip(37), recip(41) },
        .{ 1, recip(3), recip(7), recip(11) },
    });
    const value = matrix.determinantMat4d(singular);
    try testing.expect(singular.isFinite());
    // The largest entry is 1, so a determinant of order 1 is the natural scale
    // for a well-conditioned matrix of this size, and 1e-12 is twelve orders
    // below it.
    try testing.expect(@abs(value) < 1.0e-12);

    // The same matrix with the top-left entry of the last row raised by 1e-6 is
    // NOT singular, and its determinant is -2.179e-11. Without this the test
    // above would also pass for a determinant that is simply small, and it
    // would say nothing about rank.
    const perturbed = Mat4d.init(.{
        .{ 1, recip(3), recip(7), recip(11) },
        .{ recip(13), recip(17), recip(19), recip(23) },
        .{ recip(29), recip(31), recip(37), recip(41) },
        .{ 1 + 1.0e-6, recip(3), recip(7), recip(11) },
    });
    const perturbed_value = matrix.determinantMat4d(perturbed);
    try testing.expect(@abs(perturbed_value) > 1.0e-12);
    // Six orders of magnitude separate them, and the observed separation is
    // about 4.6e8, so the margin is not a guess.
    try testing.expect(@abs(perturbed_value) > 1.0e6 * @abs(value));
}

// ---------------------------------------------------------------------------
// Inverse
// ---------------------------------------------------------------------------

test "inverseMat4d inverts a well-conditioned matrix" {
    // A with determinant 99, and an inverse that is exactly
    //     [  4/33  -2/33   1/33  -4/11 ]
    //     [ -10/99  5/99  14/99  10/33 ]
    //     [ 29/99 -64/99  -1/99   4/33 ]
    //     [  5/99  47/99  -7/99  -5/33 ]
    // over the rationals, which the derivation below re-establishes: multiplying
    // that matrix by A gives the identity in exact arithmetic.
    //
    // The determinant is an exact 99 by the same argument as the triangular
    // cases: every cofactor term is a product of integers of magnitude below
    // a hundred.
    const a = Mat4d.init(.{
        .{ 1, 2, 3, 4 },
        .{ 0, 1, 0, 2 },
        .{ 5, 6, 0, 0 },
        .{ -2, 1, 1, 1 },
    });
    try testing.expectEqual(@as(f64, 99), matrix.determinantMat4d(a));

    const expected = [4][4]f64{
        .{ 4.0 / 33.0, -2.0 / 33.0, 1.0 / 33.0, -4.0 / 11.0 },
        .{ -10.0 / 99.0, 5.0 / 99.0, 14.0 / 99.0, 10.0 / 33.0 },
        .{ 29.0 / 99.0, -64.0 / 99.0, -1.0 / 99.0, 4.0 / 33.0 },
        .{ 5.0 / 99.0, 47.0 / 99.0, -7.0 / 99.0, -5.0 / 33.0 },
    };
    const inverse = matrix.inverseMat4d(a) orelse return error.TestUnexpectedResult;
    for (0..4) |row| {
        for (0..4) |column| {
            try testing.expectApproxEqAbs(
                expected[row][column],
                inverse.at(row, column),
                1.0e-12,
            );
        }
    }

    // The identity that matters, and the one that holds for ANY invertible
    // matrix rather than only for this one. Replaying the Gauss-Jordan
    // reduction of this matrix in binary64 gives a worst residual of 2.22e-16,
    // so 1e-12 has five orders of margin.
    const product = matrix.mul(4, inverse, a);
    const id = Mat4d.identity();
    for (0..4) |row| {
        for (0..4) |column| {
            try testing.expectApproxEqAbs(
                id.at(row, column),
                product.at(row, column),
                1.0e-12,
            );
        }
    }
    const other_order = matrix.mul(4, a, inverse);
    for (0..4) |row| {
        for (0..4) |column| {
            try testing.expectApproxEqAbs(
                id.at(row, column),
                other_order.at(row, column),
                1.0e-12,
            );
        }
    }
}

test "inverseMat4d rejects a singular matrix" {
    // Rows 0 and 3 are equal, so the last pivot collapses to exactly zero.
    const repeated_row = Mat4d.init(.{
        .{ 1, 2, 3, 4 },
        .{ 5, 6, 7, 8 },
        .{ 9, 10, 11, 12 },
        .{ 1, 2, 3, 4 },
    });
    try testing.expectEqual(@as(f64, 0), matrix.determinantMat4d(repeated_row));
    try testing.expect(matrix.inverseMat4d(repeated_row) == null);

    // A zero row, which is singular and is detected at the first pivot.
    const zero_row = Mat4d.init(.{
        .{ 1, 0, 0, 0 },
        .{ 0, 0, 0, 0 },
        .{ 0, 0, 1, 0 },
        .{ 0, 0, 0, 1 },
    });
    try testing.expect(matrix.inverseMat4d(zero_row) == null);

    // The zero matrix. The largest entry is zero, so the pivot test compares
    // against a zero scale and every candidate fails it.
    const zero = Mat4d.fill(0);
    try testing.expect(matrix.inverseMat4d(zero) == null);

    // A singular matrix whose cofactor determinant is NOT exactly zero, which
    // is the case the pivot test exists for. The rows are linearly dependent
    // over the rationals but the expansion leaves a rounding residue.
    const nearly = Mat4d.init(.{
        .{ 1, recip(3), recip(7), recip(11) },
        .{ recip(13), recip(17), recip(19), recip(23) },
        .{ recip(29), recip(31), recip(37), recip(41) },
        .{ 1, recip(3), recip(7), recip(11) },
    });
    try testing.expect(!nearly.isFinite() == false);
    try testing.expect(@abs(matrix.determinantMat4d(nearly)) > 0);
    try testing.expect(matrix.inverseMat4d(nearly) == null);
}

test "the pivot threshold rejects a small pivot and accepts a merely uneven one" {
    // The threshold is `2^-50` for binary64, the unit roundoff times the
    // `2^(4-1)` element growth of the reduction. This matrix has a last pivot
    // of 1e-14, which is about eleven binary orders above the threshold and
    // therefore still inverts, with an inverse entry of 1e14. The arithmetic is
    // exact here: 1e-14 is a dyadic rational in binary64 and every entry of
    // this matrix is a power of two times a power of two, so the reduction
    // divides and multiplies without losing anything until the last step.
    var uneven = Mat4d.identity();
    uneven.set(3, 3, 1.0e-14);
    const inverse = matrix.inverseMat4d(uneven) orelse return error.TestUnexpectedResult;
    try testing.expectApproxEqRel(@as(f64, 1.0e14), inverse.at(3, 3), 1.0e-9);
    try testing.expectApproxEqAbs(@as(f64, 1), inverse.at(0, 0), 1.0e-12);

    // One binary order below the threshold now, and the matrix is reported as
    // singular. The value is a power of two, so the comparison in the
    // implementation is the only thing that can decide it.
    var too_small = Mat4d.identity();
    too_small.set(3, 3, 0x1p-51);
    try testing.expect(matrix.inverseMat4d(too_small) == null);

    // And just above it, so the boundary is pinned from both sides.
    var just_big_enough = Mat4d.identity();
    just_big_enough.set(3, 3, 0x1p-49);
    try testing.expect(matrix.inverseMat4d(just_big_enough) != null);

    // The threshold is a documented public constant, not a hidden epsilon.
    try testing.expectEqual(@as(f64, 0x1p-50), matrix.inverse_pivot_threshold_d);
    try testing.expectEqual(@as(f32, 0x1p-21), matrix.inverse_pivot_threshold_f);
}

test "inverseMat4d rejects a non-finite entry" {
    const with_nan = Mat4d.init(.{
        .{ 1, 0, 0, 0 },
        .{ 0, 1, 0, 0 },
        .{ 0, 0, 1, 0 },
        .{ 0, 0, 0, nan_d },
    });
    try testing.expect(!with_nan.isFinite());
    try testing.expect(matrix.inverseMat4d(with_nan) == null);

    const with_infinity = Mat4d.init(.{
        .{ 1, 0, 0, 0 },
        .{ 0, 1, 0, 0 },
        .{ 0, 0, 1, 0 },
        .{ 0, 0, 0, inf_d },
    });
    try testing.expect(!with_infinity.isFinite());
    try testing.expect(matrix.inverseMat4d(with_infinity) == null);

    // A non-finite entry off the diagonal, where a reduction that checked only
    // the pivot column would not notice it. The contract's clause is "or
    // non-finite", so the rejection is on the whole entry set.
    const off_diagonal = Mat4d.init(.{
        .{ 1, 0, 0, 0 },
        .{ 0, 1, 0, 0 },
        .{ 0, 0, 1, inf_d },
        .{ 0, 0, 0, 1 },
    });
    try testing.expect(matrix.inverseMat4d(off_diagonal) == null);

    // What this test does NOT establish: that the check runs BEFORE the
    // arithmetic. A NaN propagates through the reduction, so a check on the
    // result alone would return null as well. The ordering is a property of
    // the implementation, not something an observation from outside can
    // distinguish.
}

test "inverseMat2f matches a hand-written inverse" {
    // [ 4 7 ] has determinant 4*6 - 7*2 = 24 - 14 = 10, and the 2 by 2 adjugate
    // formula gives its inverse as
    //     [  6/10  -7/10 ]
    //     [ -2/10   4/10 ]
    //     [  0.6   -0.7  ]
    //     [ -0.2    0.4  ]
    // None of those decimals is representable in binary32, so the comparison
    // is by tolerance; the Gauss-Jordan reduction reproduces them to one unit
    // in the last place.
    const m = Mat2f.init(.{
        .{ 4, 7 },
        .{ 2, 6 },
    });
    try testing.expectEqual(@as(f32, 10), matrix.determinantMat2f(m));

    const inverse = matrix.inverseMat2f(m) orelse return error.TestUnexpectedResult;
    try testing.expectApproxEqAbs(@as(f32, 0.6), inverse.at(0, 0), 1.0e-6);
    try testing.expectApproxEqAbs(@as(f32, -0.7), inverse.at(0, 1), 1.0e-6);
    try testing.expectApproxEqAbs(@as(f32, -0.2), inverse.at(1, 0), 1.0e-6);
    try testing.expectApproxEqAbs(@as(f32, 0.4), inverse.at(1, 1), 1.0e-6);

    const product = matrix.mul(2, inverse, m);
    const id = Mat2f.identity();
    for (0..2) |row| {
        for (0..2) |column| {
            try testing.expectApproxEqAbs(id.at(row, column), product.at(row, column), 1.0e-6);
        }
    }

    // The binary64 sibling, where the same fractions are much closer to
    // exact and the tolerance can be far tighter.
    const md = Mat2d.init(.{
        .{ 4, 7 },
        .{ 2, 6 },
    });
    const inverse_d = matrix.inverseMat2d(md) orelse return error.TestUnexpectedResult;
    try testing.expectApproxEqAbs(@as(f64, 0.6), inverse_d.at(0, 0), 1.0e-15);
    try testing.expectApproxEqAbs(@as(f64, -0.7), inverse_d.at(0, 1), 1.0e-15);
    try testing.expectApproxEqAbs(@as(f64, -0.2), inverse_d.at(1, 0), 1.0e-15);
    try testing.expectApproxEqAbs(@as(f64, 0.4), inverse_d.at(1, 1), 1.0e-15);

    // A singular 2 by 2 is still rejected, through the same pivot test.
    try testing.expect(matrix.inverseMat2f(Mat2f.init(.{
        .{ 1, 2 },
        .{ 2, 4 },
    })) == null);
}

test "the size-four and size-three inverses invert their own matrices" {
    // The same identity as for size four, on a matrix whose determinant is a
    // small integer, so the cofactor expansion is exact and the reduction has
    // nothing to lose: 1*(1*1 - 3*0) - 2*(0*1 - 3*2) + 0*(0*0 - 1*2)
    // = 1 + 12 + 0 = 13.
    const a = Mat3d.init(.{
        .{ 1, 2, 0 },
        .{ 0, 1, 3 },
        .{ 2, 0, 1 },
    });
    try testing.expectEqual(@as(f64, 13), matrix.determinantMat3d(a));
    const inverse = matrix.inverseMat3d(a) orelse return error.TestUnexpectedResult;
    const product = matrix.mul(3, inverse, a);
    const id = Mat3d.identity();
    for (0..3) |row| {
        for (0..3) |column| {
            try testing.expectApproxEqAbs(id.at(row, column), product.at(row, column), 1.0e-14);
        }
    }

    const a4 = Mat4f.init(.{
        .{ 2, 0, 0, 1 },
        .{ 1, 3, 0, 2 },
        .{ 0, 1, 4, 3 },
        .{ 1, 1, 1, 1 },
    });
    const inverse_f = matrix.inverseMat4f(a4) orelse return error.TestUnexpectedResult;
    const product_f = matrix.mul(4, inverse_f, a4);
    const id_f = Mat4f.identity();
    for (0..4) |row| {
        for (0..4) |column| {
            try testing.expectApproxEqAbs(id_f.at(row, column), product_f.at(row, column), 1.0e-5);
        }
    }
    try testing.expect(matrix.inverseMat3f(Mat3f.init(.{
        .{ 1, 2, 3 },
        .{ 2, 4, 6 },
        .{ 1, 1, 1 },
    })) == null);
}

// ---------------------------------------------------------------------------
// Homogeneous transforms
// ---------------------------------------------------------------------------

test "mulMat4fVec4f applies a homogeneous transform to a point" {
    // Vectors are COLUMNS, so the translation of a column-vector convention is
    // the last COLUMN of the matrix, and in a row-major literal that is the
    // last entry of each of the first three rows. Writing the translation into
    // the last row instead would compute `tx*x + ty*y + tz*z + w`, which is a
    // dot product and not a translation.
    const translation = Mat4f.init(.{
        .{ 1, 0, 0, 10 },
        .{ 0, 1, 0, 20 },
        .{ 0, 0, 1, 30 },
        .{ 0, 0, 0, 1 },
    });
    const moved = matrix.mulMat4fVec4f(translation, Vec4f.init(1, 2, 3, 1));
    try testing.expectEqual(@as(f32, 11), moved.x);
    try testing.expectEqual(@as(f32, 22), moved.y);
    try testing.expectEqual(@as(f32, 33), moved.z);
    try testing.expectEqual(@as(f32, 1), moved.w);

    // A direction has w = 0 and must come back untranslated.
    const direction = matrix.mulMat4fVec4f(translation, Vec4f.init(1, 2, 3, 0));
    try testing.expectEqual(@as(f32, 1), direction.x);
    try testing.expectEqual(@as(f32, 2), direction.y);
    try testing.expectEqual(@as(f32, 3), direction.z);
    try testing.expectEqual(@as(f32, 0), direction.w);

    // Scale by two, three, four: 1->2, 2->6, 3->12, and the homogeneous
    // coordinate is untouched.
    const scale = Mat4f.init(.{
        .{ 2, 0, 0, 0 },
        .{ 0, 3, 0, 0 },
        .{ 0, 0, 4, 0 },
        .{ 0, 0, 0, 1 },
    });
    const scaled = matrix.mulMat4fVec4f(scale, Vec4f.init(1, 2, 3, 1));
    try testing.expectEqual(@as(f32, 2), scaled.x);
    try testing.expectEqual(@as(f32, 6), scaled.y);
    try testing.expectEqual(@as(f32, 12), scaled.z);
    try testing.expectEqual(@as(f32, 1), scaled.w);

    // A right-handed quarter turn about z sends (1, 0, 0) to (0, 1, 0).
    const rotation = Mat4f.init(.{
        .{ 0, -1, 0, 0 },
        .{ 1, 0, 0, 0 },
        .{ 0, 0, 1, 0 },
        .{ 0, 0, 0, 1 },
    });
    const turned = matrix.mulMat4fVec4f(rotation, Vec4f.init(1, 0, 0, 1));
    try testing.expectEqual(@as(f32, 0), turned.x);
    try testing.expectEqual(@as(f32, 1), turned.y);
    try testing.expectEqual(@as(f32, 0), turned.z);
    try testing.expectEqual(@as(f32, 1), turned.w);

    // Translation after rotation, which is where the composition order shows
    // up: rotate (1, 0, 0) to (0, 1, 0) and then add (10, 20, 30) gives
    // (10, 21, 30). The other order would add first and then rotate, giving
    // (-20, 10, 30).
    const composed = matrix.mul(4, translation, rotation);
    const result = matrix.mulMat4fVec4f(composed, Vec4f.init(1, 0, 0, 1));
    try testing.expectEqual(@as(f32, 10), result.x);
    try testing.expectEqual(@as(f32, 21), result.y);
    try testing.expectEqual(@as(f32, 30), result.z);
    try testing.expectEqual(@as(f32, 1), result.w);

    // The binary64 sibling, on a translation by tenths, where the products are
    // not exact but the residual is far below the tolerance.
    const translation_d = Mat4d.init(.{
        .{ 1, 0, 0, 0.1 },
        .{ 0, 1, 0, 0.2 },
        .{ 0, 0, 1, 0.3 },
        .{ 0, 0, 0, 1 },
    });
    const moved_d = matrix.mulMat4dVec4d(translation_d, Vec4d.init(1, 2, 3, 1));
    try testing.expectApproxEqAbs(@as(f64, 1.1), moved_d.x, 1.0e-15);
    try testing.expectApproxEqAbs(@as(f64, 2.2), moved_d.y, 1.0e-15);
    try testing.expectApproxEqAbs(@as(f64, 3.3), moved_d.z, 1.0e-15);
    try testing.expectEqual(@as(f64, 1), moved_d.w);
}

// ---------------------------------------------------------------------------
// Complex entries
// ---------------------------------------------------------------------------

test "a complex instantiation flows through the same kernels" {
    //     [ 1+2i   3  ]
    //     [  0    1-i ]
    // Every operation below is exact in binary32: the naive four-multiply
    // complex product sums and differences of products of integers of magnitude
    // below ten.
    const a = Mat(2, 2, Complexf).init(.{
        .{ cf(1, 2), cf(3, 0) },
        .{ cf(0, 0), cf(1, -1) },
    });
    const b = Mat(2, 2, Complexf).init(.{
        .{ cf(0, 0), cf(1, 0) },
        .{ cf(1, 0), cf(0, 0) },
    });

    // The entrywise sum.
    const sum = a.add(b);
    try testing.expectEqual(@as(f32, 1), sum.at(0, 0).re);
    try testing.expectEqual(@as(f32, 2), sum.at(0, 0).im);
    try testing.expectEqual(@as(f32, 4), sum.at(0, 1).re);
    try testing.expectEqual(@as(f32, 0), sum.at(0, 1).im);
    try testing.expectEqual(@as(f32, 1), sum.at(1, 0).re);
    try testing.expectEqual(@as(f32, 0), sum.at(1, 0).im);
    try testing.expectEqual(@as(f32, 1), sum.at(1, 1).re);
    try testing.expectEqual(@as(f32, -1), sum.at(1, 1).im);

    // The product, contracted along the inner index. Column by column:
    // (0,0) = (1+2i)*0 + 3*1 = 3
    // (0,1) = (1+2i)*1 + 3*0 = 1+2i
    // (1,0) = 0*0 + (1-i)*1 = 1-i
    // (1,1) = 0*1 + (1-i)*0 = 0
    const product = matrix.mul(2, a, b);
    try testing.expectEqual(@as(f32, 3), product.at(0, 0).re);
    try testing.expectEqual(@as(f32, 0), product.at(0, 0).im);
    try testing.expectEqual(@as(f32, 1), product.at(0, 1).re);
    try testing.expectEqual(@as(f32, 2), product.at(0, 1).im);
    try testing.expectEqual(@as(f32, 1), product.at(1, 0).re);
    try testing.expectEqual(@as(f32, -1), product.at(1, 0).im);
    try testing.expectEqual(@as(f32, 0), product.at(1, 1).re);
    try testing.expectEqual(@as(f32, 0), product.at(1, 1).im);

    // The determinant keeps the entry type: (1+2i)*(1-i) - 3*0
    // = (1*1 - 2*(-1)) + i(1*(-1) + 2*1) = 3 + i.
    const determinant_value = a.determinant();
    try testing.expectEqual(@as(f32, 3), determinant_value.re);
    try testing.expectEqual(@as(f32, 1), determinant_value.im);

    // Negation of a complex entry flips BOTH components, and the pin is the
    // real part: `1 + 2i` negates to `-1 - 2i`, which is not the conjugate
    // `1 - 2i`. The last pair is the same distinction from the other side,
    // `1 - i` negating to `-1 + i`.
    const negated = a.neg();
    try testing.expectEqual(@as(f32, -1), negated.at(0, 0).re);
    try testing.expectEqual(@as(f32, -2), negated.at(0, 0).im);
    try testing.expectEqual(@as(f32, -3), negated.at(0, 1).re);
    try testing.expectEqual(@as(f32, 0), negated.at(0, 1).im);
    try testing.expectEqual(@as(f32, -1), negated.at(1, 1).re);
    try testing.expectEqual(@as(f32, 1), negated.at(1, 1).im);
    // Negation is an involution on a complex entry as well, which a
    // conjugation-based implementation would also pass, so it is asserted for
    // the composition and not as a substitute for the four pins above.
    const round_trip = negated.neg();
    try testing.expectEqual(@as(f32, 1), round_trip.at(0, 0).re);
    try testing.expectEqual(@as(f32, 2), round_trip.at(0, 0).im);

    // A complex matrix takes a plain array, because no `Vec*` type is complex
    // and inventing one would contradict the one-generic rule.
    const applied = a.mulVec([2]Complexf{ cf(0, 0), cf(1, 0) });
    try testing.expectEqual(@as(f32, 3), applied[0].re);
    try testing.expectEqual(@as(f32, 0), applied[0].im);
    try testing.expectEqual(@as(f32, 1), applied[1].re);
    try testing.expectEqual(@as(f32, -1), applied[1].im);

    // A complex 3 by 3 reaches the cofactor recursion, whose inner call lands
    // on a complex 2 by 2 and therefore on `det2`. The determinant is
    // (1+2i)*(1-i) - 3*0 = 3 + i, and the 3 by 3 below is upper triangular in
    // its first two rows so the recursion sees exactly that.
    const c3 = Mat(3, 3, Complexf).init(.{
        .{ cf(1, 2), cf(0, 0), cf(0, 0) },
        .{ cf(3, 0), cf(1, -1), cf(0, 0) },
        .{ cf(0, 0), cf(0, 0), cf(2, 0) },
    });
    const det3_value = c3.determinant();
    // The matrix is block upper triangular, so the determinant is the product
    // of the two diagonal blocks: (3 + i) * 2 = 6 + 2i.
    try testing.expectEqual(@as(f32, 6), det3_value.re);
    try testing.expectEqual(@as(f32, 2), det3_value.im);
    try testing.expectEqual(Complexf, @TypeOf(det3_value));
    try testing.expectEqual(f32, @TypeOf(a).element_component);
}

// ---------------------------------------------------------------------------
// The entry-type restriction and the named sizes
// ---------------------------------------------------------------------------

test "the entry type guard accepts exactly four types" {
    // A `@compileError` cannot be observed from a passing test: naming an
    // unsupported entry type would fail to compile the whole test binary, and
    // that failure is the behaviour under test. What IS observable is the
    // predicate the guard consults, and the guard consults nothing else, so
    // this pins the accepted set from both sides: a type it accepts is
    // instantiable and binds to the right component format, and a type it
    // rejects is not in the set. Widening or narrowing `entry_types` without
    // touching the guard is what this test would miss, and the guard reads
    // `entry_types` directly, so there is nothing to widen.
    try testing.expectEqual(@as(usize, 4), matrix.entry_types.len);
    try testing.expect(matrix.acceptsEntryType(f32));
    try testing.expect(matrix.acceptsEntryType(f64));
    try testing.expect(matrix.acceptsEntryType(Complexf));
    try testing.expect(matrix.acceptsEntryType(complex.Complexd));

    // Types outside the set. `f16` is in the zmath scalar surface but is not an
    // entry type, which is the check worth making: a caller must not be able to
    // smuggle it in by assuming that any float works.
    try testing.expect(!matrix.acceptsEntryType(u8));
    try testing.expect(!matrix.acceptsEntryType(i32));
    try testing.expect(!matrix.acceptsEntryType(f16));
    try testing.expect(!matrix.acceptsEntryType(bool));
    try testing.expect(!matrix.acceptsEntryType(void));
    try testing.expect(!matrix.acceptsEntryType([2]f32));
    try testing.expect(!matrix.acceptsEntryType(vector.Vec3f));
    try testing.expect(!matrix.acceptsEntryType(struct {}));

    // The positive branch of the guard: each accepted type instantiates and
    // binds its component format.
    try testing.expectEqual(f32, Mat(2, 2, f32).element_component);
    try testing.expectEqual(f64, Mat(2, 2, f64).element_component);
    try testing.expectEqual(f32, Mat(2, 2, Complexf).element_component);
    try testing.expectEqual(f64, Mat(2, 2, complex.Complexd).element_component);
    try testing.expectEqual(f32, Mat(2, 2, f32).Element);
    try testing.expectEqual(f64, Mat(2, 2, f64).Element);
    try testing.expectEqual(Complexf, Mat(2, 2, Complexf).Element);
    try testing.expectEqual(complex.Complexd, Mat(2, 2, complex.Complexd).Element);
}

comptime {
    // The eight named sizes, in the shapes `API.md:924-931` gives them. A
    // compile-time block rather than a test, because a wrong shape is a
    // compile error everywhere else and never a runtime observation.
    const expected = .{
        .{ Mat2f, 2, 2, f32 },
        .{ Mat3f, 3, 3, f32 },
        .{ Mat4f, 4, 4, f32 },
        .{ Mat5f, 5, 5, f32 },
        .{ Mat2d, 2, 2, f64 },
        .{ Mat3d, 3, 3, f64 },
        .{ Mat4d, 4, 4, f64 },
        .{ Mat5d, 5, 5, f64 },
    };
    for (expected) |entry| {
        const M = entry[0];
        if (M.rows != entry[1] or M.columns != entry[2]) {
            @compileError("a named matrix size has the wrong shape");
        }
        if (M.Element != entry[3]) {
            @compileError("a named matrix size has the wrong entry type");
        }
    }

    // `Mat4f` and `Mat4d` are the homogeneous 3D sizes and `Mat5f` and
    // `Mat5d` the homogeneous 4D ones. The arithmetic cannot see the
    // difference, so the only thing a caller can rely on is the shape, and the
    // reason those two sizes exist is recorded on the aliases themselves.
    if (Mat4f.vector_type != Vec4f or Mat4d.vector_type != Vec4d) {
        @compileError("a 4 by 4 matrix must apply itself to a four-component vector");
    }
    if (Mat3f.vector_type != Vec3f) {
        @compileError("a 3 by 3 f32 matrix must apply itself to a Vec3f");
    }
    if (Mat2f.vector_type != Vec2f) {
        @compileError("a 2 by 2 f32 matrix must apply itself to a Vec2f");
    }
    if (Mat5f.vector_type != [5]f32) {
        @compileError("a 5 by 5 f32 matrix must apply itself to a five-element array");
    }
    // The row count alone decides the vector type; the column count has nothing
    // to do with it, so a two-row matrix of any width takes a `Vec2f`.
    if (Mat(2, 3, f32).vector_type != Vec2f) {
        @compileError("a two-row f32 matrix must apply itself to a Vec2f");
    }
    // A row count matching no named vector type falls back to a plain array,
    // which is where a 5 by 3 matrix lands.
    if (Mat(5, 3, f32).vector_type != [5]f32) {
        @compileError("an unmatched row count must use a plain array");
    }
}
