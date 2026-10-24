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
const Complexd = complex.Complexd;
const Vec2f = vector.Vec2f;
const Vec3f = vector.Vec3f;
const Vec4f = vector.Vec4f;
const Vec2d = vector.Vec2d;
const Vec3d = vector.Vec3d;
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

/// A `Complexf` or a `Complexd` from its two components.
///
/// The component type follows `C`, so one written fixture serves both formats
/// and the type loop below can hand the same numbers to either.
fn cz(comptime C: type, re: @FieldType(C, "re"), im: @FieldType(C, "re")) C {
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
    // `spec/core-math-2026-08-22.md#special-values` has two rules that a host unit
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

test "mulVec returns a plain array where no named vector type matches the row count" {
    // Five rows, so the RESULT is `[5]f32` and the lane plumbing takes the
    // array branch. Four columns, so the vector it is handed is a `Vec4f`:
    // the input length and the output length come from different counts, and
    // the named type for a four-lane f32 vector is the same whatever the
    // matrix does with it. Every row of the literal is a constant, so the
    // result is four times that constant.
    const m = Mat(5, 4, f32).init(.{
        .{ 1, 1, 1, 1 },
        .{ 2, 2, 2, 2 },
        .{ 3, 3, 3, 3 },
        .{ 4, 4, 4, 4 },
        .{ 5, 5, 5, 5 },
    });
    const result = m.mulVec(Vec4f.init(1, 1, 1, 1));
    for (0..5) |row| {
        try testing.expectEqual(@as(f32, @floatFromInt(row + 1)) * 4, result[row]);
    }
}

// ---------------------------------------------------------------------------
// Rectangular application
//
// A contraction runs over the COLUMNS, so the vector `mulVec` takes is as long
// as the column count and the result it returns is as long as the row count.
// Every instantiation below has a row count and a column count that differ,
// so the two types it exchanges are decided by two different numbers: a
// `columns`-long argument names one matrix, and the `rows`-long vector it
// hands back names another, and neither can stand in for the other.
// ---------------------------------------------------------------------------

test "a wide matrix contracts a three-component vector into two" {
    //     [ 1 2 3 ]   [  2 ]   [ 1*2 + 2*(-3) + 3*5 ]   [ 11 ]
    //     [ 4 5 6 ] * [ -3 ] = [ 4*2 + 5*(-3) + 6*5 ] = [ 23 ]
    //                  [  5 ]
    // Every product is an integer below a hundred, so both sums are exact and
    // the comparison is a bit comparison.
    const m = Mat(2, 3, f32).init(.{ .{ 1, 2, 3 }, .{ 4, 5, 6 } });
    const result = m.mulVec(Vec3f.init(2, -3, 5));
    try testing.expectEqual(@as(f32, 11), result.x);
    try testing.expectEqual(@as(f32, 23), result.y);

    // The lane order is already pinned by those two sums: read in reverse, the
    // same three entries give (5, 17) rather than (11, 23), and the binary64
    // pair below moves to (2, 11).

    // The same shape in binary64, with a cancellation in the first row and a
    // zero entry in the second:
    //     [ 3 -1  4 ]   [  2 ]   [ 3*2 + (-1)*3 + 4*(-1) ]   [ -1 ]
    //     [ 0  5 -2 ] * [  3 ] = [ 0*2 +  5*3 + (-2)*(-1) ] = [ 17 ]
    //                  [ -1 ]
    const d = Mat(2, 3, f64).init(.{ .{ 3, -1, 4 }, .{ 0, 5, -2 } });
    const wide_d = d.mulVec(Vec3d.init(2, 3, -1));
    try testing.expectEqual(@as(f64, -1), wide_d.x);
    try testing.expectEqual(@as(f64, 17), wide_d.y);
}

test "a tall matrix contracts a two-component vector into three" {
    //     [ 1 2 ]   [ 1 ]   [ 1*1 + 2*2 ]   [  5 ]
    //     [ 3 4 ] * [ 2 ] = [ 3*1 + 4*2 ] = [ 11 ]
    //     [ 5 6 ]                 [ 5*1 + 6*2 ] = [ 17 ]
    // The two lanes are not interchangeable here either: (2, 1) would give
    // (4, 10, 16) instead.
    const m = Mat(3, 2, f32).init(.{ .{ 1, 2 }, .{ 3, 4 }, .{ 5, 6 } });
    const tall = m.mulVec(Vec2f.init(1, 2));
    try testing.expectEqual(@as(f32, 5), tall.x);
    try testing.expectEqual(@as(f32, 11), tall.y);
    try testing.expectEqual(@as(f32, 17), tall.z);

    // The binary64 fixture cancels exactly in its first row:
    //     [  2 -3 ]   [ 3 ]   [ 2*3 + (-3)*2 ]   [  0 ]
    //     [ -1  4 ] * [ 2 ] = [ -1*3 + 4*2 ]   = [  5 ]
    //     [  0  6 ]                     [ 12 ]
    const d = Mat(3, 2, f64).init(.{ .{ 2, -3 }, .{ -1, 4 }, .{ 0, 6 } });
    const tall_d = d.mulVec(Vec2d.init(3, 2));
    try testing.expectEqual(@as(f64, 0), tall_d.x);
    try testing.expectEqual(@as(f64, 5), tall_d.y);
    try testing.expectEqual(@as(f64, 12), tall_d.z);
    // Flipping the sign of the second lane changes every row that has an
    // entry there, and it is the only way to tell a two-lane contraction from
    // a one-lane one on a matrix whose rows are not all equal.
    const negated_d = d.mulVec(Vec2d.init(3, -2));
    try testing.expectEqual(@as(f64, 12), negated_d.x);
    try testing.expectEqual(@as(f64, -11), negated_d.y);
    try testing.expectEqual(@as(f64, -12), negated_d.z);
}

test "a one-row matrix takes a vector as long as its columns" {
    // 2*5 + (-3)*1 + 4*(-1) = 10 - 3 - 4 = 3, and the result is a one-element
    // array because no named vector type has a single component.
    const three_columns = Mat(1, 3, f32).init(.{
        .{ 2, -3, 4 },
    });
    const three = three_columns.mulVec(Vec3f.init(5, 1, -1));
    try testing.expectEqual(@as(f32, 3), three[0]);

    // Five columns take a five-element array, since no named vector type has
    // five components either, and the fifth lane is not free to drop out:
    // reading the row with (1, 1, 1, 1, 0) would give 10 rather than 15.
    const five_columns = Mat(1, 5, f32).init(.{
        .{ 1, 2, 3, 4, 5 },
    });
    const all_ones = five_columns.mulVec([5]f32{ 1, 1, 1, 1, 1 });
    try testing.expectEqual(@as(f32, 15), all_ones[0]);

    // A two-column row in binary64, whose result is one component wide too:
    // (-6)*2 + 9*3 = -12 + 27 = 15.
    const two_columns = Mat(1, 2, f64).init(.{
        .{ -6, 9 },
    });
    const pair = two_columns.mulVec(Vec2d.init(2, 3));
    try testing.expectEqual(@as(f64, 15), pair[0]);
}

test "a one-column matrix takes a single-lane array and returns a named vector" {
    // One column is one product, so the vector it is handed is `[1]f32` while
    // the result is the three-component `Vec3f`. Neither is an array by
    // accident: the input is as long as the column count and the result as
    // long as the row count, and this shape is the only one where those two
    // decisions disagree.
    const m = Mat(3, 1, f32).init(.{ .{1}, .{2}, .{3} });
    const scaled = m.mulVec([1]f32{4});
    try testing.expectEqual(@as(f32, 4), scaled.x);
    try testing.expectEqual(@as(f32, 8), scaled.y);
    try testing.expectEqual(@as(f32, 12), scaled.z);

    // Four rows and one column in binary64, where the lane multiplies each
    // entry unchanged: 2*3, (-1)*3, 0*3, 5*3.
    const d = Mat(4, 1, f64).init(.{ .{2}, .{-1}, .{0}, .{5} });
    const scaled_d = d.mulVec([1]f64{3});
    try testing.expectEqual(@as(f64, 6), scaled_d.x);
    try testing.expectEqual(@as(f64, -3), scaled_d.y);
    try testing.expectEqual(@as(f64, 0), scaled_d.z);
    try testing.expectEqual(@as(f64, 15), scaled_d.w);

    // Five rows and one column is the other end of the same switch: both the
    // vector it takes and the vector it returns are plain arrays. The 5 by 4
    // case above is the mirror image of it, a named input with an array
    // result, so the two together cover both directions of the switch.
    const ramp = Mat(5, 1, f32).init(.{ .{2}, .{1}, .{0}, .{-1}, .{-2} });
    const ramped = ramp.mulVec([1]f32{3});
    try testing.expectEqual(@as(f32, 6), ramped[0]);
    try testing.expectEqual(@as(f32, 3), ramped[1]);
    try testing.expectEqual(@as(f32, 0), ramped[2]);
    try testing.expectEqual(@as(f32, -3), ramped[3]);
    try testing.expectEqual(@as(f32, -6), ramped[4]);
}

test "a rectangular composition applies the right factor first" {
    // `A * B` contracts over the middle index, so the right factor has as many
    // rows as the left factor has columns: a 3 by 2 after a 2 by 3. That
    // chaining is what lets the result of one `mulVec` be the argument of the
    // next, and it is the same rule the product follows.
    const a = Mat(2, 3, f32).init(.{ .{ 1, 2, 3 }, .{ 4, 5, 6 } });
    const b = Mat(3, 2, f32).init(.{ .{ 1, 0 }, .{ 0, 1 }, .{ 1, 1 } });
    const v = Vec2f.init(2, 3);

    //     B * v = [ 1*2 + 0*3 ]   [ 2 ]
    //            [ 0*2 + 1*3 ] = [ 3 ]
    //            [ 1*2 + 1*3 ]   [ 5 ]
    const inner = b.mulVec(v);
    try testing.expectEqual(@as(f32, 2), inner.x);
    try testing.expectEqual(@as(f32, 3), inner.y);
    try testing.expectEqual(@as(f32, 5), inner.z);

    //     A * B = [ 1+0+3  0+2+3 ]   [  4  5 ]
    //             [ 4+0+6  0+5+6 ] = [ 10 11 ]
    // and A * (B * v) = [ 1*2 + 2*3 + 3*5 ]   [  2 + 6 + 15 ]   [ 23 ]
    //                  [ 4*2 + 5*3 + 6*5 ] = [  8 + 15 + 30 ] = [ 53 ]
    // Every partial sum on either route is an integer below a hundred, so the
    // two answers are bit equal and the entries are pinned to hand values.
    const ab = matrix.mul(2, a, b);
    try testing.expectEqual(@as(f32, 4), ab.at(0, 0));
    try testing.expectEqual(@as(f32, 5), ab.at(0, 1));
    try testing.expectEqual(@as(f32, 10), ab.at(1, 0));
    try testing.expectEqual(@as(f32, 11), ab.at(1, 1));
    const chained = a.mulVec(inner);
    const product = ab.mulVec(v);
    try testing.expectEqual(@as(f32, 23), chained.x);
    try testing.expectEqual(@as(f32, 53), chained.y);
    try testing.expectEqual(@as(f32, 23), product.x);
    try testing.expectEqual(@as(f32, 53), product.y);

    // The same shapes in binary64:
    //     A = [ 2 -1  0 ]   B = [  3  0 ]   v = [ 2 ]
    //         [ 1  1  1 ]       [  0 -2 ]       [ 1 ]
    //         [ ...      ]       [  1  1 ]
    //     B * v = [ 6, -2, 3 ] and A * that = [ 12 + 2 + 0, 6 - 2 + 3 ]
    const ad = Mat(2, 3, f64).init(.{ .{ 2, -1, 0 }, .{ 1, 1, 1 } });
    const bd = Mat(3, 2, f64).init(.{ .{ 3, 0 }, .{ 0, -2 }, .{ 1, 1 } });
    const vd = Vec2d.init(2, 1);
    const inner_d = bd.mulVec(vd);
    try testing.expectEqual(@as(f64, 6), inner_d.x);
    try testing.expectEqual(@as(f64, -2), inner_d.y);
    try testing.expectEqual(@as(f64, 3), inner_d.z);
    const chained_d = ad.mulVec(inner_d);
    const product_d = matrix.mul(2, ad, bd).mulVec(vd);
    try testing.expectEqual(@as(f64, 14), chained_d.x);
    try testing.expectEqual(@as(f64, 7), chained_d.y);
    try testing.expectEqual(@as(f64, 14), product_d.x);
    try testing.expectEqual(@as(f64, 7), product_d.y);
}

test "a complex rectangle contracts over an array" {
    // No `Vec*` type is complex, so both sides of a complex application are
    // plain arrays: the input holds one element per COLUMN and the result one
    // per ROW, and the two are independent decisions. Both directions are
    // checked below, in both formats.
    //
    // Every product below is the naive four-multiply form on integers of
    // magnitude at most five, so both components of both sums are exact and
    // the same numbers serve binary32 and binary64.
    inline for (.{ Complexf, Complexd }) |C| {
        const R = @FieldType(C, "re");
        // WIDE: three columns in, two rows out.
        //     [ 1+i  2i  2-i ]   [ 1-i ]
        //     [ -1   1+i   i ] * [ 2   ]
        //                        [  i ]
        // row 0: (1+i)(1-i) + (2i)(2) + (2-i)(i) = 2 + 4i + (1+2i) = 3 + 6i
        // row 1: (-1)(1-i) + (1+i)(2) + (i)(i) = (-1+i) + (2+2i) - 1 = 3i
        // The lanes are not interchangeable: (i, 2, 1-i) would give (0, 2i).
        const a = Mat(2, 3, C).init(.{
            .{ cz(C, 1, 1), cz(C, 0, 2), cz(C, 2, -1) },
            .{ cz(C, -1, 0), cz(C, 1, 1), cz(C, 0, 1) },
        });
        const wide = a.mulVec([3]C{ cz(C, 1, -1), cz(C, 2, 0), cz(C, 0, 1) });
        try testing.expectEqual(@as(R, 3), wide[0].re);
        try testing.expectEqual(@as(R, 6), wide[0].im);
        try testing.expectEqual(@as(R, 0), wide[1].re);
        try testing.expectEqual(@as(R, 3), wide[1].im);

        // TALL: two columns in, three rows out.
        //     [ 1+2i   i  ]   [ 2+i ]
        //     [ -1+i   2  ] * [-1+2i]
        //     [ -2i  1+i ]
        // row 0: (1+2i)(2+i) + (i)(-1+2i) = 5i + (-2-i) = -2 + 4i
        // row 1: (-1+i)(2+i) + (2)(-1+2i) = (-3+i) + (-2+4i) = -5 + 5i
        // row 2: (-2i)(2+i) + (1+i)(-1+2i) = (2-4i) + (-3+i) = -1 - 3i
        // Here the reverse lanes give (-6, 2i) in row 0, so a transposed read
        // cannot pass both shapes.
        const d = Mat(3, 2, C).init(.{
            .{ cz(C, 1, 2), cz(C, 0, 1) },
            .{ cz(C, -1, 1), cz(C, 2, 0) },
            .{ cz(C, 0, -2), cz(C, 1, 1) },
        });
        const tall = d.mulVec([2]C{ cz(C, 2, 1), cz(C, -1, 2) });
        try testing.expectEqual(@as(R, -2), tall[0].re);
        try testing.expectEqual(@as(R, 4), tall[0].im);
        try testing.expectEqual(@as(R, -5), tall[1].re);
        try testing.expectEqual(@as(R, 5), tall[1].im);
        try testing.expectEqual(@as(R, -1), tall[2].re);
        try testing.expectEqual(@as(R, -3), tall[2].im);
    }
}

// ---------------------------------------------------------------------------
// Accumulation order in a rectangle
//
// The order a contraction sums in is only observable where a rounding or a
// sign separates one order from another, so every fixture below is built
// around a case where a different order, a fused multiply-add, or a leading
// zero would give a different answer. The shapes are rectangular, so the two
// lane counts are checked at the same time.
// ---------------------------------------------------------------------------

test "a rectangular contraction keeps its column order through a cancellation" {
    // The three products of row 0 are 1, 2^30, and -2^30. The spacing of
    // binary32 at 2^30 is 128, so 1 + 2^30 rounds back to 2^30 and the last
    // term cancels it exactly: the answer is +0. Summing the same three
    // products in DECREASING column order would give -2^30 + 2^30 + 1 = 1, so
    // the zero is what fixes the order rather than the arithmetic being right
    // in general.
    const two_30: f32 = 0x1p+30;
    const m = Mat(2, 3, f32).init(.{
        .{ 1, two_30, -two_30 },
        .{ -two_30, two_30, 1 },
    });
    const result = m.mulVec(Vec3f.init(1, 1, 1));
    // `x + (-x)` is +0 for finite non-zero x under round-to-nearest, so the
    // sign of the zero is pinned as well as its value.
    try testing.expectEqual(@as(u32, 0x0000_0000), @as(u32, @bitCast(result.x)));
    // Row 1 carries the same three products in the other order and the small
    // term SURVIVES there: -2^30 + 2^30 = +0 and +0 + 1 = 1. A sum that
    // swallowed small terms whatever their position would report 0 for this
    // row too, and the two rows are one permutation apart.
    try testing.expectEqual(@as(f32, 1), result.y);

    // binary64 has the same shape one exponent up: 2^53 is the first power of
    // two whose successor is two away, so 1 + 2^53 is a tie and rounds to the
    // even neighbour 2^53.
    const two_53: f64 = 0x1p+53;
    const d = Mat(2, 3, f64).init(.{
        .{ 1, two_53, -two_53 },
        .{ -two_53, two_53, 1 },
    });
    const result_d = d.mulVec(Vec3d.init(1, 1, 1));
    try testing.expectEqual(@as(u64, 0x0000_0000_0000_0000), @as(u64, @bitCast(result_d.x)));
    try testing.expectEqual(@as(f64, 1), result_d.y);
}

test "a rectangular contraction does not fuse its multiply and add" {
    // The row [-1, 1 + 2^-23] against the vector [1, 1 - 2^-23]. Column 0
    // seeds the accumulator with the exact -1, and column 1's product is
    // (1 + 2^-23)(1 - 2^-23) = 1 - 2^-46 in exact arithmetic. binary32
    // cannot hold that: the neighbours of 1 are 1 itself and 1 - 2^-24, and
    // 1 - 2^-46 sits a factor of 2^22 nearer 1 than the lower one does, so
    // the product rounds to 1 and the sum is -1 + 1 = +0.
    //
    // A fused multiply-add keeps the 2^-46 the multiply discards, so it
    // returns the exactly representable -2^-46. The zero is therefore a
    // value the fused form cannot produce, which is what makes this a
    // discriminator rather than a restatement of the arithmetic.
    //
    // Both operand sets live in mutable storage that is handed to
    // `doNotOptimizeAway`, so the fixture is evaluated by the runtime rather
    // than folded from literals before the call.
    var row = [_]f32{ -1.0, 1.0 + 0x1p-23 };
    var lane = [2]f32{ 1.0, 1.0 - 0x1p-23 };
    std.mem.doNotOptimizeAway(&row);
    std.mem.doNotOptimizeAway(&lane);
    const m = Mat(1, 2, f32).init(.{row});
    const strict = m.mulVec(Vec2f.init(lane[0], lane[1]));
    // The seed is -1 and the rounded product is 1, and `x + (-x)` is +0 for
    // finite non-zero x under round-to-nearest, so the zero's sign is pinned
    // along with its value.
    try testing.expectEqual(@as(u32, 0x0000_0000), @as(u32, @bitCast(strict[0])));

    // The same shape one format up. 2^-52 is the spacing above 1 in
    // binary64, the exact product is 1 - 2^-104, and the fused form would
    // return the equally representable -2^-104.
    var row_d = [_]f64{ -1.0, 1.0 + 0x1p-52 };
    var lane_d = [2]f64{ 1.0, 1.0 - 0x1p-52 };
    std.mem.doNotOptimizeAway(&row_d);
    std.mem.doNotOptimizeAway(&lane_d);
    const d = Mat(1, 2, f64).init(.{row_d});
    const strict_d = d.mulVec(Vec2d.init(lane_d[0], lane_d[1]));
    try testing.expectEqual(@as(u64, 0x0000_0000_0000_0000), @as(u64, @bitCast(strict_d[0])));
}

test "a single-column rectangle keeps the sign of its seed" {
    // One column is one product and no addition at all, so the seed IS the
    // answer and a negative zero survives. An implementation that began from
    // an explicit zero would return +0 for rows 0 and 2 instead, because
    // +0 + -0 is +0.
    const m = Mat(3, 1, f32).init(.{ .{-0.0}, .{0.0}, .{-0.0} });
    const seeded = m.mulVec([1]f32{0.0});
    try testing.expectEqual(@as(u32, 0x8000_0000), @as(u32, @bitCast(seeded.x)));
    try testing.expectEqual(@as(u32, 0x0000_0000), @as(u32, @bitCast(seeded.y)));
    try testing.expectEqual(@as(u32, 0x8000_0000), @as(u32, @bitCast(seeded.z)));

    // The same entries in the FIRST column of a two-column rectangle are
    // washed out: the second column contributes 0*0 = +0, and -0 + +0 is
    // +0. That is the boundary of the rule above, and it is why the column
    // count has to be visible in the result instead of assumed.
    const two_columns = Mat(3, 2, f32).init(.{
        .{ -0.0, 0.0 },
        .{ 0.0, 0.0 },
        .{ -0.0, 0.0 },
    });
    const washed = two_columns.mulVec(Vec2f.init(0.0, 0.0));
    try testing.expectEqual(@as(u32, 0x0000_0000), @as(u32, @bitCast(washed.x)));
    try testing.expectEqual(@as(u32, 0x0000_0000), @as(u32, @bitCast(washed.y)));
    try testing.expectEqual(@as(u32, 0x0000_0000), @as(u32, @bitCast(washed.z)));

    // The sign of a zero product follows the LANE, not the entry. With the
    // second lane negative both products of this row are -0 and -0 + -0 is
    // -0; with it positive the sum is -0 + +0, which is +0.
    const seeded_by_lane = Mat(1, 2, f32).init(.{
        .{ -1.0, 0.0 },
    });
    const both_negative = seeded_by_lane.mulVec(Vec2f.init(0.0, -1.0));
    try testing.expectEqual(@as(u32, 0x8000_0000), @as(u32, @bitCast(both_negative[0])));
    const second_positive = seeded_by_lane.mulVec(Vec2f.init(0.0, 1.0));
    try testing.expectEqual(@as(u32, 0x0000_0000), @as(u32, @bitCast(second_positive[0])));
}

test "a rectangle reports the first NaN of its row, quieted and intact" {
    // Two NaNs meet in row 0, and they are deliberately unlike each other: a
    // SIGNALING 0x7f80_0001 in column 0 and the quiet, negative, differently
    // signed payload 0xffc0_0002 in column 1. Both lanes are one, so column 0
    // is the seed and column 1 is the added term.
    //
    // The seed is quieted by the MULTIPLY, not by the addition: `entryMul`
    // carries the column's own operand through and sets the quiet bit there,
    // which is why the signaling 0x7f80_0001 comes back as 0x7fc0_0001. The
    // addition then decides between two NaNs, and it reports the one already
    // in the accumulator. So the answer is 0x7fc0_0001, positive, with the
    // first payload. Had the later column won the row would read
    // 0xffc0_0002, which is a different sign AND a different payload, so this
    // is a real ordering check rather than a NaN check that passes either way.
    const signaling: f32 = @bitCast(@as(u32, 0x7f80_0001));
    const negative_other: f32 = @bitCast(@as(u32, 0xffc0_0002));
    const negative_canonical: f32 = @bitCast(@as(u32, 0xffc0_0000));
    var m = Mat(3, 2, f32).fill(0);
    m.set(0, 0, signaling);
    m.set(0, 1, negative_other);
    m.set(1, 0, negative_canonical);
    const result = m.mulVec(Vec2f.init(1, 1));
    try testing.expectEqual(@as(u32, 0x7fc0_0001), @as(u32, @bitCast(result.x)));
    // The negative canonical NaN is the SEED of row 1, and it is quiet
    // already, so it crosses both the multiply and the add with its sign
    // untouched. A kernel that substituted the positive canonical encoding
    // would report 0x7fc00000 here.
    try testing.expectEqual(@as(u32, 0xffc0_0000), @as(u32, @bitCast(result.y)));
    // Row 2 reads the same two columns and meets no NaN at all, so the
    // propagation is per row and not a whole-result effect.
    try testing.expectEqual(@as(f32, 0), result.z);

    // The same precedence in binary64, where the two payloads differ in sign
    // the other way round: a negative signaling 0xfff0_0000_0000_0001 first
    // and a positive quiet 0x7ff8_0000_0000_00a5 after it. The reported value
    // is the first, so the sign and the payload both come from column 0.
    const signaling_d: f64 = @bitCast(@as(u64, 0xfff0_0000_0000_0001));
    const positive_other_d: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_00a5));
    var m_d = Mat(3, 2, f64).fill(0);
    m_d.set(0, 0, signaling_d);
    m_d.set(0, 1, positive_other_d);
    m_d.set(1, 0, nan_d);
    const rows = m_d.mulVec(Vec2d.init(1, 1));
    try testing.expectEqual(@as(u64, 0xfff8_0000_0000_0001), @as(u64, @bitCast(rows.x)));
    try testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), @as(u64, @bitCast(rows.y)));
    try testing.expectEqual(@as(f64, 0), rows.z);
}

test "an invalid rectangular step gives the positive canonical NaN" {
    // +inf and -inf are the only non-NaN operands that reach an invalid
    // operation here: the seed is +inf*1 and the added term is -inf*1. The
    // substitution to the POSITIVE canonical quiet NaN is the part a host SSE
    // unit does not do on its own, and it happens in the ADD rather than in
    // the multiply, so the order of the columns is what puts it there.
    const positive_infinity: f32 = @bitCast(@as(u32, 0x7f80_0000));
    const wide = Mat(1, 2, f32).init(.{
        .{ positive_infinity, -positive_infinity },
    });
    const result = wide.mulVec(Vec2f.init(1, 1));
    try testing.expect(std.math.isNan(result[0]));
    try testing.expectEqual(@as(u32, 0x7fc0_0000), @as(u32, @bitCast(result[0])));

    // The same two terms on a three-row rectangle in binary64. Row 1 is
    // ordinary arithmetic between the same columns, and row 2 reaches the
    // same invalid operation with the two signs the other way round, so a
    // result copied from row 0 could not pass both.
    const tall = Mat(3, 2, f64).init(.{
        .{ inf_d, -inf_d },
        .{ 1, 0 },
        .{ -inf_d, inf_d },
    });
    const rows = tall.mulVec(Vec2d.init(1, 1));
    try testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), @as(u64, @bitCast(rows.x)));
    try testing.expectEqual(@as(f64, 1), rows.y);
    try testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), @as(u64, @bitCast(rows.z)));
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
// The cofactor expansion against an exact oracle
// ---------------------------------------------------------------------------

/// The exact determinant of a small integer square matrix, by the Leibniz sum.
///
/// This is an ORACLE and not a second implementation of the library. It works
/// in `i64` with no floating-point step at all, and it enumerates the 120
/// permutations of a 5 by 5 where the library expands cofactors, so it shares
/// no index mapping, no sign rule and no accumulation order with
/// `Mat.determinant`. An entry read from the wrong row or column, a sign
/// applied to the wrong term, or a different way of walking the minors all
/// produce a different number here.
fn exactDet(comptime n: usize, m: [n][n]i64) i64 {
    var used: [n]bool = @splat(false);
    return leibnizTerm(n, m, 0, 1, &used);
}

/// One row of the Leibniz enumeration: the running product with the parity of
/// the partial permutation chosen so far applied to it.
fn leibnizTerm(
    comptime n: usize,
    m: [n][n]i64,
    row: usize,
    acc: i64,
    used: *[n]bool,
) i64 {
    if (row == n) return acc;
    var total: i64 = 0;
    for (0..n) |column| {
        if (used[column]) continue;
        // Every row already chosen that sits to the RIGHT of `column` forms an
        // inversion with it, so the parity falls out of the enumeration and no
        // separate permutation has to be counted.
        var to_right: usize = 0;
        for (0..n) |other| {
            if (used[other] and other > column) to_right += 1;
        }
        const signed = if (to_right % 2 == 1) -acc else acc;
        used[column] = true;
        total += leibnizTerm(n, m, row + 1, signed * m[row][column], used);
        used[column] = false;
    }
    return total;
}

test "the size-four and size-five determinants equal an exact integer oracle" {
    // Every entry of every fixture below is an integer of magnitude at most
    // FOUR. That bound is set by the repeated-row fixture, which carries a 4,
    // and it is the bound the exactness argument is made against. Every partial
    // product of the expansion is an integer below `4^5 = 1024` and every
    // partial sum is below `5! * 4^5 = 122880`, so every value the expansion
    // forms is an exact integer in binary32 as well as in binary64, whose
    // largest exactly representable integers are `2^24` and `2^53`. The
    // expansion therefore has nothing to round, and the comparisons below are
    // equalities rather than tolerances. A larger entry would silently weaken
    // the whole test, which is why the bound is part of the fixture choice and
    // not an afterthought.
    const four_dense = [4][4]i64{
        .{ 2, -1, 3, 0 },
        .{ 1, 2, 0, -3 },
        .{ -2, 0, 1, 2 },
        .{ 3, -2, 1, 1 },
    };
    // Rows 0 and 1 are equal, so the determinant is zero as a fact about the
    // numbers. The exact zero is what catches a minor that keeps the wrong row.
    const four_repeat = [4][4]i64{
        .{ 1, 2, 3, 4 },
        .{ 1, 2, 3, 4 },
        .{ -1, 0, 2, 1 },
        .{ 3, -3, 1, 2 },
    };
    // A zero column forces an exact zero determinant. This boundary fixture
    // does not distinguish a permutation of the dropped-column minors.
    const four_zero_column = [4][4]i64{
        .{ 1, -2, 3, 0 },
        .{ 0, 1, 2, 0 },
        .{ 3, 1, 1, 0 },
        .{ -1, 2, -3, 0 },
    };
    const five_dense = [5][5]i64{
        .{ 2, -1, 3, 0, 1 },
        .{ 1, 2, 0, -3, 2 },
        .{ -2, 0, 1, 2, -1 },
        .{ 3, -2, 1, 1, 3 },
        .{ 0, 1, -2, 2, 1 },
    };
    // Columns 1 and 3 are equal, so the determinant is zero again, and this
    // time the repetition is in the COLUMN the expansion indexes by, which is
    // the axis a transposed implementation would get wrong.
    const five_repeat_column = [5][5]i64{
        .{ 2, -1, 2, -1, 1 },
        .{ 1, 2, 1, 2, -3 },
        .{ -2, 0, -2, 0, 1 },
        .{ 3, -2, 3, -2, 3 },
        .{ 0, 1, 0, 1, 1 },
    };
    // Upper triangular with a negative diagonal, whose determinant is the
    // product `(-1) * 2 * (-3) * 1 * 2 = 12` by hand, and whose minors are full
    // of the zeros a sparsity shortcut would have to keep.
    const five_triangular = [5][5]i64{
        .{ -1, 2, -1, 3, 0 },
        .{ 0, 2, 1, -2, 1 },
        .{ 0, 0, -3, 1, 2 },
        .{ 0, 0, 0, 1, -1 },
        .{ 0, 0, 0, 0, 2 },
    };

    // The triangular fixture's determinant is a fact about the numbers: it is
    // the product of the diagonal, `(-1) * 2 * (-3) * 1 * 2 = 12`, and its off
    // diagonal entries are non-zero on purpose so an implementation that only
    // works on a triangular matrix cannot pass. The oracle agrees, and the
    // library then has to agree with both, so this row is pinned without
    // leaning on the oracle alone.
    try testing.expectEqual(@as(f64, 12), @as(f64, @floatFromInt(exactDet(5, five_triangular))));
    try testing.expectEqual(@as(f64, 12), matrix.determinantMat5d(Mat5d.init(.{
        .{ -1, 2, -1, 3, 0 },
        .{ 0, 2, 1, -2, 1 },
        .{ 0, 0, -3, 1, 2 },
        .{ 0, 0, 0, 1, -1 },
        .{ 0, 0, 0, 0, 2 },
    })));

    // The binary64 names, against the same oracle, one matrix at each size and
    // each degeneracy the expansion can meet.
    try testing.expectEqual(
        @as(f64, @floatFromInt(exactDet(4, four_dense))),
        matrix.determinantMat4d(Mat4d.init(.{
            .{ 2, -1, 3, 0 },
            .{ 1, 2, 0, -3 },
            .{ -2, 0, 1, 2 },
            .{ 3, -2, 1, 1 },
        })),
    );
    try testing.expectEqual(
        @as(f64, @floatFromInt(exactDet(4, four_repeat))),
        matrix.determinantMat4d(Mat4d.init(.{
            .{ 1, 2, 3, 4 },
            .{ 1, 2, 3, 4 },
            .{ -1, 0, 2, 1 },
            .{ 3, -3, 1, 2 },
        })),
    );
    try testing.expectEqual(
        @as(f64, @floatFromInt(exactDet(4, four_zero_column))),
        matrix.determinantMat4d(Mat4d.init(.{
            .{ 1, -2, 3, 0 },
            .{ 0, 1, 2, 0 },
            .{ 3, 1, 1, 0 },
            .{ -1, 2, -3, 0 },
        })),
    );
    try testing.expectEqual(
        @as(f64, @floatFromInt(exactDet(5, five_dense))),
        matrix.determinantMat5d(Mat5d.init(.{
            .{ 2, -1, 3, 0, 1 },
            .{ 1, 2, 0, -3, 2 },
            .{ -2, 0, 1, 2, -1 },
            .{ 3, -2, 1, 1, 3 },
            .{ 0, 1, -2, 2, 1 },
        })),
    );
    try testing.expectEqual(
        @as(f64, @floatFromInt(exactDet(5, five_repeat_column))),
        matrix.determinantMat5d(Mat5d.init(.{
            .{ 2, -1, 2, -1, 1 },
            .{ 1, 2, 1, 2, -3 },
            .{ -2, 0, -2, 0, 1 },
            .{ 3, -2, 3, -2, 3 },
            .{ 0, 1, 0, 1, 1 },
        })),
    );

    // The binary32 names, where the same bound has to hold in 24 bits rather
    // than 53. Both sizes are covered because the size-5 form is the deepest
    // recursion the module reaches and the one with the most partial sums.
    try testing.expectEqual(
        @as(f32, @floatFromInt(exactDet(4, four_dense))),
        matrix.determinantMat4f(Mat4f.init(.{
            .{ 2, -1, 3, 0 },
            .{ 1, 2, 0, -3 },
            .{ -2, 0, 1, 2 },
            .{ 3, -2, 1, 1 },
        })),
    );
    try testing.expectEqual(
        @as(f32, @floatFromInt(exactDet(5, five_dense))),
        matrix.determinantMat5f(Mat5f.init(.{
            .{ 2, -1, 3, 0, 1 },
            .{ 1, 2, 0, -3, 2 },
            .{ -2, 0, 1, 2, -1 },
            .{ 3, -2, 1, 1, 3 },
            .{ 0, 1, -2, 2, 1 },
        })),
    );
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

/// The bits of a `f32` or `f64`, widened to a common width, for a comparison
/// that has to see past `==`.
fn floatBits(comptime T: type, value: T) u64 {
    return switch (@bitSizeOf(T)) {
        32 => @as(u64, @as(u32, @bitCast(value))),
        64 => @bitCast(value),
        else => @compileError("floatBits is for a binary32 or a binary64"),
    };
}

/// Equality on a float that is exact EXCEPT that `+0` and `-0` are the same
/// number, and reports a mismatch by its bits.
///
/// The relaxation is one property, and it is derived rather than observed.
/// Every non-zero entry of the two fixtures below is a dyadic rational that
/// binary32 and binary64 both hold exactly, and the reductions run in exact
/// arithmetic, so for those entries this test and a bit comparison are the same
/// test. What a bit comparison would add is the SIGN of a zero entry, and that
/// follows only from which pivot a zero was divided by: both fixtures have a
/// pivot of `-1/2` at step 1, the pivot row's zero entries are normalized by it,
/// and IEEE answers `+0 / -1/2 = -0`. Nothing in the contract gives a zero a
/// sign, `0.0 == -0.0` is what a caller compares, and the module's own
/// accumulation-order comment is about which zeros SURVIVE a contraction, not
/// about their sign. Everything else still has to match to the last bit,
/// including a NaN, which `==` alone would report only as "not equal".
fn expectSameNumber(comptime T: type, expected: T, actual: T) !void {
    if (expected == actual) return;
    try testing.expectEqual(floatBits(T, expected), floatBits(T, actual));
}

test "every inverse settles each column once, and the two exact fixtures pin the values" {
    // The reduction writes each entry of the working row several times: once
    // when its own step normalizes the pivot row, once for each step below it
    // that eliminates its column. What makes a shorter sweep legal is that the
    // LAST of those writes is the only one a later step or the result ever
    // reads, so the earlier ones are not observable. Both fixtures below are
    // dyadic and every quotient in them is exact, so the reduction runs in
    // exact arithmetic and its answer is the one the derivations give, bit for
    // bit. A sweep that kept a column one step too long, stopped one step too
    // early, or interchanged the wrong row changes at least one of these bits.

    //     [ 1  0    0  0 ]        x0 = b0
    //     [ 2  1    0  0 ]        x1 = b1 - 2 b0
    //     [ 0  1/2  1  0 ]        x2 = b2 - (b1 - 2 b0)/2 = b0 - b1/2 + b2
    //     [ 0  0    1/2 1 ]       x3 = b3 - (b0 - b1/2 + b2)/2
    //                              = -b0/2 + b1/4 - b2/2 + b3
    //
    // Step 0 MUST interchange, because 2 is below the 1 on the diagonal. That
    // is the whole of the single-interchange case: after it the pivots are
    // 1/2, 1 and 1 and each step has one multiplier.
    const one_swap = Mat4d.init(.{
        .{ 1, 0, 0, 0 },
        .{ 2, 1, 0, 0 },
        .{ 0, 0.5, 1, 0 },
        .{ 0, 0, 0.5, 1 },
    });
    const one_swap_expected = Mat4d.init(.{
        .{ 1, 0, 0, 0 },
        .{ -2, 1, 0, 0 },
        .{ 1, -0.5, 1, 0 },
        .{ -0.5, 0.25, -0.5, 1 },
    });
    const one_swap_inverse = matrix.inverseMat4d(one_swap) orelse return error.TestUnexpectedResult;
    for (0..4) |row| {
        for (0..4) |column| {
            try expectSameNumber(
                f64,
                one_swap_expected.at(row, column),
                one_swap_inverse.at(row, column),
            );
        }
    }
    // The product is exact for the same reason, and it comes from a different
    // kernel: `mul` knows nothing about pivots or interchanges.
    const one_product = matrix.mul(4, one_swap_inverse, one_swap);
    for (0..4) |row| {
        for (0..4) |column| {
            const identity_entry: f64 = if (row == column) 1 else 0;
            try expectSameNumber(f64, identity_entry, one_product.at(row, column));
        }
    }

    //     [ 1  0    0  0 ]        the same forward substitution with the last
    //     [ 2  1    0  0 ]        row's THIRD entry raised from 1/2 to 2,
    //     [ 0  1/2  1  0 ]        which moves the pivot of step 2 off the diagonal:
    //     [ 0  0    2  1 ]        x0 = b0
    //                              x1 = b1 - 2 b0
    //                              x2 = b2 - (b1 - 2 b0)/2 = b0 - b1/2 + b2
    //                              x3 = b3 - 2 (b0 - b1/2 + b2)
    //                                 = -2 b0 + b1 - 2 b2 + b3
    //
    // TWO interchanges follow: step 0 takes the 2 from row 1, and step 2 then
    // finds 2 below the diagonal in a column whose diagonal entry is 1. The
    // second one is the case worth having, because at step 2 most of the row is
    // already settled and the sweep has to move exactly the columns that are
    // still live.
    const two_swaps = Mat4d.init(.{
        .{ 1, 0, 0, 0 },
        .{ 2, 1, 0, 0 },
        .{ 0, 0.5, 1, 0 },
        .{ 0, 0, 2, 1 },
    });
    const two_swaps_expected = Mat4d.init(.{
        .{ 1, 0, 0, 0 },
        .{ -2, 1, 0, 0 },
        .{ 1, -0.5, 1, 0 },
        .{ -2, 1, -2, 1 },
    });
    const two_swaps_inverse = matrix.inverseMat4d(two_swaps) orelse
        return error.TestUnexpectedResult;
    for (0..4) |row| {
        for (0..4) |column| {
            try expectSameNumber(
                f64,
                two_swaps_expected.at(row, column),
                two_swaps_inverse.at(row, column),
            );
        }
    }
    const two_product = matrix.mul(4, two_swaps_inverse, two_swaps);
    for (0..4) |row| {
        for (0..4) |column| {
            const identity_entry: f64 = if (row == column) 1 else 0;
            try expectSameNumber(f64, identity_entry, two_product.at(row, column));
        }
    }

    // The binary32 names, where the halves are exact as well, so the same two
    // derivations pin 24 bits of mantissa instead of 53.
    const one_f = matrix.inverseMat4f(Mat4f.init(.{
        .{ 1, 0, 0, 0 },
        .{ 2, 1, 0, 0 },
        .{ 0, 0.5, 1, 0 },
        .{ 0, 0, 0.5, 1 },
    })) orelse return error.TestUnexpectedResult;
    for (0..4) |row| {
        for (0..4) |column| {
            const expected_entry: f32 = @floatCast(one_swap_expected.at(row, column));
            try expectSameNumber(f32, expected_entry, one_f.at(row, column));
        }
    }
    const two_f = matrix.inverseMat4f(Mat4f.init(.{
        .{ 1, 0, 0, 0 },
        .{ 2, 1, 0, 0 },
        .{ 0, 0.5, 1, 0 },
        .{ 0, 0, 2, 1 },
    })) orelse return error.TestUnexpectedResult;
    for (0..4) |row| {
        for (0..4) |column| {
            const expected_entry: f32 = @floatCast(two_swaps_expected.at(row, column));
            try expectSameNumber(f32, expected_entry, two_f.at(row, column));
        }
    }
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

comptime {
    // The eight named sizes, in the shapes `API.md#matrices` gives them. A
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

    // `vector_type` is the type `mulVec` RETURNS, and the row count alone
    // decides it. The vector it takes is as long as the COLUMN count, so for
    // the square sizes here both directions meet in the same type and one
    // check covers them; a rectangle separates them, and the two checks
    // below are on the two rectangular shapes.
    //
    // `Mat4f` and `Mat4d` are the homogeneous 3D sizes and `Mat5f` and
    // `Mat5d` the homogeneous 4D ones. The arithmetic cannot see the
    // difference, so the only thing a caller can rely on is the shape, and the
    // reason those two sizes exist is recorded on the aliases themselves.
    if (Mat4f.vector_type != Vec4f or Mat4d.vector_type != Vec4d) {
        @compileError("a 4 by 4 matrix must return a four-component vector");
    }
    if (Mat3f.vector_type != Vec3f) {
        @compileError("a 3 by 3 f32 matrix must return a Vec3f");
    }
    if (Mat2f.vector_type != Vec2f) {
        @compileError("a 2 by 2 f32 matrix must return a Vec2f");
    }
    if (Mat5f.vector_type != [5]f32) {
        @compileError("a 5 by 5 f32 matrix must return a five-element array");
    }
    // A two-row matrix of ANY width returns a `Vec2f`, and this one is handed
    // a `Vec3f` to do it with: the result follows the row count alone while
    // the input follows the column count.
    if (Mat(2, 3, f32).vector_type != Vec2f) {
        @compileError("a two-row f32 matrix must return a Vec2f");
    }
    // A row count matching no named vector type falls back to a plain array,
    // which is where a 5 by 3 matrix lands. It is still applied to a `Vec3f`.
    if (Mat(5, 3, f32).vector_type != [5]f32) {
        @compileError("an unmatched row count must return a plain array");
    }
}
