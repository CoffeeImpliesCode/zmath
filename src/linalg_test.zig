//! Tests for the `zmath.linalg` surface.
//!
//! Every expected value here is an identity or a hand derivation, never a
//! pinned bit pattern that was not derived. Two families appear, and each
//! carries the reason it is checked the way it is.
//!
//! EXACT comparisons, on operands that are small integers. Every product,
//! every partial sum, and every factor is an integer of magnitude far below
//! 2^24 in binary32 and 2^53 in binary64, so the operations are exact and a
//! bit comparison is the right comparison. This is what catches a transposed
//! index or a transposed sign.
//!
//! TOLERANT comparisons, on operands that no binary format holds exactly. A
//! factorization does divisions, so `reconstructLu` of a well-conditioned
//! matrix is not its input bit for bit, and every tolerance below is stated
//! with the magnitude actually observed.
//!
//! The two central identities are `reconstructLu(lu(A).?) == A` and
//! `reconstructCholesky(cholesky(A).?) == A`: they are what prove a
//! factorization is faithful. The determinant tests CROSS-CHECK against
//! `zmath.matrix`, which computes a determinant by cofactor expansion and is
//! an independent computation of the same quantity. The solve tests are the
//! standard construct-and-solve round trip, and they test the SOLVE, not the
//! factorization: `reconstructLu` already covers that, so a round trip alone
//! would be satisfied by a factorization and a solve that make the same
//! mistake twice.

const std = @import("std");
const testing = std.testing;
const linalg = @import("linalg.zig");
const matrix = @import("matrix.zig");
const vector = @import("vector.zig");

const Mat = matrix.Mat;
const Mat2f = matrix.Mat2f;
const Mat3f = matrix.Mat3f;
const Mat4f = matrix.Mat4f;
const Mat2d = matrix.Mat2d;
const Mat3d = matrix.Mat3d;
const Mat4d = matrix.Mat4d;
const Vec2f = vector.Vec2f;
const Vec3f = vector.Vec3f;
const Vec4f = vector.Vec4f;
const Vec2d = vector.Vec2d;
const Vec3d = vector.Vec3d;
const Vec4d = vector.Vec4d;

const LuMat2f = linalg.LuMat2f;
const LuMat3f = linalg.LuMat3f;
const LuMat4f = linalg.LuMat4f;
const LuMat2d = linalg.LuMat2d;
const LuMat3d = linalg.LuMat3d;
const LuMat4d = linalg.LuMat4d;
const CholeskyMat2f = linalg.CholeskyMat2f;
const CholeskyMat3f = linalg.CholeskyMat3f;
const CholeskyMat4f = linalg.CholeskyMat4f;
const CholeskyMat2d = linalg.CholeskyMat2d;
const CholeskyMat3d = linalg.CholeskyMat3d;
const CholeskyMat4d = linalg.CholeskyMat4d;

/// The positive canonical quiet NaN of binary32, from the special-value table
/// of `spec/core-math-2026-08-22.md:769-801`.
const nan_f: f32 = @bitCast(@as(u32, 0x7fc0_0000));
/// The positive infinity of binary32, from the same table.
const inf_f: f32 = @bitCast(@as(u32, 0x7f80_0000));
/// The positive canonical quiet NaN of binary64.
const nan_d: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_0000));
/// The positive infinity of binary64.
const inf_d: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0000));

/// Every entry of `actual` within `tolerance` of the matching entry of
/// `expected`.
fn expectClose(comptime M: type, expected: M, actual: M, tolerance: M.Element) !void {
    const n = comptime M.rows;
    inline for (0..n) |row| {
        inline for (0..n) |column| {
            try testing.expectApproxEqAbs(
                expected.data[column][row],
                actual.data[column][row],
                tolerance,
            );
        }
    }
}

/// Every entry of `actual` EQUAL to the matching entry of `expected`, which is
/// the right comparison only when the operands are exact.
fn expectExactly(comptime M: type, expected: M, actual: M) !void {
    const n = comptime M.rows;
    inline for (0..n) |row| {
        inline for (0..n) |column| {
            try testing.expectEqual(
                expected.data[column][row],
                actual.data[column][row],
            );
        }
    }
}

/// The positive canonical quiet NaN or the positive infinity of `M`'s entry
/// type, chosen by comptime `want_nan`.
fn specialOf(comptime M: type, comptime want_nan: bool) M.Element {
    const T = M.Element;
    if (T == f32) return if (want_nan) nan_f else inf_f;
    return if (want_nan) nan_d else inf_d;
}

// ---------------------------------------------------------------------------
// The matrices under test
//
// Every one is either an exactly representable small integer matrix or a
// derivation in the comment beside it. The symmetric ones are symmetric BIT for
// bit, which is what the Cholesky asymmetry policy asks for.
// ---------------------------------------------------------------------------

/// A well-conditioned 2 by 2 in small integers: determinant `4 * 3 - 3 * 6 =
/// -6`. Its first column is `[4, 6]`, so partial pivoting interchanges rows 0
/// and 1, and one interchange is an odd parity.
const well2f = Mat2f.init(.{ .{ 4, 3 }, .{ 6, 3 } });

/// A well-conditioned 3 by 3 in small integers: determinant
/// `2 * (9 - 4) - 1 * (3 - 2) + 1 * (2 - 3) = 10 - 1 - 1 = 8`. Its first
/// column is `[2, 1, 1]`, whose largest entry is already the pivot, and so is
/// every later one: the reduction needs no interchange.
const well3f = Mat3f.init(.{ .{ 2, 1, 1 }, .{ 1, 3, 2 }, .{ 1, 2, 3 } });

/// A well-conditioned 4 by 4 in small integers: determinant 12 by cofactor
/// expansion, with pivots `2, 2.5, 1.6, 1.5`, every one of them already in
/// place.
const well4f = Mat4f.init(.{
    .{ 2, 1, 1, 0 },
    .{ 1, 3, 2, 1 },
    .{ 1, 2, 3, 1 },
    .{ 0, 1, 1, 2 },
});

/// A symmetric positive definite 2 by 2 in small integers. Sylvester's
/// criterion: the leading principal minors are `4` and `4 * 3 - 1 = 11`, both
/// positive.
const spd2f = Mat2f.init(.{ .{ 4, 1 }, .{ 1, 3 } });

/// A symmetric positive definite 3 by 3 in small integers, leading principal
/// minors `4`, `19`, `97`, all positive.
const spd3f = Mat3f.init(.{ .{ 4, 1, 1 }, .{ 1, 5, 2 }, .{ 1, 2, 6 } });

/// A symmetric positive definite 4 by 4 in small integers, leading principal
/// minors `4`, `19`, `97`, `608`, all positive.
const spd4f = Mat4f.init(.{
    .{ 4, 1, 1, 0 },
    .{ 1, 5, 2, 1 },
    .{ 1, 2, 6, 2 },
    .{ 0, 1, 2, 7 },
});

/// The identity with its first two rows interchanged, in 2, 3, and 4 sizes.
///
/// Column 0 is `[0, 1, ...]`, so the largest magnitude below the diagonal is 1
/// and step 0 MUST interchange. After that the matrix IS the identity and no
/// later step moves a row, so exactly one interchange follows and the parity is
/// `-1`. Every entry is 0 or 1, so the whole factorization is exact: the
/// factors are the identity, and the reconstruction is the transposition
/// itself. A bit comparison is therefore the right comparison here, and it is
/// what catches a permutation applied to the wrong axis.
const swap2f = Mat2f.init(.{ .{ 0, 1 }, .{ 1, 0 } });
const swap3f = Mat3f.init(.{ .{ 0, 1, 0 }, .{ 1, 0, 0 }, .{ 0, 0, 1 } });
const swap4f = Mat4f.init(.{
    .{ 0, 1, 0, 0 },
    .{ 1, 0, 0, 0 },
    .{ 0, 0, 1, 0 },
    .{ 0, 0, 0, 1 },
});

/// The 4 by 4 matrix where partial pivoting genuinely matters.
///
///     [ e  1  0  0 ]      e = 4e-6, the correctly rounded binary32 for that
///     [ 1  1  1  0 ]      decimal. It is deliberately NOT a power of two,
///     [ 0  1  1  1 ]      because 1/e is then exact and an unpivoted
///     [ 0  0  1  1 ]      reduction of this matrix would be exact too
///
/// Column 0 is `[e, 1, 0, 0]`, so the pivot is not on the diagonal and step 0
/// MUST interchange. Partial pivoting takes the 1 from row 1, so the first
/// multiplier is `e / 1 = 4e-6` and the largest multiplier of the whole
/// reduction is 1: the element growth factor is ONE.
///
/// WITHOUT the interchange the first multiplier would be `1 / e`, about
/// `2.5e5`, which binary32 cannot hold exactly. The trailing row would then be
/// scaled by `2.5e5` while its own entries are of magnitude 1, every later step
/// would work at that scale, and the growth factor of the reduction would be
/// `2.5e5`. Growth costs about `log2(2.5e5) = 18` of binary32's 24 significand
/// bits, so an unpivoted reduction has about six bits of result left, near
/// `2e-2`, while this one keeps all 24. Measured on this matrix with a
/// throwaway no-pivot elimination, the unpivoted solve of
/// `b = pivot4f * (1, -2, 3, -4)` is off by `7.3e-2` and the pivoted one is
/// exact. The test below asserts the pivoted half of that, and also asserts
/// the multiplier bound that is the whole content of partial pivoting; the
/// unpivoted path is not in this library and is not tested.
///
/// The matrix is not singular: its determinant is `-e = -4e-6` by cofactor
/// expansion along the first row, whose `1 * minor` term vanishes because that
/// minor has two equal columns. Its last pivot is exactly `e`, which clears
/// `matrix.inverse_pivot_threshold_f = 2^-21`, about `4.77e-7`, by a factor of
/// 8, so the matrix is inside the threshold and the test is about accuracy
/// rather than about acceptance.
///
/// The trace is derivable by hand: step 0 interchanges rows 0 and 1, step 1
/// interchanges rows 1 and 2, step 2 sees `[-1, 1]` in the column and keeps the
/// diagonal entry, and step 3 finishes. So `permutation` is `[1, 2, 0, 3]`, TWO
/// interchanges, the parity is `+1`, and the determinant is negative because a
/// pivot of U is negative.
const pivot4f = Mat4f.init(.{
    .{ 4e-6, 1, 0, 0 },
    .{ 1, 1, 1, 0 },
    .{ 0, 1, 1, 1 },
    .{ 0, 0, 1, 1 },
});

/// A singular 3 by 3: rows 0 and 1 are EQUAL, so they are linearly dependent
/// and the determinant is zero as a fact about the numbers rather than as a
/// rounding accident. In the elimination the two equal rows also receive the
/// same multiplier at every step, so the last pivot is the difference of two
/// identical quantities and is therefore exactly zero, not nearly zero.
const singular3f = Mat3f.init(.{ .{ 1, 2, 3 }, .{ 1, 2, 3 }, .{ 4, 5, 7 } });

/// An indefinite 2 by 2: symmetric, and the factorization finds it at the second
/// pivot. The first pivot is 1, the first entry below the diagonal is 2, and
/// the second pivot is `1 - 2 * 2 = -3`, which is not positive. Its
/// eigenvalues are 3 and -1, since `det(A - tI) = (1 - t)^2 - 4 =
/// t^2 - 2t - 3 = (t - 3)(t + 1)`.
const indefinite2f = Mat2f.init(.{ .{ 1, 2 }, .{ 2, 1 } });

/// Symmetric with a negative pivot, and a POSITIVE determinant, so a test that
/// decided positivity from the sign of a computed determinant would accept it.
/// `det(A - tI)` for `diag(-1, -1, 1)` is `(-1 - t)^2 (1 - t)`, so the
/// eigenvalues are -1, -1, and 1 and the determinant is their product, which
/// is 1. This is the case that forces the positivity test to be a pivot test
/// DURING the factorization.
const twoNegative3f = Mat3f.init(.{
    .{ -1, 0, 0 },
    .{ 0, -1, 0 },
    .{ 0, 0, 1 },
});

/// An asymmetric 2 by 2: `a[0][1] = 1` and `a[1][0] = 2`. Its determinant is
/// `4 * 3 - 2 * 1 = 1`, positive, so nothing about its determinant keeps it
/// out: it is the symmetry test that refuses it.
const asymmetric2f = Mat2f.init(.{ .{ 4, 1 }, .{ 2, 3 } });

// ---------------------------------------------------------------------------
// The constructor types
// ---------------------------------------------------------------------------

test "lu and cholesky return the contracted types" {
    try testing.expect(@TypeOf(linalg.lu(well2f)) == ?LuMat2f);
    try testing.expect(@TypeOf(linalg.lu(well3f)) == ?LuMat3f);
    try testing.expect(@TypeOf(linalg.lu(well4f)) == ?LuMat4f);
    try testing.expect(@TypeOf(linalg.cholesky(spd2f)) == ?CholeskyMat2f);
    try testing.expect(@TypeOf(linalg.cholesky(spd3f)) == ?CholeskyMat3f);
    try testing.expect(@TypeOf(linalg.cholesky(spd4f)) == ?CholeskyMat4f);

    const a2d = Mat2d.init(.{ .{ 4, 3 }, .{ 6, 3 } });
    const a3d = Mat3d.init(.{ .{ 2, 1, 1 }, .{ 1, 3, 2 }, .{ 1, 2, 3 } });
    const a4d = Mat4d.init(.{
        .{ 2, 1, 1, 0 },
        .{ 1, 3, 2, 1 },
        .{ 1, 2, 3, 1 },
        .{ 0, 1, 1, 2 },
    });
    try testing.expect(@TypeOf(linalg.lu(a2d)) == ?LuMat2d);
    try testing.expect(@TypeOf(linalg.lu(a3d)) == ?LuMat3d);
    try testing.expect(@TypeOf(linalg.lu(a4d)) == ?LuMat4d);

    const s2d = Mat2d.init(.{ .{ 4, 1 }, .{ 1, 3 } });
    const s3d = Mat3d.init(.{ .{ 4, 1, 1 }, .{ 1, 5, 2 }, .{ 1, 2, 6 } });
    const s4d = Mat4d.init(.{
        .{ 4, 1, 1, 0 },
        .{ 1, 5, 2, 1 },
        .{ 1, 2, 6, 2 },
        .{ 0, 1, 2, 7 },
    });
    try testing.expect(@TypeOf(linalg.cholesky(s2d)) == ?CholeskyMat2d);
    try testing.expect(@TypeOf(linalg.cholesky(s3d)) == ?CholeskyMat3d);
    try testing.expect(@TypeOf(linalg.cholesky(s4d)) == ?CholeskyMat4d);

    // The operations take and return the sizes the contract names.
    try testing.expect(@TypeOf(linalg.solveLu(linalg.lu(well3f).?, Vec3f.init(1, 1, 1))) == ?Vec3f);
    try testing.expect(@TypeOf(linalg.reconstructLu(linalg.lu(well3f).?)) == Mat3f);
    try testing.expect(@TypeOf(linalg.determinantLu(linalg.lu(well3f).?)) == f32);
    try testing.expect(
        @TypeOf(linalg.solveCholesky(linalg.cholesky(spd4f).?, Vec4f.init(1, 1, 1, 1))) == ?Vec4f,
    );
    try testing.expect(@TypeOf(linalg.reconstructCholesky(linalg.cholesky(spd4f).?)) == Mat4f);
    try testing.expect(@TypeOf(linalg.determinantCholesky(linalg.cholesky(spd4f).?)) == f32);
}

test "the factor types carry the contracted fields" {
    const f = linalg.lu(well3f).?;
    try testing.expectEqual(Mat3f, @TypeOf(f.lu));
    try testing.expectEqual([3]u8, @TypeOf(f.permutation));
    try testing.expectEqual(i1, @TypeOf(f.parity));

    const c = linalg.cholesky(spd4f).?;
    try testing.expectEqual(Mat4f, @TypeOf(c.lower));
}

// ---------------------------------------------------------------------------
// reconstruct: the central tests
// ---------------------------------------------------------------------------

test "reconstructLu rebuilds a well-conditioned matrix in every size and precision" {
    // binary32. The factors hold ratios such as 6/4, so the reconstruction is
    // not exact; the observed entry error is of order 1e-7 against entries of
    // magnitude 1, and 1e-6 is the tolerance that holds it.
    try expectClose(Mat2f, well2f, linalg.reconstructLu(linalg.lu(well2f).?), 1e-6);
    try expectClose(Mat3f, well3f, linalg.reconstructLu(linalg.lu(well3f).?), 1e-6);
    try expectClose(Mat4f, well4f, linalg.reconstructLu(linalg.lu(well4f).?), 1e-6);

    // binary64, the same three matrices in binary64. Its unit roundoff is
    // 2^-53, so the same reconstruction is within a few times 1e-15 and 1e-12
    // is the tolerance that holds it.
    const w2d = Mat2d.init(.{ .{ 4, 3 }, .{ 6, 3 } });
    const w3d = Mat3d.init(.{ .{ 2, 1, 1 }, .{ 1, 3, 2 }, .{ 1, 2, 3 } });
    const w4d = Mat4d.init(.{
        .{ 2, 1, 1, 0 },
        .{ 1, 3, 2, 1 },
        .{ 1, 2, 3, 1 },
        .{ 0, 1, 1, 2 },
    });
    try expectClose(Mat2d, w2d, linalg.reconstructLu(linalg.lu(w2d).?), 1e-12);
    try expectClose(Mat3d, w3d, linalg.reconstructLu(linalg.lu(w3d).?), 1e-12);
    try expectClose(Mat4d, w4d, linalg.reconstructLu(linalg.lu(w4d).?), 1e-12);
}

test "reconstructLu rebuilds an exactly factorizable matrix bit for bit" {
    // The first two rows of these are interchanged and everything else is an
    // identity row, so the factorization is exact and the reconstruction can be
    // compared entry for entry. This is what catches a permutation applied to
    // the wrong axis, which a tolerance wide enough to hide a transpose on a
    // matrix of ones would not.
    try expectExactly(Mat2f, swap2f, linalg.reconstructLu(linalg.lu(swap2f).?));
    try expectExactly(Mat3f, swap3f, linalg.reconstructLu(linalg.lu(swap3f).?));
    try expectExactly(Mat4f, swap4f, linalg.reconstructLu(linalg.lu(swap4f).?));

    const swap2d = Mat2d.init(.{ .{ 0, 1 }, .{ 1, 0 } });
    const swap3d = Mat3d.init(.{ .{ 0, 1, 0 }, .{ 1, 0, 0 }, .{ 0, 0, 1 } });
    const swap4d = Mat4d.init(.{
        .{ 0, 1, 0, 0 },
        .{ 1, 0, 0, 0 },
        .{ 0, 0, 1, 0 },
        .{ 0, 0, 0, 1 },
    });
    try expectExactly(Mat2d, swap2d, linalg.reconstructLu(linalg.lu(swap2d).?));
    try expectExactly(Mat3d, swap3d, linalg.reconstructLu(linalg.lu(swap3d).?));
    try expectExactly(Mat4d, swap4d, linalg.reconstructLu(linalg.lu(swap4d).?));
}

test "reconstructCholesky rebuilds a positive definite matrix in every size and precision" {
    try expectClose(Mat2f, spd2f, linalg.reconstructCholesky(linalg.cholesky(spd2f).?), 1e-6);
    try expectClose(Mat3f, spd3f, linalg.reconstructCholesky(linalg.cholesky(spd3f).?), 1e-6);
    try expectClose(Mat4f, spd4f, linalg.reconstructCholesky(linalg.cholesky(spd4f).?), 1e-6);

    const s2d = Mat2d.init(.{ .{ 4, 1 }, .{ 1, 3 } });
    const s3d = Mat3d.init(.{ .{ 4, 1, 1 }, .{ 1, 5, 2 }, .{ 1, 2, 6 } });
    const s4d = Mat4d.init(.{
        .{ 4, 1, 1, 0 },
        .{ 1, 5, 2, 1 },
        .{ 1, 2, 6, 2 },
        .{ 0, 1, 2, 7 },
    });
    try expectClose(Mat2d, s2d, linalg.reconstructCholesky(linalg.cholesky(s2d).?), 1e-12);
    try expectClose(Mat3d, s3d, linalg.reconstructCholesky(linalg.cholesky(s3d).?), 1e-12);
    try expectClose(Mat4d, s4d, linalg.reconstructCholesky(linalg.cholesky(s4d).?), 1e-12);
}

// ---------------------------------------------------------------------------
// Determinants, cross-checked against zmath.matrix
// ---------------------------------------------------------------------------

test "determinantLu agrees with the cofactor expansion of zmath.matrix" {
    // Two independent computations of the same quantity: this one is
    // `parity * product(diag(U))` and `zmath.matrix` expands along the first
    // row by cofactors with its own fixed operation order. They share no code,
    // so agreement is a cross-check and not a tautology.
    //
    // The observed deltas against the cofactor expansion are 0 for the 2 by 2,
    // 4.8e-7 for the 3 by 3, and 9.5e-7 for the 4 by 4; the tolerances below
    // sit one decimal digit above them. The 4 by 4 tolerance reads large next
    // to the others only because its determinant is 12.
    try testing.expectApproxEqAbs(
        matrix.determinantMat2f(well2f),
        linalg.determinantLu(linalg.lu(well2f).?),
        1e-5,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat3f(well3f),
        linalg.determinantLu(linalg.lu(well3f).?),
        1e-6,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat4f(well4f),
        linalg.determinantLu(linalg.lu(well4f).?),
        2e-6,
    );
    // The near-singular `pivot4f` is the one case where the two agree only to
    // a few parts in a thousand: its last pivot is `1 - fl(1 - e)`, so it
    // carries the rounding of `1 - e` a second time. The observed delta is
    // 6.5e-9 on a determinant of -4e-6.
    try testing.expectApproxEqAbs(
        matrix.determinantMat4f(pivot4f),
        linalg.determinantLu(linalg.lu(pivot4f).?),
        1e-8,
    );

    try testing.expectApproxEqAbs(
        matrix.determinantMat2f(swap2f),
        linalg.determinantLu(linalg.lu(swap2f).?),
        1e-6,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat3f(swap3f),
        linalg.determinantLu(linalg.lu(swap3f).?),
        1e-6,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat4f(swap4f),
        linalg.determinantLu(linalg.lu(swap4f).?),
        1e-6,
    );

    const w2d = Mat2d.init(.{ .{ 4, 3 }, .{ 6, 3 } });
    const w3d = Mat3d.init(.{ .{ 2, 1, 1 }, .{ 1, 3, 2 }, .{ 1, 2, 3 } });
    const w4d = Mat4d.init(.{
        .{ 2, 1, 1, 0 },
        .{ 1, 3, 2, 1 },
        .{ 1, 2, 3, 1 },
        .{ 0, 1, 1, 2 },
    });
    const swap2d = Mat2d.init(.{ .{ 0, 1 }, .{ 1, 0 } });
    const swap3d = Mat3d.init(.{ .{ 0, 1, 0 }, .{ 1, 0, 0 }, .{ 0, 0, 1 } });
    const swap4d = Mat4d.init(.{
        .{ 0, 1, 0, 0 },
        .{ 1, 0, 0, 0 },
        .{ 0, 0, 1, 0 },
        .{ 0, 0, 0, 1 },
    });
    try testing.expectApproxEqAbs(
        matrix.determinantMat2d(w2d),
        linalg.determinantLu(linalg.lu(w2d).?),
        1e-12,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat3d(w3d),
        linalg.determinantLu(linalg.lu(w3d).?),
        1e-12,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat4d(w4d),
        linalg.determinantLu(linalg.lu(w4d).?),
        1e-12,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat2d(swap2d),
        linalg.determinantLu(linalg.lu(swap2d).?),
        1e-12,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat3d(swap3d),
        linalg.determinantLu(linalg.lu(swap3d).?),
        1e-12,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat4d(swap4d),
        linalg.determinantLu(linalg.lu(swap4d).?),
        1e-12,
    );
}

test "determinantCholesky agrees with the cofactor expansion of zmath.matrix" {
    // The hand values are the full cofactor expansions of the matrices above:
    // 11 for the 2 by 2, 97 for the 3 by 3, and 608 for the 4 by 4, which are
    // the leading principal minors of `spd4f` taken all at once.
    try testing.expectApproxEqAbs(
        @as(f32, 11),
        linalg.determinantCholesky(linalg.cholesky(spd2f).?),
        2e-6,
    );
    try testing.expectApproxEqAbs(
        @as(f32, 97),
        linalg.determinantCholesky(linalg.cholesky(spd3f).?),
        2e-6,
    );
    try testing.expectApproxEqAbs(
        @as(f32, 608),
        linalg.determinantCholesky(linalg.cholesky(spd4f).?),
        2e-4,
    );

    try testing.expectApproxEqAbs(
        matrix.determinantMat2f(spd2f),
        linalg.determinantCholesky(linalg.cholesky(spd2f).?),
        2e-6,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat3f(spd3f),
        linalg.determinantCholesky(linalg.cholesky(spd3f).?),
        2e-6,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat4f(spd4f),
        linalg.determinantCholesky(linalg.cholesky(spd4f).?),
        2e-4,
    );

    const s2d = Mat2d.init(.{ .{ 4, 1 }, .{ 1, 3 } });
    const s3d = Mat3d.init(.{ .{ 4, 1, 1 }, .{ 1, 5, 2 }, .{ 1, 2, 6 } });
    const s4d = Mat4d.init(.{
        .{ 4, 1, 1, 0 },
        .{ 1, 5, 2, 1 },
        .{ 1, 2, 6, 2 },
        .{ 0, 1, 2, 7 },
    });
    try testing.expectApproxEqAbs(
        matrix.determinantMat2d(s2d),
        linalg.determinantCholesky(linalg.cholesky(s2d).?),
        1e-11,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat3d(s3d),
        linalg.determinantCholesky(linalg.cholesky(s3d).?),
        1e-12,
    );
    try testing.expectApproxEqAbs(
        matrix.determinantMat4d(s4d),
        linalg.determinantCholesky(linalg.cholesky(s4d).?),
        1e-12,
    );
}

// ---------------------------------------------------------------------------
// The permutation parity
// ---------------------------------------------------------------------------

test "an interchange is required, recorded, and applied to the sign" {
    // Column 0 of each matrix below is `[0, 1, ...]`, so the pivot is not on
    // the diagonal and step 0 MUST interchange. Exactly one interchange
    // follows, so the parity is -1. The stored diagonal product is +1, because
    // the factors are the identity, so the `-1` the caller sees can only have
    // come from the parity: an implementation that ignored it would return +1
    // and fail here.
    const f2 = linalg.lu(swap2f).?;
    const f3 = linalg.lu(swap3f).?;
    const f4 = linalg.lu(swap4f).?;

    try testing.expectEqual(@as(i1, -1), f2.parity);
    try testing.expectEqual(@as(i1, -1), f3.parity);
    try testing.expectEqual(@as(i1, -1), f4.parity);

    // The row each step promoted. Step 0 brought original row 1 into place 0,
    // and every later step found its pivot already in place.
    try testing.expectEqual(@as(u8, 1), f2.permutation[0]);
    try testing.expectEqual(@as(u8, 0), f2.permutation[1]);
    try testing.expectEqual(@as(u8, 1), f3.permutation[0]);
    try testing.expectEqual(@as(u8, 0), f3.permutation[1]);
    try testing.expectEqual(@as(u8, 2), f3.permutation[2]);
    try testing.expectEqual(@as(u8, 1), f4.permutation[0]);
    try testing.expectEqual(@as(u8, 0), f4.permutation[1]);
    try testing.expectEqual(@as(u8, 2), f4.permutation[2]);
    try testing.expectEqual(@as(u8, 3), f4.permutation[3]);

    try testing.expectEqual(@as(f32, 1), f2.lu.data[0][0]);
    try testing.expectEqual(@as(f32, 1), f2.lu.data[1][1]);
    try testing.expectEqual(@as(f32, -1), linalg.determinantLu(f2));
    try testing.expectEqual(@as(f32, -1), linalg.determinantLu(f3));
    try testing.expectEqual(@as(f32, -1), linalg.determinantLu(f4));

    // The binary64 siblings, where the factors are exact for the same reason.
    const swap2d = Mat2d.init(.{ .{ 0, 1 }, .{ 1, 0 } });
    const swap3d = Mat3d.init(.{ .{ 0, 1, 0 }, .{ 1, 0, 0 }, .{ 0, 0, 1 } });
    const swap4d = Mat4d.init(.{
        .{ 0, 1, 0, 0 },
        .{ 1, 0, 0, 0 },
        .{ 0, 0, 1, 0 },
        .{ 0, 0, 0, 1 },
    });
    try testing.expectEqual(@as(i1, -1), linalg.lu(swap2d).?.parity);
    try testing.expectEqual(@as(i1, -1), linalg.lu(swap3d).?.parity);
    try testing.expectEqual(@as(i1, -1), linalg.lu(swap4d).?.parity);
    try testing.expectEqual(@as(f64, -1), linalg.determinantLu(linalg.lu(swap2d).?));
    try testing.expectEqual(@as(f64, -1), linalg.determinantLu(linalg.lu(swap3d).?));
    try testing.expectEqual(@as(f64, -1), linalg.determinantLu(linalg.lu(swap4d).?));
}

test "a matrix that needs no interchange has the even parity" {
    // The largest entry of column 0 of `well3f` is already the pivot and so is
    // every later one, so no step interchanges. An `i1` cannot hold `+1`, so
    // the even case is stored as zero, which is the identity of the sign
    // multiply. Without this the `parity` field would only ever be seen at -1.
    try testing.expectEqual(@as(i1, 0), linalg.lu(well3f).?.parity);
    try testing.expectEqual(@as(i1, 0), linalg.lu(well4f).?.parity);
    try testing.expectEqual(@as(u8, 0), linalg.lu(well3f).?.permutation[0]);
    try testing.expectEqual(@as(u8, 3), linalg.lu(well4f).?.permutation[3]);

    // `well2f` interchanges exactly once, at step 0, and its determinant is
    // negative, so the parity and the sign of the stored diagonal agree.
    try testing.expectEqual(@as(i1, -1), linalg.lu(well2f).?.parity);
    try testing.expect(linalg.determinantLu(linalg.lu(well2f).?) < 0);

    // `pivot4f` interchanges twice, at step 0 and at step 1, and step 2 sees the
    // tie `[-1, 1]` and keeps the diagonal entry, so the parity is even while
    // the determinant is negative because a pivot of U is negative.
    const pivoted = linalg.lu(pivot4f).?;
    try testing.expectEqual(@as(i1, 0), pivoted.parity);
    try testing.expectEqual(@as(u8, 1), pivoted.permutation[0]);
    try testing.expectEqual(@as(u8, 2), pivoted.permutation[1]);
    try testing.expectEqual(@as(u8, 0), pivoted.permutation[2]);
    try testing.expectEqual(@as(u8, 3), pivoted.permutation[3]);
    try testing.expect(linalg.determinantLu(pivoted) < 0);
}

// ---------------------------------------------------------------------------
// Pivoting that matters
// ---------------------------------------------------------------------------

test "partial pivoting bounds every multiplier on a 4x4 whose pivot is below the diagonal" {
    const f = linalg.lu(pivot4f).?;
    // The interchange happened, and it happened at step 0, which is the only
    // place a no-pivot reduction could have differed: the pivot came from the
    // row below the diagonal.
    try testing.expectEqual(@as(u8, 1), f.permutation[0]);

    // Every multiplier is at or below ONE in magnitude, which is the whole
    // content of partial pivoting and what bounds the element growth of the
    // reduction.
    for (0..4) |column| {
        for (column + 1..4) |row| {
            try testing.expect(@abs(f.lu.data[column][row]) <= 1);
        }
    }
    // The multiplier that would break that bound is the one row 0 contributes
    // at step 0: `e / 1 = 4e-6` with the interchange, `1 / e`, about `2.5e5`,
    // without it. Step 1 then interchanges row 0 into position 2, so the
    // multiplier is stored in the third slot of the first column.
    try testing.expectApproxEqAbs(@as(f32, 4e-6), f.lu.data[0][2], 1e-12);

    try expectClose(Mat4f, pivot4f, linalg.reconstructLu(f), 1e-6);

    // The construct-and-solve round trip. This tests the SOLVE, not the
    // factorization, which the reconstruction above already covers. Six bits
    // of significand, which is what a growth factor of 2.5e5 leaves, would put
    // the answer near 2e-2; the tolerance below is two orders tighter.
    const x = Vec4f.init(1, -2, 3, -4);
    const solved = linalg.solveLu(f, pivot4f.mulVec(x)).?;
    try testing.expectApproxEqAbs(@as(f32, 1), solved.x, 1e-5);
    try testing.expectApproxEqAbs(@as(f32, -2), solved.y, 1e-5);
    try testing.expectApproxEqAbs(@as(f32, 3), solved.z, 1e-5);
    try testing.expectApproxEqAbs(@as(f32, -4), solved.w, 1e-5);
}

// ---------------------------------------------------------------------------
// Solve round trips
// ---------------------------------------------------------------------------

test "solveLu returns the constructed solution in every size and precision" {
    // The standard construct-and-solve round trip: pick a solution, form the
    // right-hand side by applying the matrix to it, solve, and compare. This
    // tests the SOLVE. The factorization itself is what `reconstructLu` above
    // covers, and a round trip on its own would also be satisfied by a
    // factorization and a solve that make the same mistake twice.
    const x2 = Vec2f.init(3, -1);
    const s2 = linalg.solveLu(linalg.lu(well2f).?, well2f.mulVec(x2)).?;
    try testing.expectApproxEqAbs(@as(f32, 3), s2.x, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, -1), s2.y, 1e-6);

    const x3 = Vec3f.init(1, 2, -3);
    const s3 = linalg.solveLu(linalg.lu(well3f).?, well3f.mulVec(x3)).?;
    try testing.expectApproxEqAbs(@as(f32, 1), s3.x, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 2), s3.y, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, -3), s3.z, 1e-6);

    const x4 = Vec4f.init(1, -2, 3, -4);
    const s4 = linalg.solveLu(linalg.lu(well4f).?, well4f.mulVec(x4)).?;
    try testing.expectApproxEqAbs(@as(f32, 1), s4.x, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, -2), s4.y, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 3), s4.z, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, -4), s4.w, 1e-6);

    // The same round trip in binary64, where the tolerance follows the binary64
    // unit roundoff rather than the binary32 one.
    const w2d = Mat2d.init(.{ .{ 4, 3 }, .{ 6, 3 } });
    const w3d = Mat3d.init(.{ .{ 2, 1, 1 }, .{ 1, 3, 2 }, .{ 1, 2, 3 } });
    const w4d = Mat4d.init(.{
        .{ 2, 1, 1, 0 },
        .{ 1, 3, 2, 1 },
        .{ 1, 2, 3, 1 },
        .{ 0, 1, 1, 2 },
    });
    const x2d = Vec2d.init(3, -1);
    const y2 = linalg.solveLu(linalg.lu(w2d).?, w2d.mulVec(x2d)).?;
    try testing.expectApproxEqAbs(@as(f64, 3), y2.x, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -1), y2.y, 1e-12);

    const x3d = Vec3d.init(1, 2, -3);
    const y3 = linalg.solveLu(linalg.lu(w3d).?, w3d.mulVec(x3d)).?;
    try testing.expectApproxEqAbs(@as(f64, 1), y3.x, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 2), y3.y, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -3), y3.z, 1e-12);

    const x4d = Vec4d.init(1, -2, 3, -4);
    const y4 = linalg.solveLu(linalg.lu(w4d).?, w4d.mulVec(x4d)).?;
    try testing.expectApproxEqAbs(@as(f64, 1), y4.x, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -2), y4.y, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 3), y4.z, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -4), y4.w, 1e-12);
}

test "solveCholesky returns the constructed solution in every size and precision" {
    // The same construct-and-solve round trip, against the Cholesky solve.
    const x2 = Vec2f.init(3, -1);
    const s2 = linalg.solveCholesky(linalg.cholesky(spd2f).?, spd2f.mulVec(x2)).?;
    try testing.expectApproxEqAbs(@as(f32, 3), s2.x, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, -1), s2.y, 1e-6);

    const x3 = Vec3f.init(1, 2, -3);
    const s3 = linalg.solveCholesky(linalg.cholesky(spd3f).?, spd3f.mulVec(x3)).?;
    try testing.expectApproxEqAbs(@as(f32, 1), s3.x, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 2), s3.y, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, -3), s3.z, 1e-6);

    const x4 = Vec4f.init(1, -2, 3, -4);
    const s4 = linalg.solveCholesky(linalg.cholesky(spd4f).?, spd4f.mulVec(x4)).?;
    try testing.expectApproxEqAbs(@as(f32, 1), s4.x, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, -2), s4.y, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 3), s4.z, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, -4), s4.w, 1e-6);

    const s2d = Mat2d.init(.{ .{ 4, 1 }, .{ 1, 3 } });
    const s3d = Mat3d.init(.{ .{ 4, 1, 1 }, .{ 1, 5, 2 }, .{ 1, 2, 6 } });
    const s4d = Mat4d.init(.{
        .{ 4, 1, 1, 0 },
        .{ 1, 5, 2, 1 },
        .{ 1, 2, 6, 2 },
        .{ 0, 1, 2, 7 },
    });
    const x2d = Vec2d.init(3, -1);
    const y2 = linalg.solveCholesky(linalg.cholesky(s2d).?, s2d.mulVec(x2d)).?;
    try testing.expectApproxEqAbs(@as(f64, 3), y2.x, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -1), y2.y, 1e-12);

    const x3d = Vec3d.init(1, 2, -3);
    const y3 = linalg.solveCholesky(linalg.cholesky(s3d).?, s3d.mulVec(x3d)).?;
    try testing.expectApproxEqAbs(@as(f64, 1), y3.x, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 2), y3.y, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -3), y3.z, 1e-12);

    const x4d = Vec4d.init(1, -2, 3, -4);
    const y4 = linalg.solveCholesky(linalg.cholesky(s4d).?, s4d.mulVec(x4d)).?;
    try testing.expectApproxEqAbs(@as(f64, 1), y4.x, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -2), y4.y, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 3), y4.z, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -4), y4.w, 1e-12);
}

test "solveMultiple solves every column against the same factorization" {
    // Three right-hand sides at once. The solutions are constructed the same way
    // as the round trips above, so this checks the batched solve against known
    // answers AND against the single solve, and a batched solve that mixed two
    // columns up would fail the first comparison.
    inline for (.{ 2, 3, 4 }) |size| {
        const A = switch (size) {
            2 => well2f,
            3 => well3f,
            else => well4f,
        };
        const V = switch (size) {
            2 => Vec2f,
            3 => Vec3f,
            else => Vec4f,
        };
        const X = Mat(size, 3, f32);
        var solutions: X = undefined;
        for (0..3) |column| {
            for (0..size) |row| {
                solutions.data[column][row] = @floatFromInt(
                    @as(i32, @intCast(row + 1)) * @as(i32, @intCast(column + 1)),
                );
            }
        }
        const rhs = matrix.mul(3, A, solutions);
        const solved = linalg.solveLuMultiple(linalg.lu(A).?, rhs).?;
        for (0..3) |column| {
            for (0..size) |row| {
                try testing.expectApproxEqAbs(
                    solutions.data[column][row],
                    solved.data[column][row],
                    1e-5,
                );
            }
            // And column by column, which is what the multiple form has to
            // agree with.
            const single = linalg.solveLu(linalg.lu(A).?, V.fromArray(rhs.data[column])).?;
            for (0..size) |row| {
                try testing.expectApproxEqAbs(
                    single.toArray()[row],
                    solved.data[column][row],
                    1e-5,
                );
            }
        }
    }

    // The Cholesky batched solve, against a symmetric positive definite matrix
    // whose solutions are constructed the same way.
    const A3 = Mat3f.init(.{ .{ 4, 1, 1 }, .{ 1, 5, 2 }, .{ 1, 2, 6 } });
    const X3 = Mat(3, 2, f32).init(.{
        .{ 1, -2 },
        .{ 2, 3 },
        .{ -3, 4 },
    });
    const rhs3 = matrix.mul(2, A3, X3);
    const solved3 = linalg.solveCholeskyMultiple(linalg.cholesky(A3).?, rhs3).?;
    for (0..2) |column| {
        for (0..3) |row| {
            try testing.expectApproxEqAbs(X3.data[column][row], solved3.data[column][row], 1e-5);
        }
        const single = linalg.solveCholesky(
            linalg.cholesky(A3).?,
            Vec3f.fromArray(rhs3.data[column]),
        ).?;
        for (0..3) |row| {
            try testing.expectApproxEqAbs(single.toArray()[row], solved3.data[column][row], 1e-5);
        }
    }
}

// ---------------------------------------------------------------------------
// Refusals
// ---------------------------------------------------------------------------

test "lu refuses a singular matrix" {
    // Rows 0 and 1 of `singular3f` are equal, so they are linearly dependent
    // and the determinant is zero as a fact about the numbers rather than as a
    // rounding accident. The two rows also receive the same multiplier at every
    // step, so the last pivot is the difference of two identical quantities and
    // is exactly zero.
    //
    // It is not the SOLVE's job to make this possible, and the test asserts
    // that the CONSTRUCTOR refuses: a factorization that existed for a singular
    // matrix would have to invent a pivot somewhere.
    try testing.expect(linalg.lu(singular3f) == null);
    try testing.expect(linalg.lu(well3f) != null);
}

test "lu refuses a pivot at or below the relative threshold" {
    // The threshold is `matrix.inverse_pivot_threshold_f`, which is `2^-21`,
    // times the largest entry magnitude of the INPUT. In both matrices below
    // that largest entry is 1, so the threshold is `2^-21` itself and the
    // second pivot is `t`: a `t` above `2^-21` is accepted and a `t` below it is
    // refused. That makes the test a MAGNITUDE test rather than an equality
    // test, since the two cases differ by a factor of sixteen and not by a
    // rounding difference, and the accepted matrix has a determinant of `t` just
    // as nonzero as the refused one has of `t`.
    try testing.expect(
        linalg.lu(Mat2f.init(.{ .{ 1, 0 }, .{ 0, 8.0 * @as(f32, 0x1p-21) } })) != null,
    );
    try testing.expect(
        linalg.lu(Mat2f.init(.{ .{ 1, 0 }, .{ 0, 0.5 * @as(f32, 0x1p-21) } })) == null,
    );
    // An identity is not singular, and a zero matrix is: the threshold cannot
    // make an invertible matrix look singular, and the equality test it
    // replaces would have missed the small-but-usable pivot above.
    try testing.expect(linalg.lu(Mat2f.identity()) != null);
    try testing.expect(linalg.lu(Mat2f.fill(0)) == null);

    // The binary64 threshold is `2^-50`: the same derivation at binary64's unit
    // roundoff of `2^-53`.
    try testing.expect(
        linalg.lu(Mat2d.init(.{ .{ 1, 0 }, .{ 0, 8.0 * @as(f64, 0x1p-50) } })) != null,
    );
    try testing.expect(
        linalg.lu(Mat2d.init(.{ .{ 1, 0 }, .{ 0, 0.5 * @as(f64, 0x1p-50) } })) == null,
    );
}

test "cholesky refuses an indefinite matrix" {
    // `indefinite2f` has eigenvalues 3 and -1, and the factorization finds it
    // at the second pivot: the first pivot is 1, the entry below the diagonal
    // is 2, and the second pivot is `1 - 2 * 2 = -3`.
    try testing.expect(linalg.cholesky(indefinite2f) == null);

    // `twoNegative3f` has a POSITIVE determinant and is still refused. Its
    // eigenvalues are -1, -1, and 1 and its determinant is their product, which
    // is 1. A test that decided positivity from the sign of a computed
    // determinant would accept this matrix; a pivot test cannot.
    try testing.expect(matrix.determinantMat3f(twoNegative3f) > 0);
    try testing.expect(linalg.cholesky(twoNegative3f) == null);
}

test "cholesky refuses a matrix that is not symmetric" {
    // The asymmetry policy is refusal, not symmetrization.
    // `asymmetric2f` has `a[0][1] = 1` and `a[1][0] = 2`, and its determinant is
    // `4 * 3 - 2 * 1 = 1`, positive, so nothing about its determinant keeps it
    // out: it is the symmetry test that refuses it.
    try testing.expect(matrix.determinantMat2f(asymmetric2f) > 0);
    try testing.expect(linalg.cholesky(asymmetric2f) == null);
    try testing.expect(linalg.cholesky(spd2f) != null);

    // A matrix symmetric to within a rounding difference is refused as well,
    // because this module has no tolerance and one invented inside a
    // factorization would be the implicit epsilon the module doc rules out. A
    // caller who wants it symmetrized symmetrizes it at the call site, where the
    // decision is visible, and the result factorizes.
    const almost = Mat2f.init(.{ .{ 4, 1 }, .{ 1 + @as(f32, 0x1p-20), 3 } });
    try testing.expect(linalg.cholesky(almost) == null);
    const symmetrized = almost.add(almost.transpose()).scale(0.5);
    try testing.expect(linalg.cholesky(symmetrized) != null);

    // A matrix whose off-diagonal zeros differ only in sign compares EQUAL as
    // numbers, so it is accepted: the documented policy is a NUMERIC symmetry
    // test, not a bitwise one, and `+0 == -0` is part of that.
    const signed_zero = Mat2f.init(.{ .{ 4, -0.0 }, .{ 0.0, 3 } });
    try testing.expect(linalg.cholesky(signed_zero) != null);
    try expectClose(
        Mat2f,
        signed_zero,
        linalg.reconstructCholesky(linalg.cholesky(signed_zero).?),
        1e-6,
    );
}

/// A finite, symmetric, strictly diagonally dominant matrix of type `M`, so it
/// factorizes by both methods before anything is poisoned.
fn seedOf(comptime M: type) M {
    const n = comptime M.rows;
    var out: M = undefined;
    for (0..n) |row| {
        for (0..n) |column| {
            out.data[column][row] = if (row == column)
                @as(M.Element, @floatFromInt(n + 2))
            else
                @as(M.Element, @floatFromInt(@min(row, column) + 1));
        }
    }
    return out;
}

test "the seed matrices are accepted before anything is poisoned" {
    // The sweep below proves nothing if its seed was already refused, so the
    // seed is checked on its own here, in every size and both precisions.
    inline for (.{ Mat2f, Mat3f, Mat4f, Mat2d, Mat3d, Mat4d }) |M| {
        const base = seedOf(M);
        try testing.expect(linalg.lu(base) != null);
        try testing.expect(linalg.cholesky(base) != null);
        try expectClose(M, base, linalg.reconstructLu(linalg.lu(base).?), 1e-5);
        try expectClose(M, base, linalg.reconstructCholesky(linalg.cholesky(base).?), 1e-5);
    }
}

/// The `lu` of `base` with one entry replaced by a NaN or by an infinity is
/// `null`, at EVERY position.
fn expectLuRejects(comptime M: type, base: M, comptime want_nan: bool) !void {
    const n = comptime M.rows;
    inline for (0..n) |row| {
        inline for (0..n) |column| {
            var poisoned = base;
            poisoned.data[column][row] = specialOf(M, want_nan);
            if (linalg.lu(poisoned) != null) {
                std.debug.print("lu accepted a non-finite entry at ({d}, {d})\n", .{ row, column });
                return error.NonFiniteEntryAccepted;
            }
        }
    }
}

/// The `cholesky` of `base` with one entry replaced by a NaN or by an infinity
/// is `null`, at EVERY position.
fn expectCholeskyRejects(comptime M: type, base: M, comptime want_nan: bool) !void {
    const n = comptime M.rows;
    inline for (0..n) |row| {
        inline for (0..n) |column| {
            var poisoned = base;
            poisoned.data[column][row] = specialOf(M, want_nan);
            if (linalg.cholesky(poisoned) != null) {
                std.debug.print(
                    "cholesky accepted a non-finite entry at ({d}, {d})\n",
                    .{ row, column },
                );
                return error.NonFiniteEntryAccepted;
            }
        }
    }
}

test "both constructors refuse a non-finite entry in every size and precision" {
    // Every position in turn, for a NaN and for an infinity separately, in every
    // size and both precisions. Checking every position is what proves the
    // finiteness test covers the whole matrix and not a column or a diagonal of
    // it, and running it for Cholesky before the symmetry test is what proves
    // the order of those two checks.
    inline for (.{ Mat2f, Mat3f, Mat4f, Mat2d, Mat3d, Mat4d }) |M| {
        const base = seedOf(M);
        try expectLuRejects(M, base, false);
        try expectLuRejects(M, base, true);
        try expectCholeskyRejects(M, base, false);
        try expectCholeskyRejects(M, base, true);
    }
}

// ---------------------------------------------------------------------------
// The stored layout
// ---------------------------------------------------------------------------

test "the packed factors are L below the diagonal and U on and above it" {
    // `well2f` has a first column of `[4, 6]`, so step 0 interchanges rows 0
    // and 1. After that the matrix is `[6, 3; 4, 3]`, the pivot is 6, the
    // multiplier is `4 / 6`, and the trailing entry is `3 - (4 / 6) * 3 = 1`.
    // So U is `[6, 3; 0, 1]`, L is `[1, 0; 4/6, 1]`, and the determinant is
    // `parity * 6 * 1 = -6`, which is `4 * 3 - 3 * 6`. The exactness here is
    // what makes the layout checkable: a transposed packing would put the 3
    // where the multiplier is.
    const f = linalg.lu(well2f).?;
    try testing.expectApproxEqAbs(@as(f32, 6), f.lu.data[0][0], 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 3), f.lu.data[1][0], 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 2.0 / 3.0), f.lu.data[0][1], 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 1), f.lu.data[1][1], 1e-6);

    // The Cholesky factor of `spd2f`: `l[0][0] = sqrt(4) = 2`, `l[1][0] = 1 / 2`,
    // and `l[1][1] = sqrt(3 - 1/4)`, whose square is `3 - 1/4 = 11/4`. So
    // `L * transpose(L)` is `[[4, 1], [1, 1/4 + 11/4]] = [[4, 1], [1, 3]]`.
    const c = linalg.cholesky(spd2f).?;
    try testing.expectApproxEqAbs(@as(f32, 2), c.lower.data[0][0], 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 0.5), c.lower.data[0][1], 1e-6);
    try testing.expectApproxEqAbs(
        @as(f32, 11) / 4.0,
        c.lower.data[1][1] * c.lower.data[1][1],
        1e-6,
    );

    // The strict upper triangle of a Cholesky factor is ZERO, not a packed
    // entry that happens to be small. A caller that reads only the lower
    // triangle can rely on that, and an implementation that silently kept the
    // input's upper triangle would fail here.
    const c4 = linalg.cholesky(spd4f).?;
    for (0..4) |row| {
        for (row + 1..4) |column| {
            try testing.expectEqual(@as(f32, 0), c4.lower.data[column][row]);
        }
    }
}
