//! Fixed-size matrix factorizations: the module `API.md:983-994` specifies.
//!
//! Two factorizations, each in the sizes 2, 3, and 4 and each in `f32` and
//! `f64`: an LU with partial pivoting, and a Cholesky of a symmetric positive
//! definite matrix. Each constructor returns an optional, and each of
//! `reconstruct`, `determinant`, and `solve` takes the factorization the
//! constructor produced. `API.md:994` says "a factorization returns `null` for
//! a singular, non-positive-definite, or non-finite input, as applicable", so
//! the refusal belongs to the constructor and the three operations cannot
//! fail on a matrix that was accepted, except in the one case the contract
//! also makes an optional: a `solve` returns `null` if the stored factors do
//! not have the usable pivots, which a hand-built factorization type can
//! violate and a constructor cannot.
//!
//! ## Specialized operations are module functions
//!
//! `solveLu`, `reconstructLu`, `determinantLu`, and their Cholesky siblings
//! are on the module, not methods on the factorization type, and they are the
//! same spellings `zmath.matrix` uses for `determinantMat3f` and
//! `inverseMat4d` (`API.md:971` names that rule for the surface as a whole).
//! The constructors take the size and the precision from their ARGUMENT rather
//! than from their name, because the contract names the TYPES and not the
//! functions. `lu(a)` for a `Mat3f` has type `?LuMat3f` and `cholesky(a)` for
//! a `Mat4f` has type `?CholeskyMat4f`: one name reaches all six of each, and
//! the return type is computed from the argument's type.
//!
//! ## The LU layout
//!
//! DOOLITTLE with partial pivoting. `lu` packs both factors into one matrix:
//! the strictly lower triangle is L, whose diagonal is UNIT, and the diagonal
//! and the upper triangle are U. That is the layout LAPACK `getrf` leaves
//! behind. Doolittle and not Crout, for one reason that touches every
//! operation below: Crout puts the unit diagonal on U, so the forward
//! substitution divides by L's diagonal, and the pivots live in a part of the
//! storage that is not the diagonal. Doolittle puts the unit diagonal on L, so
//! forward substitution is a subtraction with no division, and U's diagonal is
//! the pivot sequence itself, which makes `determinantLu` a product of stored
//! entries rather than a reassembly.
//!
//! `permutation[k]` is the index, in the ORIGINAL matrix, of the row that step
//! `k` promoted to be the k-th pivot row. It is the interchange as applied: at
//! step `k` the largest magnitude in column `k` among the rows that have not
//! yet supplied a pivot is moved into position `k`, and that row is the one
//! named. Each row becomes exactly one pivot row and position `k` is never
//! touched again, so the array is always a permutation of `0..n-1`.
//!
//! `parity` is the sign of that permutation: `-1` for an odd number of
//! interchanges and `0` for an even one. An `i1` cannot hold `+1`, and zero is
//! the identity of the sign multiply, so zero is the neutral value the type
//! allows. It exists because a factorization is only a factorization of
//! `permuted(A)`: the stored factors obey `P * A = L * U`, so
//! `det(P) * det(A) = det(L) * det(U)`, which gives
//! `det(A) = parity * product(diag(U))` and `reconstructLu` puts row
//! `permutation[k]` back where it came from. Both the parity and the row order
//! therefore have to agree with the stored factors, and
//! `reconstructLu == A` is what proves they do.
//!
//! ## The LU singularity threshold
//!
//! A pivot magnitude threshold, never an equality test, and the constants are
//! `zmath.matrix`'s own: `matrix.inverse_pivot_threshold_f` and
//! `matrix.inverse_pivot_threshold_d`. They are the right constants for this
//! elimination. Each is the round-to-nearest unit roundoff of the format
//! (`2^-24` for binary32, `2^-53` for binary64) times the `2^(n-1)`
//! element-growth bound of a reduction, which for the largest size here,
//! `n = 4`, is 8: `2^-21` and `2^-50`. The test is relative to the largest
//! entry magnitude of the WHOLE INPUT, read before any arithmetic, for the
//! reason `zmath.matrix` gives: a column-relative scale can itself collapse to
//! zero by the time that column's turn comes. A second spelling of one
//! threshold in two modules is how two modules drift apart, so this module
//! names the constants `zmath.matrix` publishes rather than declaring its own.
//! There is no epsilon argument, and no implicit epsilon anywhere.
//!
//! ## The Cholesky asymmetry policy
//!
//! A matrix that is not SYMMETRIC is refused with `null`. It is not silently
//! symmetrized and it is not read from one triangle only. A matrix that is
//! symmetric to within a rounding difference is refused as well, because there
//! is no tolerance to apply here and a tolerance invented inside a
//! factorization would be the implicit epsilon this module does not have: a
//! caller who wants a symmetrized input symmetrizes it at the call site,
//! `(value + value.transpose()).scale(0.5)`, where the decision is visible.
//! The symmetry test is a numeric equality, so `+0` and `-0` compare equal and
//! a matrix whose off-diagonal zeros differ only in sign is accepted.
//!
//! ## The Cholesky pivot test
//!
//! The positive-definiteness test is a pivot test DURING the factorization,
//! never a sign on a determinant computed afterwards. A matrix can have a
//! positive determinant and still fail Cholesky: `diag(-1, -1, 1)` has
//! determinant one and two negative eigenvalues. The pivot test is the same
//! relative threshold as the LU one, which makes it strictly stronger than
//! positivity, and that strength is deliberate for the reason
//! `zmath.matrix` states: a pivot at or below the threshold has no significant
//! bits left, so dividing the column by it would amplify rounding noise
//! instead of recovering a factor. A positive-definite matrix whose smallest
//! pivot falls at or below the threshold is reported as not positive definite
//! to working precision, which is the same honest answer an inverse gives for
//! the analogous matrix.
//!
//! ## Arithmetic
//!
//! Every arithmetic body is `@setFloatMode(.strict)`, so the compiler may
//! neither reassociate a sum nor contract a multiply and an add into an FMA,
//! and every sum of products starts from the index-zero product rather than
//! from an explicit zero. That is `zmath.matrix`'s discipline, kept here so
//! the two modules agree, and the reason this module composes
//! `matrix.mul` for its products instead of writing its own contraction. An
//! entry of `f32` or `f64` is IEEE arithmetic with the NaN rules of
//! `spec/core-math-2026-08-22.md:769-801`; both constructors reject a
//! non-finite input before any arithmetic runs, because a step on a non-finite
//! entry cannot be undone by a check afterwards. zmath never reads, validates,
//! sets, or restores `fenv`; the caller must keep round-to-nearest,
//! ties-to-even active.
//!
//! ## Storage and sizes
//!
//! A matrix is `zmath.matrix`'s `Mat`, stored COLUMN-MAJOR as `data[column]
//! [row]`, and this module never reinterprets that: `data[column][row]` is
//! row `row`, column `column` everywhere below. No operation here allocates.
//! A factorization is a fixed-size struct of fixed-size arrays, and a
//! multiple solve holds its scratch in two comptime-bounded arrays on the
//! stack.

const std = @import("std");
const matrix = @import("matrix.zig");
const scalar = @import("scalar.zig");

const Mat = matrix.Mat;
const Mat2f = matrix.Mat2f;
const Mat3f = matrix.Mat3f;
const Mat4f = matrix.Mat4f;
const Mat2d = matrix.Mat2d;
const Mat3d = matrix.Mat3d;
const Mat4d = matrix.Mat4d;

// ---------------------------------------------------------------------------
// Factorization types
//
// The names and the field layout are the ones `API.md:988-991` writes. Every
// type is spelled out rather than produced by a generic factory, so the source
// and the contract can be read side by side.
// ---------------------------------------------------------------------------

/// An LU factorization of a 2 by 2 binary32 matrix, with partial pivoting.
pub const LuMat2f = struct {
    /// L below the diagonal and U on and above it. L's diagonal is UNIT and is
    /// not stored: the strictly lower triangle is L, and the diagonal and the
    /// upper triangle are U.
    lu: Mat2f,
    /// The row of the ORIGINAL matrix that each step promoted to be the pivot
    /// row of its own step. A permutation of `0..1`.
    permutation: [2]u8,
    /// The sign of that permutation: `-1` for an odd number of interchanges
    /// and `0` for an even one. An `i1` cannot hold `+1`, and zero is the
    /// identity of the sign multiply, so zero is the neutral value here.
    parity: i1,
};

/// An LU factorization of a 3 by 3 binary32 matrix, with partial pivoting.
pub const LuMat3f = struct {
    /// L below the diagonal and U on and above it. L's diagonal is UNIT and is
    /// not stored.
    lu: Mat3f,
    /// The row of the ORIGINAL matrix that each step promoted to be the pivot
    /// row of its own step. A permutation of `0..2`.
    permutation: [3]u8,
    /// The sign of that permutation: `-1` for an odd number of interchanges
    /// and `0` for an even one.
    parity: i1,
};

/// An LU factorization of a 4 by 4 binary32 matrix, with partial pivoting.
pub const LuMat4f = struct {
    /// L below the diagonal and U on and above it. L's diagonal is UNIT and is
    /// not stored.
    lu: Mat4f,
    /// The row of the ORIGINAL matrix that each step promoted to be the pivot
    /// row of its own step. A permutation of `0..3`.
    permutation: [4]u8,
    /// The sign of that permutation: `-1` for an odd number of interchanges
    /// and `0` for an even one.
    parity: i1,
};

/// An LU factorization of a 2 by 2 binary64 matrix, with partial pivoting.
pub const LuMat2d = struct {
    /// L below the diagonal and U on and above it. L's diagonal is UNIT and is
    /// not stored.
    lu: Mat2d,
    /// The row of the ORIGINAL matrix that each step promoted to be the pivot
    /// row of its own step. A permutation of `0..1`.
    permutation: [2]u8,
    /// The sign of that permutation: `-1` for an odd number of interchanges
    /// and `0` for an even one.
    parity: i1,
};

/// An LU factorization of a 3 by 3 binary64 matrix, with partial pivoting.
pub const LuMat3d = struct {
    /// L below the diagonal and U on and above it. L's diagonal is UNIT and is
    /// not stored.
    lu: Mat3d,
    /// The row of the ORIGINAL matrix that each step promoted to be the pivot
    /// row of its own step. A permutation of `0..2`.
    permutation: [3]u8,
    /// The sign of that permutation: `-1` for an odd number of interchanges
    /// and `0` for an even one.
    parity: i1,
};

/// An LU factorization of a 4 by 4 binary64 matrix, with partial pivoting.
pub const LuMat4d = struct {
    /// L below the diagonal and U on and above it. L's diagonal is UNIT and is
    /// not stored.
    lu: Mat4d,
    /// The row of the ORIGINAL matrix that each step promoted to be the pivot
    /// row of its own step. A permutation of `0..3`.
    permutation: [4]u8,
    /// The sign of that permutation: `-1` for an odd number of interchanges
    /// and `0` for an even one.
    parity: i1,
};

/// A Cholesky factorization of a 2 by 2 symmetric positive definite binary32
/// matrix.
pub const CholeskyMat2f = struct {
    /// The lower triangular factor L, so the original matrix is `L *
    /// transpose(L)`. The strict upper triangle is ZERO, which is what lets a
    /// caller see that this is one triangle and not a packed matrix.
    lower: Mat2f,
};

/// A Cholesky factorization of a 3 by 3 symmetric positive definite binary32
/// matrix.
pub const CholeskyMat3f = struct {
    /// The lower triangular factor L, so the original matrix is `L *
    /// transpose(L)`. The strict upper triangle is ZERO.
    lower: Mat3f,
};

/// A Cholesky factorization of a 4 by 4 symmetric positive definite binary32
/// matrix.
pub const CholeskyMat4f = struct {
    /// The lower triangular factor L, so the original matrix is `L *
    /// transpose(L)`. The strict upper triangle is ZERO.
    lower: Mat4f,
};

/// A Cholesky factorization of a 2 by 2 symmetric positive definite binary64
/// matrix.
pub const CholeskyMat2d = struct {
    /// The lower triangular factor L, so the original matrix is `L *
    /// transpose(L)`. The strict upper triangle is ZERO.
    lower: Mat2d,
};

/// A Cholesky factorization of a 3 by 3 symmetric positive definite binary64
/// matrix.
pub const CholeskyMat3d = struct {
    /// The lower triangular factor L, so the original matrix is `L *
    /// transpose(L)`. The strict upper triangle is ZERO.
    lower: Mat3d,
};

/// A Cholesky factorization of a 4 by 4 symmetric positive definite binary64
/// matrix.
pub const CholeskyMat4d = struct {
    /// The lower triangular factor L, so the original matrix is `L *
    /// transpose(L)`. The strict upper triangle is ZERO.
    lower: Mat4d,
};

// ---------------------------------------------------------------------------
// Type mapping
//
// The sizes are named, not generic, so each public operation takes `anytype`
// and computes its return type from the argument. That is the shape
// `zmath.matrix.mul` uses for the same reason: the contract calls these
// operations by matrix type, and Zig cannot infer a trailing type parameter
// from an argument. The switches below are also the module's guard: a matrix
// that is not one of these six, or a value that is not one of these twelve
// factorization types, is a compile error naming what is wrong.
// ---------------------------------------------------------------------------

/// The LU type for a matrix type.
fn LuOf(comptime M: type) type {
    return switch (M) {
        Mat2f => LuMat2f,
        Mat3f => LuMat3f,
        Mat4f => LuMat4f,
        Mat2d => LuMat2d,
        Mat3d => LuMat3d,
        Mat4d => LuMat4d,
        else => @compileError(
            "zmath.linalg: '" ++ @typeName(M) ++ "' is not a matrix this module factorizes. " ++
                "The sizes are Mat2f, Mat3f, Mat4f, Mat2d, Mat3d, and Mat4d.",
        ),
    };
}

/// The Cholesky type for a matrix type.
fn CholeskyOf(comptime M: type) type {
    return switch (M) {
        Mat2f => CholeskyMat2f,
        Mat3f => CholeskyMat3f,
        Mat4f => CholeskyMat4f,
        Mat2d => CholeskyMat2d,
        Mat3d => CholeskyMat3d,
        Mat4d => CholeskyMat4d,
        else => @compileError(
            "zmath.linalg: '" ++ @typeName(M) ++ "' is not a matrix this module factorizes. " ++
                "The sizes are Mat2f, Mat3f, Mat4f, Mat2d, Mat3d, and Mat4d.",
        ),
    };
}

/// The matrix a factorization type factorizes.
fn MatrixOf(comptime F: type) type {
    return switch (F) {
        LuMat2f, CholeskyMat2f => Mat2f,
        LuMat3f, CholeskyMat3f => Mat3f,
        LuMat4f, CholeskyMat4f => Mat4f,
        LuMat2d, CholeskyMat2d => Mat2d,
        LuMat3d, CholeskyMat3d => Mat3d,
        LuMat4d, CholeskyMat4d => Mat4d,
        else => @compileError(
            "zmath.linalg: '" ++ @typeName(F) ++ "' is not a factorization type of this module. " ++
                "Build one with zmath.linalg.lu or zmath.linalg.cholesky.",
        ),
    };
}

/// The matrix type a set of right-hand sides is stored in, which is the shape
/// the matching solve returns.
fn RightHandSidesOf(comptime F: type, comptime R: type) type {
    const M = MatrixOf(F);
    if (@typeInfo(R) != .@"struct") {
        @compileError("zmath.linalg: the right-hand sides must be a Mat, not '" ++ @typeName(R) ++ "'");
    }
    if (@field(R, "rows") != M.rows or @field(R, "columns") < 1) {
        @compileError(
            "zmath.linalg: the right-hand sides must have " ++
                std.fmt.comptimePrint("{d}", .{M.rows}) ++
                " rows and at least one column",
        );
    }
    if (@field(R, "Element") != M.Element) {
        @compileError(
            "zmath.linalg: the right-hand sides have entry type '" ++
                @typeName(@field(R, "Element")) ++ "' but the factorization has '" ++
                @typeName(M.Element) ++ "'.",
        );
    }
    return Mat(M.rows, R.columns, M.Element);
}

// ---------------------------------------------------------------------------
// Shared entry arithmetic
// ---------------------------------------------------------------------------

/// The smallest pivot, relative to the largest entry magnitude of the input,
/// that a factorization will divide by: `zmath.matrix`'s constant, reused for
/// the reason the module comment gives.
fn pivotThresholdOf(comptime T: type) T {
    return switch (T) {
        f32 => matrix.inverse_pivot_threshold_f,
        f64 => matrix.inverse_pivot_threshold_d,
        else => @compileError("zmath.linalg: only a real binary32 or binary64 entry has a pivot magnitude"),
    };
}

/// The accepted square root kernel, so a Cholesky diagonal comes from the same
/// `zmath.scalar` rounding as every other square root in the library.
fn sqrtOf(comptime T: type, value: T) T {
    return switch (T) {
        f32 => scalar.sqrtf(value),
        f64 => scalar.sqrtd(value),
        else => @compileError("zmath.linalg: only a real binary32 or binary64 entry has a square root"),
    };
}

/// The largest entry magnitude of `value`, or `null` when any entry is not
/// finite.
///
/// Read over every entry BEFORE any arithmetic, because a step on a non-finite
/// entry cannot be undone by a check afterwards. `zmath.matrix` checks in the
/// same order, so the two modules accept the same inputs.
fn finiteScaleOf(comptime M: type, value: M) ?M.Element {
    @setFloatMode(.strict);
    const n = comptime M.rows;
    var scale: M.Element = 0;
    inline for (0..n) |row| {
        inline for (0..n) |column| {
            const entry = value.data[column][row];
            if (!std.math.isFinite(entry)) return null;
            const magnitude = @abs(entry);
            if (magnitude > scale) scale = magnitude;
        }
    }
    return scale;
}

/// `null` when any entry of `value` is not finite, and nothing otherwise.
///
/// This is the post-condition check on the FACTORS rather than on the input,
/// so it runs last: a matrix that passed every pivot test can still drive an
/// entry of a factor past the format's range, and reporting that as "no
/// factorization" is the honest answer while returning an infinity is not.
fn allFinite(comptime M: type, value: M) ?void {
    @setFloatMode(.strict);
    const n = comptime M.rows;
    inline for (0..n) |row| {
        inline for (0..n) |column| {
            if (!std.math.isFinite(value.data[column][row])) return null;
        }
    }
    return {};
}

/// The row count of a square matrix type, for the scratch array types below.
fn n_of(comptime M: type) usize {
    return comptime M.rows;
}

// ---------------------------------------------------------------------------
// LU
// ---------------------------------------------------------------------------

/// The LU factorization of a square real matrix by Doolittle elimination with
/// PARTIAL PIVOTING, or `null`.
///
/// `null` means one of three things: an entry was not finite, a pivot fell at
/// or below `matrix.inverse_pivot_threshold_f` or
/// `matrix.inverse_pivot_threshold_d` times the largest entry magnitude of the
/// input, or the factors themselves ran past the format's range. Pivoting is
/// not optional here. It is what keeps the multipliers at or below one, which
/// is what bounds the element growth at `2^(n-1)` and therefore what makes the
/// threshold above meaningful; it is also what puts a row interchange in the
/// result, and so what gives `permutation` and `parity` a purpose at all.
fn luFactor(comptime M: type, value: M) ?LuOf(M) {
    @setFloatMode(.strict);
    const n = comptime M.rows;
    const T = M.Element;
    const threshold = pivotThresholdOf(T);

    const scale = finiteScaleOf(M, value) orelse return null;

    var work = value;
    // The row of the ORIGINAL matrix at each position, carried alongside the
    // rows themselves. Swapping the two together is what keeps `permutation`
    // and the stored factors describing the same matrix.
    var permutation: [n]u8 = undefined;
    inline for (0..n) |row| permutation[row] = @intCast(row);
    // An `i1` cannot hold `+1`, so the sign of the permutation is stored the
    // only way the type allows: `-1` for an odd number of interchanges and `0`
    // for an even one. Zero is the identity of the sign multiply, so it is the
    // right neutral value, and the flip below is a bitwise NOT rather than a
    // negation because negating a one-bit integer cannot represent the result.
    var parity: i1 = 0;

    inline for (0..n) |k| {
        // Partial pivoting: the largest magnitude in column `k` among the rows
        // that have not yet supplied a pivot row. A tie keeps the current row,
        // so a pivot that is already in place costs no interchange and no
        // parity flip.
        var pivot_row: usize = k;
        var pivot_magnitude: T = @abs(work.data[k][k]);
        for (k + 1..n) |row| {
            const magnitude = @abs(work.data[k][row]);
            if (magnitude > pivot_magnitude) {
                pivot_magnitude = magnitude;
                pivot_row = row;
            }
        }
        if (!(pivot_magnitude > threshold * scale)) return null;

        if (pivot_row != k) {
            // Swap whole rows: the stored L entries below the diagonal move
            // with them, which is the whole point of packing the two factors
            // into one matrix.
            for (0..n) |column| {
                const held = work.data[column][k];
                work.data[column][k] = work.data[column][pivot_row];
                work.data[column][pivot_row] = held;
            }
            const held_row = permutation[k];
            permutation[k] = permutation[pivot_row];
            permutation[pivot_row] = held_row;
            // One transposition flips the sign of a permutation.
            parity = ~parity;
        }

        // Doolittle: the unit diagonal of L is implicit, so the multipliers go
        // straight into the strictly lower triangle and the diagonal stays U's.
        const pivot = work.data[k][k];
        for (k + 1..n) |row| {
            const factor = work.data[k][row] / pivot;
            work.data[k][row] = factor;
            if (factor == 0) continue;
            for (k + 1..n) |column| {
                work.data[column][row] = work.data[column][row] - factor * work.data[column][k];
            }
        }
    }

    // A matrix that passes the pivot test can still drive an entry of a factor
    // past the format's range. Reporting that as "no factorization" is the
    // honest answer; returning an infinity is not.
    allFinite(M, work) orelse return null;

    return LuOf(M){
        .lu = work,
        .permutation = permutation,
        .parity = parity,
    };
}

/// L with its unit diagonal written in and its strict upper triangle zeroed,
/// which is the factor `f` packs into the lower triangle.
fn unitLowerOf(comptime M: type, f: anytype) M {
    const n = comptime M.rows;
    var lower: M = undefined;
    inline for (0..n) |row| {
        inline for (0..n) |column| {
            lower.data[column][row] = if (row >= column) f.lu.data[column][row] else 0;
        }
    }
    inline for (0..n) |index| lower.data[index][index] = 1;
    return lower;
}

/// `f.lu` with its strict lower triangle zeroed, which is U on and above the
/// diagonal and nothing else.
fn upperOf(comptime M: type, f: anytype) M {
    const n = comptime M.rows;
    var upper: M = undefined;
    inline for (0..n) |row| {
        inline for (0..n) |column| {
            upper.data[column][row] = if (row <= column) f.lu.data[column][row] else 0;
        }
    }
    return upper;
}

/// Solves `A * x = b` for one column from an LU factorization: the rows of `b`
/// moved to where the factors put them, then forward substitution against L,
/// whose unit diagonal needs no division, then back substitution against U.
///
/// `null` means a stored pivot of U was zero or not finite, which a
/// constructor cannot produce and a hand-built factorization type can.
fn luSolveColumn(comptime M: type, f: anytype, given: [n_of(M)]M.Element) ?[n_of(M)]M.Element {
    @setFloatMode(.strict);
    const n = comptime M.rows;
    const T = M.Element;

    // The stored factors satisfy `P * A = L * U`, so solving `A * x = b` is
    // solving `L * U * x = P * b` and the RIGHT-HAND SIDE moves with the rows.
    // Skipping this is the failure that a solve of a matrix needing an
    // interchange gets, and it is why `reconstructLu` alone is not a solve
    // test.
    var b: [n]T = undefined;
    inline for (0..n) |position| b[position] = given[f.permutation[position]];
    var forward: [n]T = undefined;
    inline for (0..n) |row| {
        if (row == 0) {
            // L has a unit diagonal, so the first row is b itself.
            forward[0] = b[0];
        } else {
            var sum = f.lu.data[0][row] * forward[0];
            for (1..row) |column| sum = sum + f.lu.data[column][row] * forward[column];
            forward[row] = b[row] - sum;
        }
    }

    var solution: [n]T = undefined;
    // Back substitution visits the rows from the last up, but each row's sum
    // over the columns to its right still accumulates in INCREASING column
    // order, from the index-zero product.
    for (0..n) |step| {
        const row = n - 1 - step;
        const pivot = f.lu.data[row][row];
        if (!std.math.isFinite(pivot) or pivot == 0) return null;
        if (row + 1 == n) {
            solution[row] = forward[row] / pivot;
        } else {
            var sum = f.lu.data[row + 1][row] * solution[row + 1];
            for (row + 2..n) |column| sum = sum + f.lu.data[column][row] * solution[column];
            solution[row] = (forward[row] - sum) / pivot;
        }
    }
    return solution;
}

// ---------------------------------------------------------------------------
// Cholesky
// ---------------------------------------------------------------------------

/// The Cholesky factorization `A = L * transpose(L)` of a symmetric positive
/// definite real matrix, or `null`.
///
/// `null` means one of four things: an entry was not finite, the matrix was not
/// symmetric, a pivot was not positive to working precision, or a factor ran
/// past the format's range. The positivity test is a pivot test DURING the
/// factorization, so a matrix with a positive determinant and a negative
/// eigenvalue is refused, which a sign on a computed determinant would not do.
fn choleskyFactor(comptime M: type, value: M) ?CholeskyOf(M) {
    @setFloatMode(.strict);
    const n = comptime M.rows;
    const T = M.Element;
    const threshold = pivotThresholdOf(T);

    const scale = finiteScaleOf(M, value) orelse return null;

    // Asymmetry policy: exact symmetry or `null`. The finiteness check is
    // first, so a non-finite entry is reported as non-finite rather than as an
    // asymmetry it would also fail.
    inline for (0..n) |row| {
        inline for (0..n) |column| {
            if (row < column) {
                if (!(value.data[column][row] == value.data[row][column])) return null;
            }
        }
    }

    // Only the lower triangle is read, and only after the matrix is known to be
    // symmetric, which is what stops this from being the "use one triangle"
    // defect.
    var lower: M = undefined;
    inline for (0..n) |row| {
        inline for (0..n) |column| {
            lower.data[column][row] = if (row >= column) value.data[column][row] else 0;
        }
    }

    inline for (0..n) |k| {
        // The diagonal pivot is the entry minus the column's own contribution,
        // accumulated in increasing column order from the index-zero product.
        var subtracted: T = 0;
        if (k != 0) {
            subtracted = lower.data[0][k] * lower.data[0][k];
            for (1..k) |column| {
                subtracted = subtracted + lower.data[column][k] * lower.data[column][k];
            }
        }
        const pivot = lower.data[k][k] - subtracted;
        // A non-positive pivot means not positive definite, and the relative
        // threshold also refuses a pivot that is positive but has no
        // significant bits left in it. See the module comment.
        if (!(pivot > threshold * scale)) return null;
        const root = sqrtOf(T, pivot);
        lower.data[k][k] = root;

        for (k + 1..n) |row| {
            var total = lower.data[k][row];
            if (k != 0) {
                total = total - lower.data[0][row] * lower.data[0][k];
                for (1..k) |column| total = total - lower.data[column][row] * lower.data[column][k];
            }
            lower.data[k][row] = total / root;
        }
    }

    allFinite(M, lower) orelse return null;

    return CholeskyOf(M){ .lower = lower };
}

/// Solves `A * x = b` for one column from a Cholesky factorization, as two
/// substitutions against L and against its transpose.
fn choleskySolveColumn(comptime M: type, f: anytype, b: [n_of(M)]M.Element) ?[n_of(M)]M.Element {
    @setFloatMode(.strict);
    const n = comptime M.rows;
    const T = M.Element;

    var forward: [n]T = undefined;
    inline for (0..n) |row| {
        const pivot = f.lower.data[row][row];
        if (!std.math.isFinite(pivot) or pivot == 0) return null;
        if (row == 0) {
            forward[0] = b[0] / pivot;
        } else {
            var sum = f.lower.data[0][row] * forward[0];
            for (1..row) |column| sum = sum + f.lower.data[column][row] * forward[column];
            forward[row] = (b[row] - sum) / pivot;
        }
    }

    // The second substitution is against the TRANSPOSE of L, so row `row` of
    // that system sums the entries L[column][row] for the columns to its right:
    // `data[row][column]`, the mirror of the first substitution's
    // `data[column][row]`.
    var solution: [n]T = undefined;
    for (0..n) |step| {
        const row = n - 1 - step;
        const pivot = f.lower.data[row][row];
        if (!std.math.isFinite(pivot) or pivot == 0) return null;
        if (row + 1 == n) {
            solution[row] = forward[row] / pivot;
        } else {
            var sum = f.lower.data[row][row + 1] * solution[row + 1];
            for (row + 2..n) |column| sum = sum + f.lower.data[row][column] * solution[column];
            solution[row] = (forward[row] - sum) / pivot;
        }
    }
    return solution;
}

// ---------------------------------------------------------------------------
// Named operations: LU
//
// The contract names `LuMat2f`, `LuMat3f`, and `LuMat4f`, and the size-paired
// binary64 siblings are here too: `zmath.linalg` consumes the whole family and
// a caller should not have to reach for `Mat` and factor it by hand.
// ---------------------------------------------------------------------------

/// The LU factorization of a square real matrix, by Doolittle elimination with
/// partial pivoting, or `null`.
///
/// The size and the precision are carried by the argument and the return type
/// is computed from it, which is how one name reaches the six factorization
/// types: `lu` applied to a `Mat2f`, `Mat3f`, `Mat4f`, `Mat2d`, `Mat3d`, or
/// `Mat4d` has type `?LuMat2f`, `?LuMat3f`, `?LuMat4f`, `?LuMat2d`, `?LuMat3d`,
/// or `?LuMat4d` respectively. Any other matrix type is a compile error.
///
/// `null` means one of three things: an entry was not finite, a pivot fell at
/// or below `matrix.inverse_pivot_threshold_f` or
/// `matrix.inverse_pivot_threshold_d` times the largest entry magnitude of the
/// input, or the factors themselves ran past the format's range. See the
/// module comment for the threshold's derivation and for why there is no
/// epsilon argument.
pub fn lu(value: anytype) ?LuOf(@TypeOf(value)) {
    return luFactor(@TypeOf(value), value);
}

/// The matrix `f` was factorized from: `L * U` with row `permutation[k]` of the
/// result taken from row `k` of the product.
///
/// This is the faithfulness check of the whole module. It does not return the
/// input it was handed, because the input is not stored: it returns what the
/// stored factors say the input was, and a caller that gets the original back
/// has proved that `permutation`, `parity`, and the packed factors all agree.
pub fn reconstructLu(f: anytype) MatrixOf(@TypeOf(f)) {
    const M = MatrixOf(@TypeOf(f));
    const n = comptime M.rows;
    @setFloatMode(.strict);
    const product = matrix.mul(n, unitLowerOf(M, f), upperOf(M, f));
    var out: M = undefined;
    inline for (0..n) |position| {
        inline for (0..n) |column| {
            // The stored factors satisfy `P * A = L * U`, so row `position` of
            // the product is row `permutation[position]` of `A`.
            out.data[column][f.permutation[position]] = product.data[column][position];
        }
    }
    return out;
}

/// The determinant of the matrix `f` was factorized from:
/// `parity * product(diag(U))`, accumulated in increasing diagonal index from
/// the index-zero entry.
///
/// `det(P * A) = det(L) * det(U)` gives `parity * det(A) = product(diag(U))`,
/// and the sign of a permutation is its own inverse, so the factor is `parity`
/// and not its reciprocal. An odd number of interchanges therefore contributes
/// a factor of `-1`, and a matrix that needed no interchange the neutral one.
pub fn determinantLu(f: anytype) MatrixOf(@TypeOf(f)).Element {
    const M = MatrixOf(@TypeOf(f));
    const n = comptime M.rows;
    @setFloatMode(.strict);
    var product = f.lu.data[0][0];
    for (1..n) |index| product = product * f.lu.data[index][index];
    return if (f.parity == -1) -product else product;
}

/// The solution `x` of `A * x = b`, where `A` is the matrix `f` was factorized
/// from, or `null` if a stored pivot of U is zero or not finite.
pub fn solveLu(f: anytype, rhs: anytype) ?MatrixOf(@TypeOf(f)).vector_type {
    const M = MatrixOf(@TypeOf(f));
    const V = M.vector_type;
    if (@TypeOf(rhs) != V) {
        @compileError("zmath.linalg.solveLu: the right-hand side must be a " ++ @typeName(V));
    }
    const b: [M.rows]M.Element = rhs.toArray();
    return V.fromArray(luSolveColumn(M, f, b) orelse return null);
}

/// The solutions of `A * X = B` for every COLUMN of `B` at once, where `A` is
/// the matrix `f` was factorized from.
///
/// One substitution pair against the factors serves every column, which is the
/// whole reason to ask for several right-hand sides rather than calling
/// `solveLu` in a loop. `B` is an n by m matrix whose j-th column is the j-th
/// right-hand side, so the result has the same shape.
pub fn solveLuMultiple(f: anytype, rhs: anytype) ?RightHandSidesOf(@TypeOf(f), @TypeOf(rhs)) {
    const F = @TypeOf(f);
    const M = MatrixOf(F);
    const R = @TypeOf(rhs);
    const Out = RightHandSidesOf(F, R);
    const n = comptime M.rows;
    @setFloatMode(.strict);
    var out: Out = undefined;
    inline for (0..R.columns) |column| {
        var b: [n]M.Element = undefined;
        inline for (0..n) |row| b[row] = rhs.data[column][row];
        const solved = luSolveColumn(M, f, b) orelse return null;
        inline for (0..n) |row| out.data[column][row] = solved[row];
    }
    return out;
}

// ---------------------------------------------------------------------------
// Named operations: Cholesky
// ---------------------------------------------------------------------------

/// The Cholesky factorization `A = L * transpose(L)` of a square real matrix,
/// or `null`.
///
/// The size and the precision are carried by the argument and the return type
/// is computed from it, exactly as for `lu`.
///
/// `null` means one of four things: an entry was not finite, the matrix was not
/// SYMMETRIC, a pivot was not positive to working precision, or a factor ran
/// past the format's range. The matrix is refused rather than symmetrized and
/// rather than read from one triangle only; the module comment gives the
/// policy and its reason. The positivity test is a pivot test DURING the
/// factorization, so a matrix with a positive determinant and a negative
/// eigenvalue is refused, which a sign on a determinant computed afterwards
/// would not catch.
pub fn cholesky(value: anytype) ?CholeskyOf(@TypeOf(value)) {
    return choleskyFactor(@TypeOf(value), value);
}

/// The matrix `f` was factorized from: `L * transpose(L)`.
///
/// As with `reconstructLu`, this returns what the stored factor says the input
/// was, and returning the input is what proves the factor is faithful.
pub fn reconstructCholesky(f: anytype) MatrixOf(@TypeOf(f)) {
    const M = MatrixOf(@TypeOf(f));
    const n = comptime M.rows;
    @setFloatMode(.strict);
    return matrix.mul(n, f.lower, f.lower.transpose());
}

/// The determinant of the matrix `f` was factorized from: the product of the
/// squared diagonal of L, accumulated in increasing diagonal index from the
/// index-zero square.
///
/// `A = L * transpose(L)` gives `det(A) = det(L) * det(transpose(L)) =
/// det(L)^2`, and `det(L)` is the product of L's diagonal. Every factor of L is
/// positive, so the product of the squares is nonnegative whatever the sign of
/// the determinant would be, which is correct: a positive-definite matrix has a
/// positive determinant.
pub fn determinantCholesky(f: anytype) MatrixOf(@TypeOf(f)).Element {
    const M = MatrixOf(@TypeOf(f));
    const n = comptime M.rows;
    @setFloatMode(.strict);
    var product = f.lower.data[0][0] * f.lower.data[0][0];
    for (1..n) |index| {
        const entry = f.lower.data[index][index];
        product = product * (entry * entry);
    }
    return product;
}

/// The solution `x` of `A * x = b`, where `A` is the symmetric positive
/// definite matrix `f` was factorized from, or `null` if a stored pivot of L is
/// zero or not finite.
pub fn solveCholesky(f: anytype, rhs: anytype) ?MatrixOf(@TypeOf(f)).vector_type {
    const M = MatrixOf(@TypeOf(f));
    const V = M.vector_type;
    if (@TypeOf(rhs) != V) {
        @compileError("zmath.linalg.solveCholesky: the right-hand side must be a " ++ @typeName(V));
    }
    const b: [M.rows]M.Element = rhs.toArray();
    return V.fromArray(choleskySolveColumn(M, f, b) orelse return null);
}

/// The solutions of `A * X = B` for every COLUMN of `B` at once, where `A` is
/// the symmetric positive definite matrix `f` was factorized from.
pub fn solveCholeskyMultiple(f: anytype, rhs: anytype) ?RightHandSidesOf(@TypeOf(f), @TypeOf(rhs)) {
    const F = @TypeOf(f);
    const M = MatrixOf(F);
    const R = @TypeOf(rhs);
    const Out = RightHandSidesOf(F, R);
    const n = comptime M.rows;
    @setFloatMode(.strict);
    var out: Out = undefined;
    inline for (0..R.columns) |column| {
        var b: [n]M.Element = undefined;
        inline for (0..n) |row| b[row] = rhs.data[column][row];
        const solved = choleskySolveColumn(M, f, b) orelse return null;
        inline for (0..n) |row| out.data[column][row] = solved[row];
    }
    return out;
}
