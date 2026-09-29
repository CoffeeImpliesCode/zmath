//! Vector (`@Vector`) overloads for the remaining elementwise IEEE and
//! miscellaneous scalar operations.
//!
//! Every operation here applies the *identical* accepted scalar kernel to each
//! lane, so lane `i` of the result is bit-for-bit equal to the scalar result for
//! lane `i` of the input. That is what preserves the schema-3 `max_ulp = 0`
//! contract and the immutable corpus without new acceptance evidence: the vector
//! path is a per-lane application of the same scalar code, never a re-derived
//! polynomial.
//!
//! Most of these operations are float-to-float, so the result element type is
//! the lane type. `ilogb` and `llogb` instead return `@Vector(N, i32)` and
//! `@Vector(N, i64)`, and `totalOrder` and `totalOrderMag` return a
//! `@Vector(N, bool)` mask. Those are why the shared helpers in `vector_core`
//! take an explicit result type `R` rather than assuming the lane type.
//!
//! Result structs in `scalar_c23_types.zig` are bare `struct` values with auto
//! layout. They are never `extern` or `packed`, so `@bitCast` and `@offsetOf`
//! do not apply. Each field is assigned explicitly, never reinterpreted. No
//! operation in this module returns a struct, but the constraint holds for the
//! struct-returning modules built on the same core.
//!
//! `nan` is deliberately absent. Its type argument is a `comptime` type and its
//! tag is a slice, so it has no lane to map over; the vector spelling of a NaN
//! is `@splat` at the call site.
//!
//! Names follow `op + type + v`, exactly as in `vector.zig`. None of the
//! operations here collides with the hyperbolic generics `sinhv`, `coshv`,
//! `tanhv`, `asinhv`, `acoshv`, and `atanhv`, so every operation also has a
//! binary16 specialization. The six extrema entry points capitalise their
//! word parts, as the `totalOrderMag` entry points do, while their public
//! `root.zig` names keep the underscores. Both spellings stay distinct from
//! the already vectorized `fmaximumv` and `fminimumv`.

const std = @import("std");
const core = @import("vector_core.zig");
const scalar = @import("scalar.zig");
const ieee = @import("ieee.zig");

// ---------------------------------------------------------------------------
// Local width helper
// ---------------------------------------------------------------------------

/// Returns the lane count of the vector type `X`.
///
/// The integer- and boolean-returning entry points need the width in their
/// return type, because their result element type is not the lane type. The
/// width is derived from the type alone, so it is known at comptime.
fn widthOf(comptime X: type) usize {
    return switch (@typeInfo(X)) {
        .vector => |vec| vec.len,
        else => @compileError("zmath vector operations require a @Vector argument"),
    };
}

// ---------------------------------------------------------------------------
// Rounding
// ---------------------------------------------------------------------------

/// Computes the round of every lane in `x` with the scalar kernel for its lane type.
pub fn roundv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return core.mapUnaryR(f16, f16, scalar.roundh, x);
    if (T == f32) return core.mapUnaryR(f32, f32, scalar.roundf, x);
    if (T == f64) return core.mapUnaryR(f64, f64, scalar.roundd, x);
    @compileError("zmath.roundv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the round of every lane in `x` with the binary16 scalar kernel.
pub fn roundhv(x: anytype) @TypeOf(x) {
    return core.mapUnaryR(f16, f16, scalar.roundh, x);
}

/// Computes the round of every lane in `x` with the binary32 scalar kernel.
pub fn roundfv(x: anytype) @TypeOf(x) {
    return core.mapUnaryR(f32, f32, scalar.roundf, x);
}

/// Computes the round of every lane in `x` with the binary64 scalar kernel.
pub fn rounddv(x: anytype) @TypeOf(x) {
    return core.mapUnaryR(f64, f64, scalar.roundd, x);
}

// ---------------------------------------------------------------------------
// IEEE representation: nextUp, nextDown, getPayload
// ---------------------------------------------------------------------------

/// Computes the nextUp of every lane in `x` with the scalar kernel for its lane type.
pub fn nextUpv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return core.mapUnaryR(f16, f16, ieee.nextUph, x);
    if (T == f32) return core.mapUnaryR(f32, f32, ieee.nextUpf, x);
    if (T == f64) return core.mapUnaryR(f64, f64, ieee.nextUpd, x);
    @compileError("zmath.nextUpv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the nextUp of every lane in `x` with the binary16 scalar kernel.
pub fn nextUphv(x: anytype) @TypeOf(x) {
    return core.mapUnaryR(f16, f16, ieee.nextUph, x);
}

/// Computes the nextUp of every lane in `x` with the binary32 scalar kernel.
pub fn nextUpfv(x: anytype) @TypeOf(x) {
    return core.mapUnaryR(f32, f32, ieee.nextUpf, x);
}

/// Computes the nextUp of every lane in `x` with the binary64 scalar kernel.
pub fn nextUpdv(x: anytype) @TypeOf(x) {
    return core.mapUnaryR(f64, f64, ieee.nextUpd, x);
}

/// Computes the nextDown of every lane in `x` with the scalar kernel for its lane type.
pub fn nextDownv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return core.mapUnaryR(f16, f16, ieee.nextDownh, x);
    if (T == f32) return core.mapUnaryR(f32, f32, ieee.nextDownf, x);
    if (T == f64) return core.mapUnaryR(f64, f64, ieee.nextDownd, x);
    @compileError("zmath.nextDownv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the nextDown of every lane in `x` with the binary16 scalar kernel.
pub fn nextDownhv(x: anytype) @TypeOf(x) {
    return core.mapUnaryR(f16, f16, ieee.nextDownh, x);
}

/// Computes the nextDown of every lane in `x` with the binary32 scalar kernel.
pub fn nextDownfv(x: anytype) @TypeOf(x) {
    return core.mapUnaryR(f32, f32, ieee.nextDownf, x);
}

/// Computes the nextDown of every lane in `x` with the binary64 scalar kernel.
pub fn nextDowndv(x: anytype) @TypeOf(x) {
    return core.mapUnaryR(f64, f64, ieee.nextDownd, x);
}

/// Computes the getPayload of every lane in `x` with the scalar kernel for its lane type.
pub fn getPayloadv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return core.mapUnaryR(f16, f16, ieee.getPayloadh, x);
    if (T == f32) return core.mapUnaryR(f32, f32, ieee.getPayloadf, x);
    if (T == f64) return core.mapUnaryR(f64, f64, ieee.getPayloadd, x);
    @compileError("zmath.getPayloadv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the getPayload of every lane in `x` with the binary16 scalar kernel.
pub fn getPayloadhv(x: anytype) @TypeOf(x) {
    return core.mapUnaryR(f16, f16, ieee.getPayloadh, x);
}

/// Computes the getPayload of every lane in `x` with the binary32 scalar kernel.
pub fn getPayloadfv(x: anytype) @TypeOf(x) {
    return core.mapUnaryR(f32, f32, ieee.getPayloadf, x);
}

/// Computes the getPayload of every lane in `x` with the binary64 scalar kernel.
pub fn getPayloaddv(x: anytype) @TypeOf(x) {
    return core.mapUnaryR(f64, f64, ieee.getPayloadd, x);
}

// ---------------------------------------------------------------------------
// IEEE representation: copySign, nextAfter
// ---------------------------------------------------------------------------

/// Computes the copySign of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn copySignv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return core.mapBinaryR(f16, f16, ieee.copySignh, a, b);
    if (T == f32) return core.mapBinaryR(f32, f32, ieee.copySignf, a, b);
    if (T == f64) return core.mapBinaryR(f64, f64, ieee.copySignd, a, b);
    @compileError("zmath.copySignv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the copySign of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn copySignhv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f16, f16, ieee.copySignh, a, b);
}

/// Computes the copySign of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn copySignfv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f32, f32, ieee.copySignf, a, b);
}

/// Computes the copySign of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn copySigndv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f64, f64, ieee.copySignd, a, b);
}

/// Computes the nextAfter of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn nextAfterv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return core.mapBinaryR(f16, f16, ieee.nextAfterh, a, b);
    if (T == f32) return core.mapBinaryR(f32, f32, ieee.nextAfterf, a, b);
    if (T == f64) return core.mapBinaryR(f64, f64, ieee.nextAfterd, a, b);
    @compileError("zmath.nextAfterv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the nextAfter of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn nextAfterhv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f16, f16, ieee.nextAfterh, a, b);
}

/// Computes the nextAfter of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn nextAfterfv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f32, f32, ieee.nextAfterf, a, b);
}

/// Computes the nextAfter of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn nextAfterdv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f64, f64, ieee.nextAfterd, a, b);
}

// ---------------------------------------------------------------------------
// Extrema aliases
// ---------------------------------------------------------------------------

/// Computes the fmaximum_mag of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn fmaximumMagv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return core.mapBinaryR(f16, f16, scalar.fmaximum_magh, a, b);
    if (T == f32) return core.mapBinaryR(f32, f32, scalar.fmaximum_magf, a, b);
    if (T == f64) return core.mapBinaryR(f64, f64, scalar.fmaximum_magd, a, b);
    @compileError("zmath.fmaximumMagv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the fmaximum_mag of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn fmaximumMaghv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f16, f16, scalar.fmaximum_magh, a, b);
}

/// Computes the fmaximum_mag of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn fmaximumMagfv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f32, f32, scalar.fmaximum_magf, a, b);
}

/// Computes the fmaximum_mag of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn fmaximumMagdv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f64, f64, scalar.fmaximum_magd, a, b);
}

/// Computes the fminimum_mag of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn fminimumMagv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return core.mapBinaryR(f16, f16, scalar.fminimum_magh, a, b);
    if (T == f32) return core.mapBinaryR(f32, f32, scalar.fminimum_magf, a, b);
    if (T == f64) return core.mapBinaryR(f64, f64, scalar.fminimum_magd, a, b);
    @compileError("zmath.fminimumMagv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the fminimum_mag of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn fminimumMaghv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f16, f16, scalar.fminimum_magh, a, b);
}

/// Computes the fminimum_mag of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn fminimumMagfv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f32, f32, scalar.fminimum_magf, a, b);
}

/// Computes the fminimum_mag of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn fminimumMagdv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f64, f64, scalar.fminimum_magd, a, b);
}

/// Computes the fmaximum_num of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn fmaximumNumv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return core.mapBinaryR(f16, f16, scalar.fmaximum_numh, a, b);
    if (T == f32) return core.mapBinaryR(f32, f32, scalar.fmaximum_numf, a, b);
    if (T == f64) return core.mapBinaryR(f64, f64, scalar.fmaximum_numd, a, b);
    @compileError("zmath.fmaximumNumv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the fmaximum_num of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn fmaximumNumhv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f16, f16, scalar.fmaximum_numh, a, b);
}

/// Computes the fmaximum_num of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn fmaximumNumfv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f32, f32, scalar.fmaximum_numf, a, b);
}

/// Computes the fmaximum_num of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn fmaximumNumdv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f64, f64, scalar.fmaximum_numd, a, b);
}

/// Computes the fminimum_num of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn fminimumNumv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return core.mapBinaryR(f16, f16, scalar.fminimum_numh, a, b);
    if (T == f32) return core.mapBinaryR(f32, f32, scalar.fminimum_numf, a, b);
    if (T == f64) return core.mapBinaryR(f64, f64, scalar.fminimum_numd, a, b);
    @compileError("zmath.fminimumNumv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the fminimum_num of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn fminimumNumhv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f16, f16, scalar.fminimum_numh, a, b);
}

/// Computes the fminimum_num of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn fminimumNumfv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f32, f32, scalar.fminimum_numf, a, b);
}

/// Computes the fminimum_num of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn fminimumNumdv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f64, f64, scalar.fminimum_numd, a, b);
}

/// Computes the fmaximum_mag_num of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn fmaximumMagNumv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return core.mapBinaryR(f16, f16, scalar.fmaximum_mag_numh, a, b);
    if (T == f32) return core.mapBinaryR(f32, f32, scalar.fmaximum_mag_numf, a, b);
    if (T == f64) return core.mapBinaryR(f64, f64, scalar.fmaximum_mag_numd, a, b);
    @compileError("zmath.fmaximumMagNumv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the fmaximum_mag_num of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn fmaximumMagNumhv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f16, f16, scalar.fmaximum_mag_numh, a, b);
}

/// Computes the fmaximum_mag_num of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn fmaximumMagNumfv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f32, f32, scalar.fmaximum_mag_numf, a, b);
}

/// Computes the fmaximum_mag_num of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn fmaximumMagNumdv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f64, f64, scalar.fmaximum_mag_numd, a, b);
}

/// Computes the fminimum_mag_num of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn fminimumMagNumv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return core.mapBinaryR(f16, f16, scalar.fminimum_mag_numh, a, b);
    if (T == f32) return core.mapBinaryR(f32, f32, scalar.fminimum_mag_numf, a, b);
    if (T == f64) return core.mapBinaryR(f64, f64, scalar.fminimum_mag_numd, a, b);
    @compileError("zmath.fminimumMagNumv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the fminimum_mag_num of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn fminimumMagNumhv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f16, f16, scalar.fminimum_mag_numh, a, b);
}

/// Computes the fminimum_mag_num of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn fminimumMagNumfv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f32, f32, scalar.fminimum_mag_numf, a, b);
}

/// Computes the fminimum_mag_num of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn fminimumMagNumdv(a: anytype, b: anytype) @TypeOf(a) {
    return core.mapBinaryR(f64, f64, scalar.fminimum_mag_numd, a, b);
}

// ---------------------------------------------------------------------------
// Integer-returning exponent queries
// ---------------------------------------------------------------------------

/// Computes the binary exponent of every lane in `x` as an i32 with the scalar kernel for its lane type.
pub fn ilogbv(x: anytype) @Vector(widthOf(@TypeOf(x)), i32) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return core.mapUnaryR(f16, i32, scalar.ilogbh, x);
    if (T == f32) return core.mapUnaryR(f32, i32, scalar.ilogbf, x);
    if (T == f64) return core.mapUnaryR(f64, i32, scalar.ilogbd, x);
    @compileError("zmath.ilogbv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the binary exponent of every lane in `x` as an i32 with the binary16 scalar kernel.
pub fn ilogbhv(x: anytype) @Vector(widthOf(@TypeOf(x)), i32) {
    return core.mapUnaryR(f16, i32, scalar.ilogbh, x);
}

/// Computes the binary exponent of every lane in `x` as an i32 with the binary32 scalar kernel.
pub fn ilogbfv(x: anytype) @Vector(widthOf(@TypeOf(x)), i32) {
    return core.mapUnaryR(f32, i32, scalar.ilogbf, x);
}

/// Computes the binary exponent of every lane in `x` as an i32 with the binary64 scalar kernel.
pub fn ilogbdv(x: anytype) @Vector(widthOf(@TypeOf(x)), i32) {
    return core.mapUnaryR(f64, i32, scalar.ilogbd, x);
}

/// Computes the binary exponent of every lane in `x` as an i64 with the scalar kernel for its lane type.
pub fn llogbv(x: anytype) @Vector(widthOf(@TypeOf(x)), i64) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return core.mapUnaryR(f16, i64, scalar.llogbh, x);
    if (T == f32) return core.mapUnaryR(f32, i64, scalar.llogbf, x);
    if (T == f64) return core.mapUnaryR(f64, i64, scalar.llogbd, x);
    @compileError("zmath.llogbv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the binary exponent of every lane in `x` as an i64 with the binary16 scalar kernel.
pub fn llogbhv(x: anytype) @Vector(widthOf(@TypeOf(x)), i64) {
    return core.mapUnaryR(f16, i64, scalar.llogbh, x);
}

/// Computes the binary exponent of every lane in `x` as an i64 with the binary32 scalar kernel.
pub fn llogbfv(x: anytype) @Vector(widthOf(@TypeOf(x)), i64) {
    return core.mapUnaryR(f32, i64, scalar.llogbf, x);
}

/// Computes the binary exponent of every lane in `x` as an i64 with the binary64 scalar kernel.
pub fn llogbdv(x: anytype) @Vector(widthOf(@TypeOf(x)), i64) {
    return core.mapUnaryR(f64, i64, scalar.llogbd, x);
}

// ---------------------------------------------------------------------------
// Boolean-returning ordering comparisons
// ---------------------------------------------------------------------------

/// Computes the totalOrder of corresponding lanes of `a` and `b` as a bool
/// mask with the scalar kernel for their lane type.
pub fn totalOrderv(a: anytype, b: anytype) @Vector(widthOf(@TypeOf(a)), bool) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return core.mapBinaryR(f16, bool, ieee.totalOrderh, a, b);
    if (T == f32) return core.mapBinaryR(f32, bool, ieee.totalOrderf, a, b);
    if (T == f64) return core.mapBinaryR(f64, bool, ieee.totalOrderd, a, b);
    @compileError("zmath.totalOrderv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the totalOrder of corresponding lanes of `a` and `b` as a bool mask with the binary16 scalar kernel.
pub fn totalOrderhv(a: anytype, b: anytype) @Vector(widthOf(@TypeOf(a)), bool) {
    return core.mapBinaryR(f16, bool, ieee.totalOrderh, a, b);
}

/// Computes the totalOrder of corresponding lanes of `a` and `b` as a bool mask with the binary32 scalar kernel.
pub fn totalOrderfv(a: anytype, b: anytype) @Vector(widthOf(@TypeOf(a)), bool) {
    return core.mapBinaryR(f32, bool, ieee.totalOrderf, a, b);
}

/// Computes the totalOrder of corresponding lanes of `a` and `b` as a bool mask with the binary64 scalar kernel.
pub fn totalOrderdv(a: anytype, b: anytype) @Vector(widthOf(@TypeOf(a)), bool) {
    return core.mapBinaryR(f64, bool, ieee.totalOrderd, a, b);
}

/// Computes the totalOrderMag of corresponding lanes of `a` and `b` as a bool
/// mask with the scalar kernel for their lane type.
pub fn totalOrderMagv(a: anytype, b: anytype) @Vector(widthOf(@TypeOf(a)), bool) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return core.mapBinaryR(f16, bool, ieee.totalOrderMagh, a, b);
    if (T == f32) return core.mapBinaryR(f32, bool, ieee.totalOrderMagf, a, b);
    if (T == f64) return core.mapBinaryR(f64, bool, ieee.totalOrderMagd, a, b);
    @compileError("zmath.totalOrderMagv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the totalOrderMag of corresponding lanes of `a` and `b` as a bool mask with the binary16 scalar kernel.
pub fn totalOrderMaghv(a: anytype, b: anytype) @Vector(widthOf(@TypeOf(a)), bool) {
    return core.mapBinaryR(f16, bool, ieee.totalOrderMagh, a, b);
}

/// Computes the totalOrderMag of corresponding lanes of `a` and `b` as a bool mask with the binary32 scalar kernel.
pub fn totalOrderMagfv(a: anytype, b: anytype) @Vector(widthOf(@TypeOf(a)), bool) {
    return core.mapBinaryR(f32, bool, ieee.totalOrderMagf, a, b);
}

/// Computes the totalOrderMag of corresponding lanes of `a` and `b` as a bool mask with the binary64 scalar kernel.
pub fn totalOrderMagdv(a: anytype, b: anytype) @Vector(widthOf(@TypeOf(a)), bool) {
    return core.mapBinaryR(f64, bool, ieee.totalOrderMagd, a, b);
}

comptime {
    _ = roundv;
    _ = roundhv;
    _ = roundfv;
    _ = rounddv;
    _ = nextUpv;
    _ = nextUphv;
    _ = nextUpfv;
    _ = nextUpdv;
    _ = nextDownv;
    _ = nextDownhv;
    _ = nextDownfv;
    _ = nextDowndv;
    _ = getPayloadv;
    _ = getPayloadhv;
    _ = getPayloadfv;
    _ = getPayloaddv;
    _ = copySignv;
    _ = copySignhv;
    _ = copySignfv;
    _ = copySigndv;
    _ = nextAfterv;
    _ = nextAfterhv;
    _ = nextAfterfv;
    _ = nextAfterdv;
    _ = fmaximumMagv;
    _ = fmaximumMaghv;
    _ = fmaximumMagfv;
    _ = fmaximumMagdv;
    _ = fminimumMagv;
    _ = fminimumMaghv;
    _ = fminimumMagfv;
    _ = fminimumMagdv;
    _ = fmaximumNumv;
    _ = fmaximumNumhv;
    _ = fmaximumNumfv;
    _ = fmaximumNumdv;
    _ = fminimumNumv;
    _ = fminimumNumhv;
    _ = fminimumNumfv;
    _ = fminimumNumdv;
    _ = fmaximumMagNumv;
    _ = fmaximumMagNumhv;
    _ = fmaximumMagNumfv;
    _ = fmaximumMagNumdv;
    _ = fminimumMagNumv;
    _ = fminimumMagNumhv;
    _ = fminimumMagNumfv;
    _ = fminimumMagNumdv;
    _ = ilogbv;
    _ = ilogbhv;
    _ = ilogbfv;
    _ = ilogbdv;
    _ = llogbv;
    _ = llogbhv;
    _ = llogbfv;
    _ = llogbdv;
    _ = totalOrderv;
    _ = totalOrderhv;
    _ = totalOrderfv;
    _ = totalOrderdv;
    _ = totalOrderMagv;
    _ = totalOrderMaghv;
    _ = totalOrderMagfv;
    _ = totalOrderMagdv;
}
