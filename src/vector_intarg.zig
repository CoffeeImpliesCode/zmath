//! Vector (`@Vector`) overloads for the scalar operations that take a
//! non-float argument.
//!
//! Every operation here applies the *identical* accepted scalar kernel to each
//! lane, so lane `i` of the result is bit-for-bit equal to the scalar result for
//! lane `i` of the input and of the extra argument. That is what preserves the
//! schema-3 `max_ulp = 0` contract and the immutable corpus without new
//! acceptance evidence: the vector path is a per-lane application of the same
//! scalar code, never a re-derived polynomial.
//!
//! ## Result shape
//!
//! Zig forbids a struct as a `@Vector` element type, so the result shape is
//! not uniform across the eleven operations. This is a language limitation, not
//! an oversight.
//!
//! - The eight operations whose result is a plain float — `ldexp`, `scalbn`,
//!   `scalbln`, `nextToward`, `fromfp`, `ufromfp`, `fromfpx` and `ufromfpx` —
//!   return a real `@Vector(N, T)` and go through the shared helpers in
//!   `vector_core.zig`.
//! - The three operations whose result is a struct — `pown`, `compoundn` and
//!   `rootn`, all of which yield `PowerResult(T)` — return a lane-ordered array
//!   `[N]PowerResult(T)` instead, filled by the private per-lane helpers in this
//!   module. A `vector_core` helper cannot be used, because every one of them
//!   returns a `@Vector`.
//!
//! The per-lane bit-identity guarantee is identical in both shapes, so the
//! `max_ulp = 0` contract and the immutable corpus hold with no new acceptance
//! evidence either way. The array shape yields no throughput advantage over a
//! scalar loop, because the kernel is applied once per lane either way; its
//! value is a uniform lane-parallel API, not SIMD acceleration.
//!
//! `PowerResult(T)`, like every type in `scalar_c23_types.zig`, is a bare struct
//! with auto layout, never `extern` or `packed`, so `@bitCast` and `@offsetOf`
//! do not apply. Each lane is the kernel's own returned value, copied whole;
//! every field is assigned explicitly, never reinterpreted.
//!
//! ## Extra argument
//!
//! - A plain runtime scalar broadcast to every lane: the `i32` exponent of
//!   `ldexp` and `scalbn`, the `i64` exponent of `scalbln` and the power
//!   family, the `f64` `toward` of `nextToward`, and the `FromFpRound` direction
//!   of the `fromfp` family.
//! - A per-lane integer vector, spelled with the `Lanes` suffix, holding a
//!   `@Vector(N, i32)` or `@Vector(N, i64)` of exponents.
//! - Two broadcast scalars: the `FromFpRound` direction and the `u32` width of
//!   the `fromfp` family.
//!
//! `nextToward` has no `Lanes` form. Its `toward` argument is always `f64`, never
//! the lane type, and each per-format kernel applies its own internal tolerance
//! (binary16 42, binary32 29, binary64 0), so a per-lane spelling would have no
//! single argument type to check.
//!
//! ## Naming
//!
//! Naming follows `vector.zig`: `OPv` picks the kernel from the lane type, and
//! `OPfv`, `OPdv` and `OPhv` pin binary32, binary64 and binary16. None of these
//! eleven operations collides with the hyperbolic generics `sinhv`, `coshv`,
//! `tanhv`, `asinhv`, `acoshv` and `atanhv`, so all four spellings are present
//! for each of them.

const core = @import("vector_core.zig");
const ieee = @import("ieee.zig");
const scalar = @import("scalar.zig");
const types = @import("scalar_c23_types.zig");

// ---------------------------------------------------------------------------
// Width and lane-ordered power result
//
// Zig forbids a struct as a `@Vector` element type, so the three power
// operations cannot return a `@Vector`. They return a lane-ordered array
// instead, filled by the per-lane helpers below. Those helpers repeat what
// `vector_core` does, because every `vector_core` helper returns a `@Vector`
// and therefore cannot carry a struct result.
//
// `widthOf`, `LaneOf` and `PowerVector` are public so that the top-level
// dispatch in `root.zig` can spell the array return type of `pown`, `compoundn`
// and `rootn` for a vector argument. The existing `PowerResult(@TypeOf(value))`
// there is the scalar shape and is wrong for a vector.
// ---------------------------------------------------------------------------

/// Returns the lane count of the vector type `X`.
pub fn widthOf(comptime X: type) usize {
    return switch (@typeInfo(X)) {
        .vector => |v| v.len,
        else => @compileError("zmath vector operations require a @Vector argument"),
    };
}

/// Returns the lane type of the vector type `X`.
pub fn LaneOf(comptime X: type) type {
    return switch (@typeInfo(X)) {
        .vector => |v| v.child,
        else => @compileError("zmath vector operations require a @Vector argument"),
    };
}

/// The lane-ordered `[N]PowerResult(T)` result type for the lane type `T` and
/// the lane count of the vector type `X`.
pub fn PowerVector(comptime T: type, comptime X: type) type {
    return [widthOf(X)]types.PowerResult(T);
}

/// Applies `f` to every lane of `x` with the broadcast runtime scalar `e`, in
/// lane order.
///
/// Each lane is the kernel's own returned struct, copied whole. Nothing is
/// reinterpreted: these result types are bare structs with auto layout, so
/// `@bitCast` and `@offsetOf` do not apply.
fn mapPowerBroadcast(
    comptime T: type,
    comptime I: type,
    comptime f: fn (T, I) types.PowerResult(T),
    e: I,
    x: anytype,
) PowerVector(T, @TypeOf(x)) {
    var out: PowerVector(T, @TypeOf(x)) = undefined;
    const n = comptime widthOf(@TypeOf(x));
    inline for (0..n) |i| out[i] = f(x[i], e);
    return out;
}

/// Applies `f` to corresponding lanes of `x` and the per-lane integer `e`, in
/// lane order.
fn mapPowerLanes(
    comptime T: type,
    comptime I: type,
    comptime f: fn (T, I) types.PowerResult(T),
    x: anytype,
    e: anytype,
) PowerVector(T, @TypeOf(x)) {
    const n = comptime widthOf(@TypeOf(x));
    if (comptime widthOf(@TypeOf(e)) != n) {
        @compileError("zmath vector per-lane operations require both arguments to have the same width");
    }
    if (comptime LaneOf(@TypeOf(e)) != I) {
        @compileError("zmath vector lane type does not match the selected scalar kernel");
    }
    var out: PowerVector(T, @TypeOf(x)) = undefined;
    inline for (0..n) |i| out[i] = f(x[i], e[i]);
    return out;
}

// ---------------------------------------------------------------------------
// ldexp
// ---------------------------------------------------------------------------

/// Scales every lane in `x` by `2^exponent` with the scalar kernel for its lane
/// type, broadcasting the `i32` `exponent` to every lane.
pub fn ldexpv(x: anytype, exponent: i32) @TypeOf(x) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return core.mapUnaryExtra(f16, f16, i32, scalar.ldexph, exponent, x);
    if (lane == f32) return core.mapUnaryExtra(f32, f32, i32, scalar.ldexpf, exponent, x);
    if (lane == f64) return core.mapUnaryExtra(f64, f64, i32, scalar.ldexpd, exponent, x);
    @compileError("zmath.ldexpv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Scales every lane in `x` by `2^exponent` with the binary16 scalar kernel,
/// broadcasting the `i32` `exponent` to every lane.
pub fn ldexphv(x: anytype, exponent: i32) @TypeOf(x) {
    return core.mapUnaryExtra(f16, f16, i32, scalar.ldexph, exponent, x);
}

/// Scales every lane in `x` by `2^exponent` with the binary32 scalar kernel,
/// broadcasting the `i32` `exponent` to every lane.
pub fn ldexpfv(x: anytype, exponent: i32) @TypeOf(x) {
    return core.mapUnaryExtra(f32, f32, i32, scalar.ldexpf, exponent, x);
}

/// Scales every lane in `x` by `2^exponent` with the binary64 scalar kernel,
/// broadcasting the `i32` `exponent` to every lane.
pub fn ldexpdv(x: anytype, exponent: i32) @TypeOf(x) {
    return core.mapUnaryExtra(f64, f64, i32, scalar.ldexpd, exponent, x);
}

/// Scales every lane in `x` by its own per-lane `i32` exponent with the scalar
/// kernel for the lane type of `x`.
pub fn ldexpvLanes(x: anytype, e: anytype) @TypeOf(x) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return core.mapIntLanes(f16, f16, i32, scalar.ldexph, x, e);
    if (lane == f32) return core.mapIntLanes(f32, f32, i32, scalar.ldexpf, x, e);
    if (lane == f64) return core.mapIntLanes(f64, f64, i32, scalar.ldexpd, x, e);
    @compileError("zmath.ldexpvLanes supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Scales every lane in `x` by its own per-lane `i32` exponent with the
/// binary16 scalar kernel.
pub fn ldexphvLanes(x: anytype, e: anytype) @TypeOf(x) {
    return core.mapIntLanes(f16, f16, i32, scalar.ldexph, x, e);
}

/// Scales every lane in `x` by its own per-lane `i32` exponent with the
/// binary32 scalar kernel.
pub fn ldexpfvLanes(x: anytype, e: anytype) @TypeOf(x) {
    return core.mapIntLanes(f32, f32, i32, scalar.ldexpf, x, e);
}

/// Scales every lane in `x` by its own per-lane `i32` exponent with the
/// binary64 scalar kernel.
pub fn ldexpdvLanes(x: anytype, e: anytype) @TypeOf(x) {
    return core.mapIntLanes(f64, f64, i32, scalar.ldexpd, x, e);
}

// ---------------------------------------------------------------------------
// scalbn
// ---------------------------------------------------------------------------

/// Scales every lane in `x` by an exact power of two with the scalar kernel
/// for its lane type, broadcasting the `i32` `exponent` to every lane.
pub fn scalbnv(x: anytype, exponent: i32) @TypeOf(x) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return core.mapUnaryExtra(f16, f16, i32, scalar.scalbnh, exponent, x);
    if (lane == f32) return core.mapUnaryExtra(f32, f32, i32, scalar.scalbnf, exponent, x);
    if (lane == f64) return core.mapUnaryExtra(f64, f64, i32, scalar.scalbnd, exponent, x);
    @compileError("zmath.scalbnv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Scales every lane in `x` by an exact power of two with the binary16 scalar
/// kernel, broadcasting the `i32` `exponent` to every lane.
pub fn scalbnhv(x: anytype, exponent: i32) @TypeOf(x) {
    return core.mapUnaryExtra(f16, f16, i32, scalar.scalbnh, exponent, x);
}

/// Scales every lane in `x` by an exact power of two with the binary32 scalar
/// kernel, broadcasting the `i32` `exponent` to every lane.
pub fn scalbnfv(x: anytype, exponent: i32) @TypeOf(x) {
    return core.mapUnaryExtra(f32, f32, i32, scalar.scalbnf, exponent, x);
}

/// Scales every lane in `x` by an exact power of two with the binary64 scalar
/// kernel, broadcasting the `i32` `exponent` to every lane.
pub fn scalbndv(x: anytype, exponent: i32) @TypeOf(x) {
    return core.mapUnaryExtra(f64, f64, i32, scalar.scalbnd, exponent, x);
}

/// Scales every lane in `x` by its own per-lane `i32` exponent with the scalar
/// kernel for the lane type of `x`.
pub fn scalbnvLanes(x: anytype, e: anytype) @TypeOf(x) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return core.mapIntLanes(f16, f16, i32, scalar.scalbnh, x, e);
    if (lane == f32) return core.mapIntLanes(f32, f32, i32, scalar.scalbnf, x, e);
    if (lane == f64) return core.mapIntLanes(f64, f64, i32, scalar.scalbnd, x, e);
    @compileError("zmath.scalbnvLanes supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Scales every lane in `x` by its own per-lane `i32` exponent with the
/// binary16 scalar kernel.
pub fn scalbnhvLanes(x: anytype, e: anytype) @TypeOf(x) {
    return core.mapIntLanes(f16, f16, i32, scalar.scalbnh, x, e);
}

/// Scales every lane in `x` by its own per-lane `i32` exponent with the
/// binary32 scalar kernel.
pub fn scalbnfvLanes(x: anytype, e: anytype) @TypeOf(x) {
    return core.mapIntLanes(f32, f32, i32, scalar.scalbnf, x, e);
}

/// Scales every lane in `x` by its own per-lane `i32` exponent with the
/// binary64 scalar kernel.
pub fn scalbndvLanes(x: anytype, e: anytype) @TypeOf(x) {
    return core.mapIntLanes(f64, f64, i32, scalar.scalbnd, x, e);
}

// ---------------------------------------------------------------------------
// scalbln
// ---------------------------------------------------------------------------

/// Scales every lane in `x` by an exact long power of two with the scalar
/// kernel for its lane type, broadcasting the `i64` `exponent` to every lane.
pub fn scalblnv(x: anytype, exponent: i64) @TypeOf(x) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return core.mapUnaryExtra(f16, f16, i64, scalar.scalblnh, exponent, x);
    if (lane == f32) return core.mapUnaryExtra(f32, f32, i64, scalar.scalblnf, exponent, x);
    if (lane == f64) return core.mapUnaryExtra(f64, f64, i64, scalar.scalblnd, exponent, x);
    @compileError("zmath.scalblnv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Scales every lane in `x` by an exact long power of two with the binary16
/// scalar kernel, broadcasting the `i64` `exponent` to every lane.
pub fn scalblnhv(x: anytype, exponent: i64) @TypeOf(x) {
    return core.mapUnaryExtra(f16, f16, i64, scalar.scalblnh, exponent, x);
}

/// Scales every lane in `x` by an exact long power of two with the binary32
/// scalar kernel, broadcasting the `i64` `exponent` to every lane.
pub fn scalblnfv(x: anytype, exponent: i64) @TypeOf(x) {
    return core.mapUnaryExtra(f32, f32, i64, scalar.scalblnf, exponent, x);
}

/// Scales every lane in `x` by an exact long power of two with the binary64
/// scalar kernel, broadcasting the `i64` `exponent` to every lane.
pub fn scalblndv(x: anytype, exponent: i64) @TypeOf(x) {
    return core.mapUnaryExtra(f64, f64, i64, scalar.scalblnd, exponent, x);
}

/// Scales every lane in `x` by its own per-lane `i64` exponent with the scalar
/// kernel for the lane type of `x`.
pub fn scalblnvLanes(x: anytype, e: anytype) @TypeOf(x) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return core.mapIntLanes(f16, f16, i64, scalar.scalblnh, x, e);
    if (lane == f32) return core.mapIntLanes(f32, f32, i64, scalar.scalblnf, x, e);
    if (lane == f64) return core.mapIntLanes(f64, f64, i64, scalar.scalblnd, x, e);
    @compileError("zmath.scalblnvLanes supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Scales every lane in `x` by its own per-lane `i64` exponent with the
/// binary16 scalar kernel.
pub fn scalblnhvLanes(x: anytype, e: anytype) @TypeOf(x) {
    return core.mapIntLanes(f16, f16, i64, scalar.scalblnh, x, e);
}

/// Scales every lane in `x` by its own per-lane `i64` exponent with the
/// binary32 scalar kernel.
pub fn scalblnfvLanes(x: anytype, e: anytype) @TypeOf(x) {
    return core.mapIntLanes(f32, f32, i64, scalar.scalblnf, x, e);
}

/// Scales every lane in `x` by its own per-lane `i64` exponent with the
/// binary64 scalar kernel.
pub fn scalblndvLanes(x: anytype, e: anytype) @TypeOf(x) {
    return core.mapIntLanes(f64, f64, i64, scalar.scalblnd, x, e);
}

// ---------------------------------------------------------------------------
// pown
// ---------------------------------------------------------------------------

/// Raises every lane in `x` to the `i64` `exponent` with the scalar kernel for
/// its lane type, broadcasting `exponent` to every lane.
///
/// The result is a lane-ordered `[N]PowerResult(T)`, not a `@Vector`: Zig
/// forbids a struct as a `@Vector` element type.
pub fn pownv(x: anytype, exponent: i64) PowerVector(LaneOf(@TypeOf(x)), @TypeOf(x)) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return mapPowerBroadcast(f16, i64, scalar.pownh, exponent, x);
    if (lane == f32) return mapPowerBroadcast(f32, i64, scalar.pownf, exponent, x);
    if (lane == f64) return mapPowerBroadcast(f64, i64, scalar.pownd, exponent, x);
    @compileError("zmath.pownv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Raises every lane in `x` to the `i64` `exponent` with the binary16 scalar
/// kernel, broadcasting `exponent` to every lane.
///
/// The result is a lane-ordered `[N]PowerResult(f16)`, not a `@Vector`.
pub fn pownhv(x: anytype, exponent: i64) PowerVector(f16, @TypeOf(x)) {
    return mapPowerBroadcast(f16, i64, scalar.pownh, exponent, x);
}

/// Raises every lane in `x` to the `i64` `exponent` with the binary32 scalar
/// kernel, broadcasting `exponent` to every lane.
///
/// The result is a lane-ordered `[N]PowerResult(f32)`, not a `@Vector`.
pub fn pownfv(x: anytype, exponent: i64) PowerVector(f32, @TypeOf(x)) {
    return mapPowerBroadcast(f32, i64, scalar.pownf, exponent, x);
}

/// Raises every lane in `x` to the `i64` `exponent` with the binary64 scalar
/// kernel, broadcasting `exponent` to every lane.
///
/// The result is a lane-ordered `[N]PowerResult(f64)`, not a `@Vector`.
pub fn powndv(x: anytype, exponent: i64) PowerVector(f64, @TypeOf(x)) {
    return mapPowerBroadcast(f64, i64, scalar.pownd, exponent, x);
}

/// Raises every lane in `x` to its own per-lane `i64` exponent with the scalar
/// kernel for the lane type of `x`.
///
/// The result is a lane-ordered `[N]PowerResult(T)`, not a `@Vector`.
pub fn pownvLanes(x: anytype, e: anytype) PowerVector(LaneOf(@TypeOf(x)), @TypeOf(x)) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return mapPowerLanes(f16, i64, scalar.pownh, x, e);
    if (lane == f32) return mapPowerLanes(f32, i64, scalar.pownf, x, e);
    if (lane == f64) return mapPowerLanes(f64, i64, scalar.pownd, x, e);
    @compileError("zmath.pownvLanes supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Raises every lane in `x` to its own per-lane `i64` exponent with the
/// binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]PowerResult(f16)`, not a `@Vector`.
pub fn pownhvLanes(x: anytype, e: anytype) PowerVector(f16, @TypeOf(x)) {
    return mapPowerLanes(f16, i64, scalar.pownh, x, e);
}

/// Raises every lane in `x` to its own per-lane `i64` exponent with the
/// binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]PowerResult(f32)`, not a `@Vector`.
pub fn pownfvLanes(x: anytype, e: anytype) PowerVector(f32, @TypeOf(x)) {
    return mapPowerLanes(f32, i64, scalar.pownf, x, e);
}

/// Raises every lane in `x` to its own per-lane `i64` exponent with the
/// binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]PowerResult(f64)`, not a `@Vector`.
pub fn powndvLanes(x: anytype, e: anytype) PowerVector(f64, @TypeOf(x)) {
    return mapPowerLanes(f64, i64, scalar.pownd, x, e);
}

// ---------------------------------------------------------------------------
// compoundn
// ---------------------------------------------------------------------------

/// Raises one plus every lane in `x` to the `i64` `exponent` with the scalar
/// kernel for its lane type, broadcasting `exponent` to every lane.
///
/// The result is a lane-ordered `[N]PowerResult(T)`, not a `@Vector`: Zig
/// forbids a struct as a `@Vector` element type.
pub fn compoundnv(x: anytype, exponent: i64) PowerVector(LaneOf(@TypeOf(x)), @TypeOf(x)) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return mapPowerBroadcast(f16, i64, scalar.compoundnh, exponent, x);
    if (lane == f32) return mapPowerBroadcast(f32, i64, scalar.compoundnf, exponent, x);
    if (lane == f64) return mapPowerBroadcast(f64, i64, scalar.compoundnd, exponent, x);
    @compileError("zmath.compoundnv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Raises one plus every lane in `x` to the `i64` `exponent` with the binary16
/// scalar kernel, broadcasting `exponent` to every lane.
///
/// The result is a lane-ordered `[N]PowerResult(f16)`, not a `@Vector`.
pub fn compoundnhv(x: anytype, exponent: i64) PowerVector(f16, @TypeOf(x)) {
    return mapPowerBroadcast(f16, i64, scalar.compoundnh, exponent, x);
}

/// Raises one plus every lane in `x` to the `i64` `exponent` with the binary32
/// scalar kernel, broadcasting `exponent` to every lane.
///
/// The result is a lane-ordered `[N]PowerResult(f32)`, not a `@Vector`.
pub fn compoundnfv(x: anytype, exponent: i64) PowerVector(f32, @TypeOf(x)) {
    return mapPowerBroadcast(f32, i64, scalar.compoundnf, exponent, x);
}

/// Raises one plus every lane in `x` to the `i64` `exponent` with the binary64
/// scalar kernel, broadcasting `exponent` to every lane.
///
/// The result is a lane-ordered `[N]PowerResult(f64)`, not a `@Vector`.
pub fn compoundndv(x: anytype, exponent: i64) PowerVector(f64, @TypeOf(x)) {
    return mapPowerBroadcast(f64, i64, scalar.compoundnd, exponent, x);
}

/// Raises one plus every lane in `x` to its own per-lane `i64` exponent with
/// the scalar kernel for the lane type of `x`.
///
/// The result is a lane-ordered `[N]PowerResult(T)`, not a `@Vector`.
pub fn compoundnvLanes(x: anytype, e: anytype) PowerVector(LaneOf(@TypeOf(x)), @TypeOf(x)) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return mapPowerLanes(f16, i64, scalar.compoundnh, x, e);
    if (lane == f32) return mapPowerLanes(f32, i64, scalar.compoundnf, x, e);
    if (lane == f64) return mapPowerLanes(f64, i64, scalar.compoundnd, x, e);
    @compileError("zmath.compoundnvLanes supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Raises one plus every lane in `x` to its own per-lane `i64` exponent with
/// the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]PowerResult(f16)`, not a `@Vector`.
pub fn compoundnhvLanes(x: anytype, e: anytype) PowerVector(f16, @TypeOf(x)) {
    return mapPowerLanes(f16, i64, scalar.compoundnh, x, e);
}

/// Raises one plus every lane in `x` to its own per-lane `i64` exponent with
/// the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]PowerResult(f32)`, not a `@Vector`.
pub fn compoundnfvLanes(x: anytype, e: anytype) PowerVector(f32, @TypeOf(x)) {
    return mapPowerLanes(f32, i64, scalar.compoundnf, x, e);
}

/// Raises one plus every lane in `x` to its own per-lane `i64` exponent with
/// the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]PowerResult(f64)`, not a `@Vector`.
pub fn compoundndvLanes(x: anytype, e: anytype) PowerVector(f64, @TypeOf(x)) {
    return mapPowerLanes(f64, i64, scalar.compoundnd, x, e);
}

// ---------------------------------------------------------------------------
// rootn
// ---------------------------------------------------------------------------

/// Computes the `exponent`-th root of every lane in `x` with the scalar kernel
/// for its lane type, broadcasting the `i64` `exponent` to every lane.
///
/// The result is a lane-ordered `[N]PowerResult(T)`, not a `@Vector`: Zig
/// forbids a struct as a `@Vector` element type.
pub fn rootnv(x: anytype, exponent: i64) PowerVector(LaneOf(@TypeOf(x)), @TypeOf(x)) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return mapPowerBroadcast(f16, i64, scalar.rootnh, exponent, x);
    if (lane == f32) return mapPowerBroadcast(f32, i64, scalar.rootnf, exponent, x);
    if (lane == f64) return mapPowerBroadcast(f64, i64, scalar.rootnd, exponent, x);
    @compileError("zmath.rootnv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the `exponent`-th root of every lane in `x` with the binary16
/// scalar kernel, broadcasting the `i64` `exponent` to every lane.
///
/// The result is a lane-ordered `[N]PowerResult(f16)`, not a `@Vector`.
pub fn rootnhv(x: anytype, exponent: i64) PowerVector(f16, @TypeOf(x)) {
    return mapPowerBroadcast(f16, i64, scalar.rootnh, exponent, x);
}

/// Computes the `exponent`-th root of every lane in `x` with the binary32
/// scalar kernel, broadcasting the `i64` `exponent` to every lane.
///
/// The result is a lane-ordered `[N]PowerResult(f32)`, not a `@Vector`.
pub fn rootnfv(x: anytype, exponent: i64) PowerVector(f32, @TypeOf(x)) {
    return mapPowerBroadcast(f32, i64, scalar.rootnf, exponent, x);
}

/// Computes the `exponent`-th root of every lane in `x` with the binary64
/// scalar kernel, broadcasting the `i64` `exponent` to every lane.
///
/// The result is a lane-ordered `[N]PowerResult(f64)`, not a `@Vector`.
pub fn rootndv(x: anytype, exponent: i64) PowerVector(f64, @TypeOf(x)) {
    return mapPowerBroadcast(f64, i64, scalar.rootnd, exponent, x);
}

/// Computes the `exponent`-th root of every lane in `x` with the scalar kernel
/// for the lane type of `x`, using each lane's own per-lane `i64` exponent.
///
/// The result is a lane-ordered `[N]PowerResult(T)`, not a `@Vector`.
pub fn rootnvLanes(x: anytype, e: anytype) PowerVector(LaneOf(@TypeOf(x)), @TypeOf(x)) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return mapPowerLanes(f16, i64, scalar.rootnh, x, e);
    if (lane == f32) return mapPowerLanes(f32, i64, scalar.rootnf, x, e);
    if (lane == f64) return mapPowerLanes(f64, i64, scalar.rootnd, x, e);
    @compileError("zmath.rootnvLanes supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes each lane's own per-lane `i64`-th root with the binary16 scalar
/// kernel.
///
/// The result is a lane-ordered `[N]PowerResult(f16)`, not a `@Vector`.
pub fn rootnhvLanes(x: anytype, e: anytype) PowerVector(f16, @TypeOf(x)) {
    return mapPowerLanes(f16, i64, scalar.rootnh, x, e);
}

/// Computes each lane's own per-lane `i64`-th root with the binary32 scalar
/// kernel.
///
/// The result is a lane-ordered `[N]PowerResult(f32)`, not a `@Vector`.
pub fn rootnfvLanes(x: anytype, e: anytype) PowerVector(f32, @TypeOf(x)) {
    return mapPowerLanes(f32, i64, scalar.rootnf, x, e);
}

/// Computes each lane's own per-lane `i64`-th root with the binary64 scalar
/// kernel.
///
/// The result is a lane-ordered `[N]PowerResult(f64)`, not a `@Vector`.
pub fn rootndvLanes(x: anytype, e: anytype) PowerVector(f64, @TypeOf(x)) {
    return mapPowerLanes(f64, i64, scalar.rootnd, x, e);
}

// ---------------------------------------------------------------------------
// nextToward
// ---------------------------------------------------------------------------

/// Returns the next representable value of each lane in `x` toward the `f64`
/// `toward` with the scalar kernel for its lane type.
///
/// `toward` is always `f64`, never the lane type, so there is no per-lane form.
pub fn nextTowardv(x: anytype, toward: f64) @TypeOf(x) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return core.mapUnaryExtra(f16, f16, f64, ieee.nextTowardh, toward, x);
    if (lane == f32) return core.mapUnaryExtra(f32, f32, f64, ieee.nextTowardf, toward, x);
    if (lane == f64) return core.mapUnaryExtra(f64, f64, f64, ieee.nextTowardd, toward, x);
    @compileError("zmath.nextTowardv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Returns the next representable value of each lane in `x` toward the `f64`
/// `toward` with the binary16 scalar kernel.
pub fn nextTowardhv(x: anytype, toward: f64) @TypeOf(x) {
    return core.mapUnaryExtra(f16, f16, f64, ieee.nextTowardh, toward, x);
}

/// Returns the next representable value of each lane in `x` toward the `f64`
/// `toward` with the binary32 scalar kernel.
pub fn nextTowardfv(x: anytype, toward: f64) @TypeOf(x) {
    return core.mapUnaryExtra(f32, f32, f64, ieee.nextTowardf, toward, x);
}

/// Returns the next representable value of each lane in `x` toward the `f64`
/// `toward` with the binary64 scalar kernel.
pub fn nextTowarddv(x: anytype, toward: f64) @TypeOf(x) {
    return core.mapUnaryExtra(f64, f64, f64, ieee.nextTowardd, toward, x);
}

// ---------------------------------------------------------------------------
// fromfp
//
// `fromfpx` and `ufromfpx` are the C23 inexact-reporting spellings. They keep
// their own public names to mirror the root surface, but each kernel body is
// byte-identical to the corresponding non-`x` kernel in `scalar_fromfp.zig`.
// ---------------------------------------------------------------------------

/// Converts every lane in `x` to an integral value with the scalar kernel for
/// its lane type, broadcasting `round` and `width` to every lane.
pub fn fromfpv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return core.mapUnaryExtra2(
        f16,
        f16,
        types.FromFpRound,
        u32,
        scalar.fromfph,
        round,
        width,
        x,
    );
    if (lane == f32) return core.mapUnaryExtra2(
        f32,
        f32,
        types.FromFpRound,
        u32,
        scalar.fromfpf,
        round,
        width,
        x,
    );
    if (lane == f64) return core.mapUnaryExtra2(
        f64,
        f64,
        types.FromFpRound,
        u32,
        scalar.fromfpd,
        round,
        width,
        x,
    );
    @compileError("zmath.fromfpv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Converts every lane in `x` to an integral value with the binary16 scalar
/// kernel, broadcasting `round` and `width` to every lane.
pub fn fromfphv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    return core.mapUnaryExtra2(f16, f16, types.FromFpRound, u32, scalar.fromfph, round, width, x);
}

/// Converts every lane in `x` to an integral value with the binary32 scalar
/// kernel, broadcasting `round` and `width` to every lane.
pub fn fromfpfv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    return core.mapUnaryExtra2(f32, f32, types.FromFpRound, u32, scalar.fromfpf, round, width, x);
}

/// Converts every lane in `x` to an integral value with the binary64 scalar
/// kernel, broadcasting `round` and `width` to every lane.
pub fn fromfpdv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    return core.mapUnaryExtra2(f64, f64, types.FromFpRound, u32, scalar.fromfpd, round, width, x);
}

/// Converts every lane in `x` to an unsigned integral value with the scalar
/// kernel for its lane type, broadcasting `round` and `width` to every lane.
pub fn ufromfpv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return core.mapUnaryExtra2(
        f16,
        f16,
        types.FromFpRound,
        u32,
        scalar.ufromfph,
        round,
        width,
        x,
    );
    if (lane == f32) return core.mapUnaryExtra2(
        f32,
        f32,
        types.FromFpRound,
        u32,
        scalar.ufromfpf,
        round,
        width,
        x,
    );
    if (lane == f64) return core.mapUnaryExtra2(
        f64,
        f64,
        types.FromFpRound,
        u32,
        scalar.ufromfpd,
        round,
        width,
        x,
    );
    @compileError("zmath.ufromfpv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Converts every lane in `x` to an unsigned integral value with the binary16
/// scalar kernel, broadcasting `round` and `width` to every lane.
pub fn ufromfphv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    return core.mapUnaryExtra2(f16, f16, types.FromFpRound, u32, scalar.ufromfph, round, width, x);
}

/// Converts every lane in `x` to an unsigned integral value with the binary32
/// scalar kernel, broadcasting `round` and `width` to every lane.
pub fn ufromfpfv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    return core.mapUnaryExtra2(f32, f32, types.FromFpRound, u32, scalar.ufromfpf, round, width, x);
}

/// Converts every lane in `x` to an unsigned integral value with the binary64
/// scalar kernel, broadcasting `round` and `width` to every lane.
pub fn ufromfpdv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    return core.mapUnaryExtra2(f64, f64, types.FromFpRound, u32, scalar.ufromfpd, round, width, x);
}

/// Converts every lane in `x` to an integral value with the C23
/// inexact-reporting spelling, broadcasting `round` and `width` to every lane.
///
/// Byte-identical to `fromfpv`; only the kernel spelling differs.
pub fn fromfpxv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return core.mapUnaryExtra2(
        f16,
        f16,
        types.FromFpRound,
        u32,
        scalar.fromfpxh,
        round,
        width,
        x,
    );
    if (lane == f32) return core.mapUnaryExtra2(
        f32,
        f32,
        types.FromFpRound,
        u32,
        scalar.fromfpxf,
        round,
        width,
        x,
    );
    if (lane == f64) return core.mapUnaryExtra2(
        f64,
        f64,
        types.FromFpRound,
        u32,
        scalar.fromfpxd,
        round,
        width,
        x,
    );
    @compileError("zmath.fromfpxv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Converts every lane in `x` to an integral value with the C23
/// inexact-reporting spelling and the binary16 scalar kernel.
pub fn fromfpxhv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    return core.mapUnaryExtra2(f16, f16, types.FromFpRound, u32, scalar.fromfpxh, round, width, x);
}

/// Converts every lane in `x` to an integral value with the C23
/// inexact-reporting spelling and the binary32 scalar kernel.
pub fn fromfpxfv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    return core.mapUnaryExtra2(f32, f32, types.FromFpRound, u32, scalar.fromfpxf, round, width, x);
}

/// Converts every lane in `x` to an integral value with the C23
/// inexact-reporting spelling and the binary64 scalar kernel.
pub fn fromfpxdv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    return core.mapUnaryExtra2(f64, f64, types.FromFpRound, u32, scalar.fromfpxd, round, width, x);
}

/// Converts every lane in `x` to an unsigned integral value with the C23
/// inexact-reporting spelling, broadcasting `round` and `width` to every lane.
///
/// Byte-identical to `ufromfpv`; only the kernel spelling differs.
pub fn ufromfpxv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    const lane = LaneOf(@TypeOf(x));
    if (lane == f16) return core.mapUnaryExtra2(
        f16,
        f16,
        types.FromFpRound,
        u32,
        scalar.ufromfpxh,
        round,
        width,
        x,
    );
    if (lane == f32) return core.mapUnaryExtra2(
        f32,
        f32,
        types.FromFpRound,
        u32,
        scalar.ufromfpxf,
        round,
        width,
        x,
    );
    if (lane == f64) return core.mapUnaryExtra2(
        f64,
        f64,
        types.FromFpRound,
        u32,
        scalar.ufromfpxd,
        round,
        width,
        x,
    );
    @compileError("zmath.ufromfpxv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Converts every lane in `x` to an unsigned integral value with the C23
/// inexact-reporting spelling and the binary16 scalar kernel.
pub fn ufromfpxhv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    return core.mapUnaryExtra2(f16, f16, types.FromFpRound, u32, scalar.ufromfpxh, round, width, x);
}

/// Converts every lane in `x` to an unsigned integral value with the C23
/// inexact-reporting spelling and the binary32 scalar kernel.
pub fn ufromfpxfv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    return core.mapUnaryExtra2(f32, f32, types.FromFpRound, u32, scalar.ufromfpxf, round, width, x);
}

/// Converts every lane in `x` to an unsigned integral value with the C23
/// inexact-reporting spelling and the binary64 scalar kernel.
pub fn ufromfpxdv(x: anytype, round: types.FromFpRound, width: u32) @TypeOf(x) {
    return core.mapUnaryExtra2(f64, f64, types.FromFpRound, u32, scalar.ufromfpxd, round, width, x);
}

// ---------------------------------------------------------------------------
// Self check: every entry point above is referenced, so none is dead code.
// ---------------------------------------------------------------------------

comptime {
    _ = widthOf;
    _ = LaneOf;
    _ = PowerVector;
    _ = mapPowerBroadcast;
    _ = mapPowerLanes;

    _ = ldexpv;
    _ = ldexphv;
    _ = ldexpfv;
    _ = ldexpdv;
    _ = ldexpvLanes;
    _ = ldexphvLanes;
    _ = ldexpfvLanes;
    _ = ldexpdvLanes;

    _ = scalbnv;
    _ = scalbnhv;
    _ = scalbnfv;
    _ = scalbndv;
    _ = scalbnvLanes;
    _ = scalbnhvLanes;
    _ = scalbnfvLanes;
    _ = scalbndvLanes;

    _ = scalblnv;
    _ = scalblnhv;
    _ = scalblnfv;
    _ = scalblndv;
    _ = scalblnvLanes;
    _ = scalblnhvLanes;
    _ = scalblnfvLanes;
    _ = scalblndvLanes;

    _ = pownv;
    _ = pownhv;
    _ = pownfv;
    _ = powndv;
    _ = pownvLanes;
    _ = pownhvLanes;
    _ = pownfvLanes;
    _ = powndvLanes;

    _ = compoundnv;
    _ = compoundnhv;
    _ = compoundnfv;
    _ = compoundndv;
    _ = compoundnvLanes;
    _ = compoundnhvLanes;
    _ = compoundnfvLanes;
    _ = compoundndvLanes;

    _ = rootnv;
    _ = rootnhv;
    _ = rootnfv;
    _ = rootndv;
    _ = rootnvLanes;
    _ = rootnhvLanes;
    _ = rootnfvLanes;
    _ = rootndvLanes;

    _ = nextTowardv;
    _ = nextTowardhv;
    _ = nextTowardfv;
    _ = nextTowarddv;

    _ = fromfpv;
    _ = fromfphv;
    _ = fromfpfv;
    _ = fromfpdv;

    _ = ufromfpv;
    _ = ufromfphv;
    _ = ufromfpfv;
    _ = ufromfpdv;

    _ = fromfpxv;
    _ = fromfpxhv;
    _ = fromfpxfv;
    _ = fromfpxdv;

    _ = ufromfpxv;
    _ = ufromfpxhv;
    _ = ufromfpxfv;
    _ = ufromfpxdv;
}
