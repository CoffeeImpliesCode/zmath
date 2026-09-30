//! Vector (`@Vector`) overloads for the elementwise scalar operations.
//!
//! Every operation here applies the *identical* accepted scalar kernel to each
//! lane, so lane `i` of the result is bit-for-bit equal to the scalar result for
//! lane `i` of the input. This preserves the schema-3 `max_ulp = 0` contract and
//! the immutable corpus without new evidence: the vector path is a per-lane
//! application of the same scalar code, never a re-derived polynomial.
//!
//! Three naming forms are provided for each elementwise float operation:
//!
//! - `zmath.sin(x)` where `x: @Vector(N, T)` routes through the generic top-level
//!   dispatch in `root.zig`.
//! - The generic vector entry point `sinv(x: @Vector(N, T))`, which selects the
//!   scalar kernel from the lane type `T`.
//! - The specialized entry points `sinfv` and `sindv`, which pin the lane type to
//!   `f32` and `f64` respectively, plus `sinhv` for binary16 where that name is
//!   not already used by another generic entry point (for example `sin` has no
//!   `sinhv`, because `sinh` already defines the generic `sinhv`).
//!
//! ## Families
//!
//! This module holds the float-to-float surface, where a lane maps to a lane of
//! the same float type: `@Vector(N, T) -> @Vector(N, T)`. The operations whose
//! lane maps to something else live in four sibling modules, re-exported here so
//! the whole `@Vector` surface is reachable as `zmath.vector.*`.
//!
//! - `narrow` holds the ternary float-to-float operation `fma`, and the six C23
//!   narrowing operations `fadd`, `fsub`, `fmul`, `fdiv`, `ffma` and `fsqrt`,
//!   which take `f64` lanes and produce `f32` lanes. Their lane type and their
//!   result element type are both fixed, so they carry no per-format forms.
//! - `ieee_family` holds the remaining IEEE and miscellaneous float-to-float
//!   operations, the integer-returning `ilogb` and `llogb`, and the
//!   boolean-returning `totalOrder` and `totalOrderMag`.
//! - `intarg` holds the eleven operations that take a non-float argument: the
//!   `i32` exponent of `ldexp` and `scalbn`, the `i64` exponent of `scalbln`
//!   and of the `pown`, `compoundn` and `rootn` power family, the `f64` target
//!   of `nextToward`, and the `FromFpRound` direction and `u32` width of the
//!   `fromfp` family. Each extra argument is a broadcast runtime scalar; the six
//!   exponent operations also take a per-lane integer vector under a `Lanes`
//!   suffix.
//! - `struct_family` holds the thirteen operations whose scalar kernel returns a
//!   result struct.
//! - `core` is the shared per-lane mapping core the other four build on. It is
//!   public so a caller can write a new lane-parallel entry point from the same
//!   helpers, with the same per-lane bit-identity.
//!
//! ## Result shape
//!
//! Most results are a real `@Vector(N, R)`: the lane type for float-to-float,
//! `f32` for the narrowing operations, `i32` or `i64` for `ilogb` and `llogb`,
//! and `bool` for `totalOrder` and `totalOrderMag`.
//!
//! The result-struct operations cannot be. Zig forbids a struct as a `@Vector`
//! element type and reports "expected integer, float, bool, or pointer for the
//! vector element type", so `struct_family` returns a lane-ordered array
//! `[N]Result` instead, and so do `pown`, `compoundn` and `rootn` in `intarg`.
//! Element `i` of that array is the value the scalar kernel itself returned for
//! lane `i`, moved as a whole, in lane order, so the per-lane bit-identity
//! guarantee is unchanged. Note what this shape does NOT buy: the kernel runs
//! once per lane, exactly as a scalar loop would, so there is no throughput
//! advantage and no SIMD acceleration over a scalar loop. The value is a
//! uniform lane-parallel API, nothing more.

const std = @import("std");
const scalar = @import("scalar.zig");
const value = @import("vector_value.zig");

/// The shared per-lane mapping core, behind `narrow`, `ieee_family` and
/// `intarg`.
pub const core = @import("vector_core.zig");

/// The thirteen operations whose scalar kernel returns a result struct.
pub const struct_family = @import("vector_struct.zig");

/// The eleven operations that take a non-float extra argument.
pub const intarg = @import("vector_intarg.zig");

/// The ternary `fma` and the six C23 narrowing operations.
pub const narrow = @import("vector_narrow.zig");

/// The remaining IEEE and miscellaneous operations, including the
/// integer-returning and boolean-returning ones.
pub const ieee_family = @import("vector_ieee.zig");

/// The named-component vector value types. This is the value-type layer named
/// by `API.md`, not the `@Vector` layer above; the two share this module path
/// and nothing else. Every `@Vector` entry point here ends in `v`, and no
/// `Vec*` name does, so the surfaces stay unambiguous at the call site.
pub const Vec2f = value.Vec2f;
pub const Vec3f = value.Vec3f;
pub const Vec4f = value.Vec4f;
pub const Vec2d = value.Vec2d;
pub const Vec3d = value.Vec3d;
pub const Vec4d = value.Vec4d;

// ---------------------------------------------------------------------------
// Core per-lane application
// ---------------------------------------------------------------------------

// Names follow `op + type + v`: `expfv` applies the binary32 kernel, `expdv`
// the binary64 one, and `expv` picks the kernel from the lane type. Where an
// operation name already ends in a format letter the two forms collide, because
// `sinhv`, `coshv`, `tanhv`, `asinhv`, `acoshv`, and `atanhv` are the GENERIC
// entry points for `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, and `atanh`. The
// binary16 specializations for `sin`, `cos`, `tan`, `asin`, `acos`, and `atan`
// are therefore intentionally absent, and binary16 callers use the generic
// `sinv`, `cosv`, `tanv`, `asinv`, `acosv`, and `atanv`.
//
// The binary16 kernel of a bare operation is not its hyperbolic twin: `root.sin`
// binds f16 to `scalar.sinh`, the binary16 SINE, and not `scalar.sinhh`, the
// binary16 HYPERBOLIC SINE. The same holds for `cos` -> `scalar.cosh`,
// `tan` -> `scalar.tanh`, and the inverse-hyperbolic equivalents. The surface
// below is machine-derived from the dispatch lines in `root.zig`, and
// `vector_surface_test.zig` is the guard over that derivation.
//
// The two helpers in this section carry a lane to a lane of the same float
// type, so they serve only the float-to-float surface declared above. An
// integer, boolean, or narrowing result needs the explicit result element type
// of `core`, and a result struct needs the lane-ordered array helpers of
// `struct_family`.

/// Applies the unary scalar kernel `f` to every lane of `x`.
///
/// The result is bit-identical to calling `f` on each lane, because each lane
/// is computed by the same scalar function on the same value. The vector width
/// is inferred from the argument, so this works for every `@Vector(N, T)`.
fn mapUnary(comptime T: type, comptime f: fn (T) T, x: anytype) @TypeOf(x) {
    const info = @typeInfo(@TypeOf(x));
    const n = switch (info) {
        .vector => |v| v.len,
        else => @compileError("zmath vector operations require a @Vector argument"),
    };
    if (info.vector.child != T) {
        @compileError("zmath vector lane type does not match the selected scalar kernel");
    }
    var out: @TypeOf(x) = undefined;
    inline for (0..n) |i| out[i] = f(x[i]);
    return out;
}

/// Applies the binary scalar kernel `f` to corresponding lanes of `a` and `b`.
fn mapBinary(comptime T: type, comptime f: fn (T, T) T, a: anytype, b: anytype) @TypeOf(a) {
    const info = @typeInfo(@TypeOf(a));
    const n = switch (info) {
        .vector => |v| v.len,
        else => @compileError("zmath vector operations require a @Vector argument"),
    };
    if (info.vector.child != T) {
        @compileError("zmath vector lane type does not match the selected scalar kernel");
    }
    if (@TypeOf(b) != @TypeOf(a)) {
        @compileError("zmath vector binary operations require both arguments to have the same type");
    }
    var out: @TypeOf(a) = undefined;
    inline for (0..n) |i| out[i] = f(a[i], b[i]);
    return out;
}

// ---------------------------------------------------------------------------
// Generic and specialized vector entry points
// ---------------------------------------------------------------------------

/// Computes the abs of every lane in `x` with the scalar kernel for its lane type.
pub fn absv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.absh, x);
    if (T == f32) return mapUnary(f32, scalar.absf, x);
    if (T == f64) return mapUnary(f64, scalar.absd, x);
    @compileError("zmath.absv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the abs of every lane in `x` with the binary16 scalar kernel.
pub fn abshv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.absh, x);
}

/// Computes the abs of every lane in `x` with the binary32 scalar kernel.
pub fn absfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.absf, x);
}

/// Computes the abs of every lane in `x` with the binary64 scalar kernel.
pub fn absdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.absd, x);
}

/// Computes the exp of every lane in `x` with the scalar kernel for its lane type.
pub fn expv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.exph, x);
    if (T == f32) return mapUnary(f32, scalar.expf, x);
    if (T == f64) return mapUnary(f64, scalar.expd, x);
    @compileError("zmath.expv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the exp of every lane in `x` with the binary16 scalar kernel.
pub fn exphv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.exph, x);
}

/// Computes the exp of every lane in `x` with the binary32 scalar kernel.
pub fn expfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.expf, x);
}

/// Computes the exp of every lane in `x` with the binary64 scalar kernel.
pub fn expdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.expd, x);
}

/// Computes the exp10 of every lane in `x` with the scalar kernel for its lane type.
pub fn exp10v(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.exp10h, x);
    if (T == f32) return mapUnary(f32, scalar.exp10f, x);
    if (T == f64) return mapUnary(f64, scalar.exp10d, x);
    @compileError("zmath.exp10v supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the exp10 of every lane in `x` with the binary16 scalar kernel.
pub fn exp10hv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.exp10h, x);
}

/// Computes the exp10 of every lane in `x` with the binary32 scalar kernel.
pub fn exp10fv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.exp10f, x);
}

/// Computes the exp10 of every lane in `x` with the binary64 scalar kernel.
pub fn exp10dv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.exp10d, x);
}

/// Computes the exp2m1 of every lane in `x` with the scalar kernel for its lane type.
pub fn exp2m1v(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.exp2m1h, x);
    if (T == f32) return mapUnary(f32, scalar.exp2m1f, x);
    if (T == f64) return mapUnary(f64, scalar.exp2m1d, x);
    @compileError("zmath.exp2m1v supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the exp2m1 of every lane in `x` with the binary16 scalar kernel.
pub fn exp2m1hv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.exp2m1h, x);
}

/// Computes the exp2m1 of every lane in `x` with the binary32 scalar kernel.
pub fn exp2m1fv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.exp2m1f, x);
}

/// Computes the exp2m1 of every lane in `x` with the binary64 scalar kernel.
pub fn exp2m1dv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.exp2m1d, x);
}

/// Computes the exp10m1 of every lane in `x` with the scalar kernel for its lane type.
pub fn exp10m1v(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.exp10m1h, x);
    if (T == f32) return mapUnary(f32, scalar.exp10m1f, x);
    if (T == f64) return mapUnary(f64, scalar.exp10m1d, x);
    @compileError("zmath.exp10m1v supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the exp10m1 of every lane in `x` with the binary16 scalar kernel.
pub fn exp10m1hv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.exp10m1h, x);
}

/// Computes the exp10m1 of every lane in `x` with the binary32 scalar kernel.
pub fn exp10m1fv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.exp10m1f, x);
}

/// Computes the exp10m1 of every lane in `x` with the binary64 scalar kernel.
pub fn exp10m1dv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.exp10m1d, x);
}

/// Computes the erfc of every lane in `x` with the scalar kernel for its lane type.
pub fn erfcv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.erfch, x);
    if (T == f32) return mapUnary(f32, scalar.erfcf, x);
    if (T == f64) return mapUnary(f64, scalar.erfcd, x);
    @compileError("zmath.erfcv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the erfc of every lane in `x` with the binary16 scalar kernel.
pub fn erfchv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.erfch, x);
}

/// Computes the erfc of every lane in `x` with the binary32 scalar kernel.
pub fn erfcfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.erfcf, x);
}

/// Computes the erfc of every lane in `x` with the binary64 scalar kernel.
pub fn erfcdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.erfcd, x);
}

/// Computes the tgamma of every lane in `x` with the scalar kernel for its lane type.
pub fn tgammav(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.tgammah, x);
    if (T == f32) return mapUnary(f32, scalar.tgammaf, x);
    if (T == f64) return mapUnary(f64, scalar.tgammad, x);
    @compileError("zmath.tgammav supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the tgamma of every lane in `x` with the binary16 scalar kernel.
pub fn tgammahv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.tgammah, x);
}

/// Computes the tgamma of every lane in `x` with the binary32 scalar kernel.
pub fn tgammafv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.tgammaf, x);
}

/// Computes the tgamma of every lane in `x` with the binary64 scalar kernel.
pub fn tgammadv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.tgammad, x);
}

/// Computes the log of every lane in `x` with the scalar kernel for its lane type.
pub fn logv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.logh, x);
    if (T == f32) return mapUnary(f32, scalar.logf, x);
    if (T == f64) return mapUnary(f64, scalar.logd, x);
    @compileError("zmath.logv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the log of every lane in `x` with the binary16 scalar kernel.
pub fn loghv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.logh, x);
}

/// Computes the log of every lane in `x` with the binary32 scalar kernel.
pub fn logfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.logf, x);
}

/// Computes the log of every lane in `x` with the binary64 scalar kernel.
pub fn logdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.logd, x);
}

/// Computes the log2 of every lane in `x` with the scalar kernel for its lane type.
pub fn log2v(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.log2h, x);
    if (T == f32) return mapUnary(f32, scalar.log2f, x);
    if (T == f64) return mapUnary(f64, scalar.log2d, x);
    @compileError("zmath.log2v supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the log2 of every lane in `x` with the binary16 scalar kernel.
pub fn log2hv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.log2h, x);
}

/// Computes the log2 of every lane in `x` with the binary32 scalar kernel.
pub fn log2fv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.log2f, x);
}

/// Computes the log2 of every lane in `x` with the binary64 scalar kernel.
pub fn log2dv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.log2d, x);
}

/// Computes the log2p1 of every lane in `x` with the scalar kernel for its lane type.
pub fn log2p1v(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.log2p1h, x);
    if (T == f32) return mapUnary(f32, scalar.log2p1f, x);
    if (T == f64) return mapUnary(f64, scalar.log2p1d, x);
    @compileError("zmath.log2p1v supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the log2p1 of every lane in `x` with the binary16 scalar kernel.
pub fn log2p1hv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.log2p1h, x);
}

/// Computes the log2p1 of every lane in `x` with the binary32 scalar kernel.
pub fn log2p1fv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.log2p1f, x);
}

/// Computes the log2p1 of every lane in `x` with the binary64 scalar kernel.
pub fn log2p1dv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.log2p1d, x);
}

/// Computes the log10 of every lane in `x` with the scalar kernel for its lane type.
pub fn log10v(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.log10h, x);
    if (T == f32) return mapUnary(f32, scalar.log10f, x);
    if (T == f64) return mapUnary(f64, scalar.log10d, x);
    @compileError("zmath.log10v supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the log10 of every lane in `x` with the binary16 scalar kernel.
pub fn log10hv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.log10h, x);
}

/// Computes the log10 of every lane in `x` with the binary32 scalar kernel.
pub fn log10fv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.log10f, x);
}

/// Computes the log10 of every lane in `x` with the binary64 scalar kernel.
pub fn log10dv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.log10d, x);
}

/// Computes the log10p1 of every lane in `x` with the scalar kernel for its lane type.
pub fn log10p1v(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.log10p1h, x);
    if (T == f32) return mapUnary(f32, scalar.log10p1f, x);
    if (T == f64) return mapUnary(f64, scalar.log10p1d, x);
    @compileError("zmath.log10p1v supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the log10p1 of every lane in `x` with the binary16 scalar kernel.
pub fn log10p1hv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.log10p1h, x);
}

/// Computes the log10p1 of every lane in `x` with the binary32 scalar kernel.
pub fn log10p1fv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.log10p1f, x);
}

/// Computes the log10p1 of every lane in `x` with the binary64 scalar kernel.
pub fn log10p1dv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.log10p1d, x);
}

/// Computes the erf of every lane in `x` with the scalar kernel for its lane type.
pub fn erfv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.erfh, x);
    if (T == f32) return mapUnary(f32, scalar.erff, x);
    if (T == f64) return mapUnary(f64, scalar.erfd, x);
    @compileError("zmath.erfv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the erf of every lane in `x` with the binary16 scalar kernel.
pub fn erfhv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.erfh, x);
}

/// Computes the erf of every lane in `x` with the binary32 scalar kernel.
pub fn erffv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.erff, x);
}

/// Computes the erf of every lane in `x` with the binary64 scalar kernel.
pub fn erfdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.erfd, x);
}

/// Computes the sin of every lane in `x` with the scalar kernel for its lane type.
pub fn sinv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.sinh, x);
    if (T == f32) return mapUnary(f32, scalar.sinf, x);
    if (T == f64) return mapUnary(f64, scalar.sind, x);
    @compileError("zmath.sinv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the sin of every lane in `x` with the binary32 scalar kernel.
pub fn sinfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.sinf, x);
}

/// Computes the sin of every lane in `x` with the binary64 scalar kernel.
pub fn sindv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.sind, x);
}

/// Computes the sinh of every lane in `x` with the scalar kernel for its lane type.
pub fn sinhv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.sinhh, x);
    if (T == f32) return mapUnary(f32, scalar.sinhf, x);
    if (T == f64) return mapUnary(f64, scalar.sinhd, x);
    @compileError("zmath.sinhv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the sinh of every lane in `x` with the binary16 scalar kernel.
pub fn sinhhv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.sinhh, x);
}

/// Computes the sinh of every lane in `x` with the binary32 scalar kernel.
pub fn sinhfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.sinhf, x);
}

/// Computes the sinh of every lane in `x` with the binary64 scalar kernel.
pub fn sinhdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.sinhd, x);
}

/// Computes the sinpi of every lane in `x` with the scalar kernel for its lane type.
pub fn sinpiv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.sinpih, x);
    if (T == f32) return mapUnary(f32, scalar.sinpif, x);
    if (T == f64) return mapUnary(f64, scalar.sinpid, x);
    @compileError("zmath.sinpiv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the sinpi of every lane in `x` with the binary16 scalar kernel.
pub fn sinpihv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.sinpih, x);
}

/// Computes the sinpi of every lane in `x` with the binary32 scalar kernel.
pub fn sinpifv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.sinpif, x);
}

/// Computes the sinpi of every lane in `x` with the binary64 scalar kernel.
pub fn sinpidv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.sinpid, x);
}

/// Computes the cos of every lane in `x` with the scalar kernel for its lane type.
pub fn cosv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.cosh, x);
    if (T == f32) return mapUnary(f32, scalar.cosf, x);
    if (T == f64) return mapUnary(f64, scalar.cosd, x);
    @compileError("zmath.cosv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the cos of every lane in `x` with the binary32 scalar kernel.
pub fn cosfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.cosf, x);
}

/// Computes the cos of every lane in `x` with the binary64 scalar kernel.
pub fn cosdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.cosd, x);
}

/// Computes the cosh of every lane in `x` with the scalar kernel for its lane type.
pub fn coshv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.coshh, x);
    if (T == f32) return mapUnary(f32, scalar.coshf, x);
    if (T == f64) return mapUnary(f64, scalar.coshd, x);
    @compileError("zmath.coshv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the cosh of every lane in `x` with the binary16 scalar kernel.
pub fn coshhv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.coshh, x);
}

/// Computes the cosh of every lane in `x` with the binary32 scalar kernel.
pub fn coshfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.coshf, x);
}

/// Computes the cosh of every lane in `x` with the binary64 scalar kernel.
pub fn coshdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.coshd, x);
}

/// Computes the cospi of every lane in `x` with the scalar kernel for its lane type.
pub fn cospiv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.cospih, x);
    if (T == f32) return mapUnary(f32, scalar.cospif, x);
    if (T == f64) return mapUnary(f64, scalar.cospid, x);
    @compileError("zmath.cospiv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the cospi of every lane in `x` with the binary16 scalar kernel.
pub fn cospihv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.cospih, x);
}

/// Computes the cospi of every lane in `x` with the binary32 scalar kernel.
pub fn cospifv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.cospif, x);
}

/// Computes the cospi of every lane in `x` with the binary64 scalar kernel.
pub fn cospidv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.cospid, x);
}

/// Computes the acosh of every lane in `x` with the scalar kernel for its lane type.
pub fn acoshv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.acoshh, x);
    if (T == f32) return mapUnary(f32, scalar.acoshf, x);
    if (T == f64) return mapUnary(f64, scalar.acoshd, x);
    @compileError("zmath.acoshv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the acosh of every lane in `x` with the binary16 scalar kernel.
pub fn acoshhv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.acoshh, x);
}

/// Computes the acosh of every lane in `x` with the binary32 scalar kernel.
pub fn acoshfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.acoshf, x);
}

/// Computes the acosh of every lane in `x` with the binary64 scalar kernel.
pub fn acoshdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.acoshd, x);
}

/// Computes the tan of every lane in `x` with the scalar kernel for its lane type.
pub fn tanv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.tanh, x);
    if (T == f32) return mapUnary(f32, scalar.tanf, x);
    if (T == f64) return mapUnary(f64, scalar.tand, x);
    @compileError("zmath.tanv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the tan of every lane in `x` with the binary32 scalar kernel.
pub fn tanfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.tanf, x);
}

/// Computes the tan of every lane in `x` with the binary64 scalar kernel.
pub fn tandv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.tand, x);
}

/// Computes the tanpi of every lane in `x` with the scalar kernel for its lane type.
pub fn tanpiv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.tanpih, x);
    if (T == f32) return mapUnary(f32, scalar.tanpif, x);
    if (T == f64) return mapUnary(f64, scalar.tanpid, x);
    @compileError("zmath.tanpiv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the tanpi of every lane in `x` with the binary16 scalar kernel.
pub fn tanpihv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.tanpih, x);
}

/// Computes the tanpi of every lane in `x` with the binary32 scalar kernel.
pub fn tanpifv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.tanpif, x);
}

/// Computes the tanpi of every lane in `x` with the binary64 scalar kernel.
pub fn tanpidv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.tanpid, x);
}

/// Computes the tanh of every lane in `x` with the scalar kernel for its lane type.
pub fn tanhv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.tanhh, x);
    if (T == f32) return mapUnary(f32, scalar.tanhf, x);
    if (T == f64) return mapUnary(f64, scalar.tanhd, x);
    @compileError("zmath.tanhv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the tanh of every lane in `x` with the binary16 scalar kernel.
pub fn tanhhv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.tanhh, x);
}

/// Computes the tanh of every lane in `x` with the binary32 scalar kernel.
pub fn tanhfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.tanhf, x);
}

/// Computes the tanh of every lane in `x` with the binary64 scalar kernel.
pub fn tanhdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.tanhd, x);
}

/// Computes the rsqrt of every lane in `x` with the scalar kernel for its lane type.
pub fn rsqrtv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.rsqrth, x);
    if (T == f32) return mapUnary(f32, scalar.rsqrtf, x);
    if (T == f64) return mapUnary(f64, scalar.rsqrtd, x);
    @compileError("zmath.rsqrtv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the rsqrt of every lane in `x` with the binary16 scalar kernel.
pub fn rsqrthv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.rsqrth, x);
}

/// Computes the rsqrt of every lane in `x` with the binary32 scalar kernel.
pub fn rsqrtfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.rsqrtf, x);
}

/// Computes the rsqrt of every lane in `x` with the binary64 scalar kernel.
pub fn rsqrtdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.rsqrtd, x);
}

/// Computes the sqrt of every lane in `x` with the scalar kernel for its lane type.
pub fn sqrtv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.sqrth, x);
    if (T == f32) return mapUnary(f32, scalar.sqrtf, x);
    if (T == f64) return mapUnary(f64, scalar.sqrtd, x);
    @compileError("zmath.sqrtv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the sqrt of every lane in `x` with the binary16 scalar kernel.
pub fn sqrthv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.sqrth, x);
}

/// Computes the sqrt of every lane in `x` with the binary32 scalar kernel.
pub fn sqrtfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.sqrtf, x);
}

/// Computes the sqrt of every lane in `x` with the binary64 scalar kernel.
pub fn sqrtdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.sqrtd, x);
}

/// Computes the acos of every lane in `x` with the scalar kernel for its lane type.
pub fn acosv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.acosh, x);
    if (T == f32) return mapUnary(f32, scalar.acosf, x);
    if (T == f64) return mapUnary(f64, scalar.acosd, x);
    @compileError("zmath.acosv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the acos of every lane in `x` with the binary32 scalar kernel.
pub fn acosfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.acosf, x);
}

/// Computes the acos of every lane in `x` with the binary64 scalar kernel.
pub fn acosdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.acosd, x);
}

/// Computes the asin of every lane in `x` with the scalar kernel for its lane type.
pub fn asinv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.asinh, x);
    if (T == f32) return mapUnary(f32, scalar.asinf, x);
    if (T == f64) return mapUnary(f64, scalar.asind, x);
    @compileError("zmath.asinv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the asin of every lane in `x` with the binary32 scalar kernel.
pub fn asinfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.asinf, x);
}

/// Computes the asin of every lane in `x` with the binary64 scalar kernel.
pub fn asindv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.asind, x);
}

/// Computes the asinh of every lane in `x` with the scalar kernel for its lane type.
pub fn asinhv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.asinhh, x);
    if (T == f32) return mapUnary(f32, scalar.asinhf, x);
    if (T == f64) return mapUnary(f64, scalar.asinhd, x);
    @compileError("zmath.asinhv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the asinh of every lane in `x` with the binary16 scalar kernel.
pub fn asinhhv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.asinhh, x);
}

/// Computes the asinh of every lane in `x` with the binary32 scalar kernel.
pub fn asinhfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.asinhf, x);
}

/// Computes the asinh of every lane in `x` with the binary64 scalar kernel.
pub fn asinhdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.asinhd, x);
}

/// Computes the atanh of every lane in `x` with the scalar kernel for its lane type.
pub fn atanhv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.atanhh, x);
    if (T == f32) return mapUnary(f32, scalar.atanhf, x);
    if (T == f64) return mapUnary(f64, scalar.atanhd, x);
    @compileError("zmath.atanhv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the atanh of every lane in `x` with the binary16 scalar kernel.
pub fn atanhhv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.atanhh, x);
}

/// Computes the atanh of every lane in `x` with the binary32 scalar kernel.
pub fn atanhfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.atanhf, x);
}

/// Computes the atanh of every lane in `x` with the binary64 scalar kernel.
pub fn atanhdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.atanhd, x);
}

/// Computes the atan of every lane in `x` with the scalar kernel for its lane type.
pub fn atanv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.atanh, x);
    if (T == f32) return mapUnary(f32, scalar.atanf, x);
    if (T == f64) return mapUnary(f64, scalar.atand, x);
    @compileError("zmath.atanv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the atan of every lane in `x` with the binary32 scalar kernel.
pub fn atanfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.atanf, x);
}

/// Computes the atan of every lane in `x` with the binary64 scalar kernel.
pub fn atandv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.atand, x);
}

/// Computes the cbrt of every lane in `x` with the scalar kernel for its lane type.
pub fn cbrtv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.cbrth, x);
    if (T == f32) return mapUnary(f32, scalar.cbrtf, x);
    if (T == f64) return mapUnary(f64, scalar.cbrtd, x);
    @compileError("zmath.cbrtv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the cbrt of every lane in `x` with the binary16 scalar kernel.
pub fn cbrthv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.cbrth, x);
}

/// Computes the cbrt of every lane in `x` with the binary32 scalar kernel.
pub fn cbrtfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.cbrtf, x);
}

/// Computes the cbrt of every lane in `x` with the binary64 scalar kernel.
pub fn cbrtdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.cbrtd, x);
}

/// Computes the exp2 of every lane in `x` with the scalar kernel for its lane type.
pub fn exp2v(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.exp2h, x);
    if (T == f32) return mapUnary(f32, scalar.exp2f, x);
    if (T == f64) return mapUnary(f64, scalar.exp2d, x);
    @compileError("zmath.exp2v supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the exp2 of every lane in `x` with the binary16 scalar kernel.
pub fn exp2hv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.exp2h, x);
}

/// Computes the exp2 of every lane in `x` with the binary32 scalar kernel.
pub fn exp2fv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.exp2f, x);
}

/// Computes the exp2 of every lane in `x` with the binary64 scalar kernel.
pub fn exp2dv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.exp2d, x);
}

/// Computes the expm1 of every lane in `x` with the scalar kernel for its lane type.
pub fn expm1v(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.expm1h, x);
    if (T == f32) return mapUnary(f32, scalar.expm1f, x);
    if (T == f64) return mapUnary(f64, scalar.expm1d, x);
    @compileError("zmath.expm1v supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the expm1 of every lane in `x` with the binary16 scalar kernel.
pub fn expm1hv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.expm1h, x);
}

/// Computes the expm1 of every lane in `x` with the binary32 scalar kernel.
pub fn expm1fv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.expm1f, x);
}

/// Computes the expm1 of every lane in `x` with the binary64 scalar kernel.
pub fn expm1dv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.expm1d, x);
}

/// Computes the log1p of every lane in `x` with the scalar kernel for its lane type.
pub fn log1pv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.log1ph, x);
    if (T == f32) return mapUnary(f32, scalar.log1pf, x);
    if (T == f64) return mapUnary(f64, scalar.log1pd, x);
    @compileError("zmath.log1pv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the log1p of every lane in `x` with the binary16 scalar kernel.
pub fn log1phv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.log1ph, x);
}

/// Computes the log1p of every lane in `x` with the binary32 scalar kernel.
pub fn log1pfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.log1pf, x);
}

/// Computes the log1p of every lane in `x` with the binary64 scalar kernel.
pub fn log1pdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.log1pd, x);
}

/// Computes the logp1 of every lane in `x` with the scalar kernel for its lane type.
pub fn logp1v(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.logp1h, x);
    if (T == f32) return mapUnary(f32, scalar.logp1f, x);
    if (T == f64) return mapUnary(f64, scalar.logp1d, x);
    @compileError("zmath.logp1v supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the logp1 of every lane in `x` with the binary16 scalar kernel.
pub fn logp1hv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.logp1h, x);
}

/// Computes the logp1 of every lane in `x` with the binary32 scalar kernel.
pub fn logp1fv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.logp1f, x);
}

/// Computes the logp1 of every lane in `x` with the binary64 scalar kernel.
pub fn logp1dv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.logp1d, x);
}

/// Computes the ceil of every lane in `x` with the scalar kernel for its lane type.
pub fn ceilv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.ceilh, x);
    if (T == f32) return mapUnary(f32, scalar.ceilf, x);
    if (T == f64) return mapUnary(f64, scalar.ceild, x);
    @compileError("zmath.ceilv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the ceil of every lane in `x` with the binary16 scalar kernel.
pub fn ceilhv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.ceilh, x);
}

/// Computes the ceil of every lane in `x` with the binary32 scalar kernel.
pub fn ceilfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.ceilf, x);
}

/// Computes the ceil of every lane in `x` with the binary64 scalar kernel.
pub fn ceildv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.ceild, x);
}

/// Computes the floor of every lane in `x` with the scalar kernel for its lane type.
pub fn floorv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.floorh, x);
    if (T == f32) return mapUnary(f32, scalar.floorf, x);
    if (T == f64) return mapUnary(f64, scalar.floord, x);
    @compileError("zmath.floorv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the floor of every lane in `x` with the binary16 scalar kernel.
pub fn floorhv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.floorh, x);
}

/// Computes the floor of every lane in `x` with the binary32 scalar kernel.
pub fn floorfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.floorf, x);
}

/// Computes the floor of every lane in `x` with the binary64 scalar kernel.
pub fn floordv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.floord, x);
}

/// Computes the nearbyint of every lane in `x` with the scalar kernel for its lane type.
pub fn nearbyintv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.nearbyinth, x);
    if (T == f32) return mapUnary(f32, scalar.nearbyintf, x);
    if (T == f64) return mapUnary(f64, scalar.nearbyintd, x);
    @compileError("zmath.nearbyintv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the nearbyint of every lane in `x` with the binary16 scalar kernel.
pub fn nearbyinthv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.nearbyinth, x);
}

/// Computes the nearbyint of every lane in `x` with the binary32 scalar kernel.
pub fn nearbyintfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.nearbyintf, x);
}

/// Computes the nearbyint of every lane in `x` with the binary64 scalar kernel.
pub fn nearbyintdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.nearbyintd, x);
}

/// Computes the rint of every lane in `x` with the scalar kernel for its lane type.
pub fn rintv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.rinth, x);
    if (T == f32) return mapUnary(f32, scalar.rintf, x);
    if (T == f64) return mapUnary(f64, scalar.rintd, x);
    @compileError("zmath.rintv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the rint of every lane in `x` with the binary16 scalar kernel.
pub fn rinthv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.rinth, x);
}

/// Computes the rint of every lane in `x` with the binary32 scalar kernel.
pub fn rintfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.rintf, x);
}

/// Computes the rint of every lane in `x` with the binary64 scalar kernel.
pub fn rintdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.rintd, x);
}

/// Computes the roundeven of every lane in `x` with the scalar kernel for its lane type.
pub fn roundevenv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.roundevenh, x);
    if (T == f32) return mapUnary(f32, scalar.roundevenf, x);
    if (T == f64) return mapUnary(f64, scalar.roundevend, x);
    @compileError("zmath.roundevenv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the roundeven of every lane in `x` with the binary16 scalar kernel.
pub fn roundevenhv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.roundevenh, x);
}

/// Computes the roundeven of every lane in `x` with the binary32 scalar kernel.
pub fn roundevenfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.roundevenf, x);
}

/// Computes the roundeven of every lane in `x` with the binary64 scalar kernel.
pub fn roundevendv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.roundevend, x);
}

/// Computes the trunc of every lane in `x` with the scalar kernel for its lane type.
pub fn truncv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.trunch, x);
    if (T == f32) return mapUnary(f32, scalar.truncf, x);
    if (T == f64) return mapUnary(f64, scalar.truncd, x);
    @compileError("zmath.truncv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the trunc of every lane in `x` with the binary16 scalar kernel.
pub fn trunchv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.trunch, x);
}

/// Computes the trunc of every lane in `x` with the binary32 scalar kernel.
pub fn truncfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.truncf, x);
}

/// Computes the trunc of every lane in `x` with the binary64 scalar kernel.
pub fn truncdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.truncd, x);
}

/// Computes the asinpi of every lane in `x` with the scalar kernel for its lane type.
pub fn asinpiv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.asinpih, x);
    if (T == f32) return mapUnary(f32, scalar.asinpif, x);
    if (T == f64) return mapUnary(f64, scalar.asinpid, x);
    @compileError("zmath.asinpiv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the asinpi of every lane in `x` with the binary16 scalar kernel.
pub fn asinpihv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.asinpih, x);
}

/// Computes the asinpi of every lane in `x` with the binary32 scalar kernel.
pub fn asinpifv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.asinpif, x);
}

/// Computes the asinpi of every lane in `x` with the binary64 scalar kernel.
pub fn asinpidv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.asinpid, x);
}

/// Computes the acospi of every lane in `x` with the scalar kernel for its lane type.
pub fn acospiv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.acospih, x);
    if (T == f32) return mapUnary(f32, scalar.acospif, x);
    if (T == f64) return mapUnary(f64, scalar.acospid, x);
    @compileError("zmath.acospiv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the acospi of every lane in `x` with the binary16 scalar kernel.
pub fn acospihv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.acospih, x);
}

/// Computes the acospi of every lane in `x` with the binary32 scalar kernel.
pub fn acospifv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.acospif, x);
}

/// Computes the acospi of every lane in `x` with the binary64 scalar kernel.
pub fn acospidv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.acospid, x);
}

/// Computes the atanpi of every lane in `x` with the scalar kernel for its lane type.
pub fn atanpiv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.atanpih, x);
    if (T == f32) return mapUnary(f32, scalar.atanpif, x);
    if (T == f64) return mapUnary(f64, scalar.atanpid, x);
    @compileError("zmath.atanpiv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the atanpi of every lane in `x` with the binary16 scalar kernel.
pub fn atanpihv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.atanpih, x);
}

/// Computes the atanpi of every lane in `x` with the binary32 scalar kernel.
pub fn atanpifv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.atanpif, x);
}

/// Computes the atanpi of every lane in `x` with the binary64 scalar kernel.
pub fn atanpidv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.atanpid, x);
}

/// Computes the logb of every lane in `x` with the scalar kernel for its lane type.
pub fn logbv(x: anytype) @TypeOf(x) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnary(f16, scalar.logbh, x);
    if (T == f32) return mapUnary(f32, scalar.logbf, x);
    if (T == f64) return mapUnary(f64, scalar.logbd, x);
    @compileError("zmath.logbv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the logb of every lane in `x` with the binary16 scalar kernel.
pub fn logbhv(x: anytype) @TypeOf(x) {
    return mapUnary(f16, scalar.logbh, x);
}

/// Computes the logb of every lane in `x` with the binary32 scalar kernel.
pub fn logbfv(x: anytype) @TypeOf(x) {
    return mapUnary(f32, scalar.logbf, x);
}

/// Computes the logb of every lane in `x` with the binary64 scalar kernel.
pub fn logbdv(x: anytype) @TypeOf(x) {
    return mapUnary(f64, scalar.logbd, x);
}

/// Computes pow of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn powv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return mapBinary(f16, scalar.powh, a, b);
    if (T == f32) return mapBinary(f32, scalar.powf, a, b);
    if (T == f64) return mapBinary(f64, scalar.powd, a, b);
    @compileError("zmath.powv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes pow of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn powhv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f16, scalar.powh, a, b);
}

/// Computes pow of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn powfv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f32, scalar.powf, a, b);
}

/// Computes pow of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn powdv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f64, scalar.powd, a, b);
}

/// Computes powr of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn powrv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return mapBinary(f16, scalar.powrh, a, b);
    if (T == f32) return mapBinary(f32, scalar.powrf, a, b);
    if (T == f64) return mapBinary(f64, scalar.powrd, a, b);
    @compileError("zmath.powrv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes powr of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn powrhv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f16, scalar.powrh, a, b);
}

/// Computes powr of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn powrfv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f32, scalar.powrf, a, b);
}

/// Computes powr of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn powrdv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f64, scalar.powrd, a, b);
}

/// Computes hypot of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn hypotv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return mapBinary(f16, scalar.hypoth, a, b);
    if (T == f32) return mapBinary(f32, scalar.hypotf, a, b);
    if (T == f64) return mapBinary(f64, scalar.hypotd, a, b);
    @compileError("zmath.hypotv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes hypot of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn hypothv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f16, scalar.hypoth, a, b);
}

/// Computes hypot of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn hypotfv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f32, scalar.hypotf, a, b);
}

/// Computes hypot of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn hypotdv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f64, scalar.hypotd, a, b);
}

/// Computes fmod of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn fmodv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return mapBinary(f16, scalar.fmodh, a, b);
    if (T == f32) return mapBinary(f32, scalar.fmodf, a, b);
    if (T == f64) return mapBinary(f64, scalar.fmodd, a, b);
    @compileError("zmath.fmodv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes fmod of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn fmodhv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f16, scalar.fmodh, a, b);
}

/// Computes fmod of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn fmodfv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f32, scalar.fmodf, a, b);
}

/// Computes fmod of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn fmoddv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f64, scalar.fmodd, a, b);
}

/// Computes remainder of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn remainderv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return mapBinary(f16, scalar.remainderh, a, b);
    if (T == f32) return mapBinary(f32, scalar.remainderf, a, b);
    if (T == f64) return mapBinary(f64, scalar.remainderd, a, b);
    @compileError("zmath.remainderv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes remainder of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn remainderhv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f16, scalar.remainderh, a, b);
}

/// Computes remainder of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn remainderfv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f32, scalar.remainderf, a, b);
}

/// Computes remainder of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn remainderdv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f64, scalar.remainderd, a, b);
}

/// Computes fdim of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn fdimv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return mapBinary(f16, scalar.fdimh, a, b);
    if (T == f32) return mapBinary(f32, scalar.fdimf, a, b);
    if (T == f64) return mapBinary(f64, scalar.fdimd, a, b);
    @compileError("zmath.fdimv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes fdim of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn fdimhv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f16, scalar.fdimh, a, b);
}

/// Computes fdim of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn fdimfv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f32, scalar.fdimf, a, b);
}

/// Computes fdim of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn fdimdv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f64, scalar.fdimd, a, b);
}

/// Computes fmax of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn fmaxv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return mapBinary(f16, scalar.fmaxh, a, b);
    if (T == f32) return mapBinary(f32, scalar.fmaxf, a, b);
    if (T == f64) return mapBinary(f64, scalar.fmaxd, a, b);
    @compileError("zmath.fmaxv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes fmax of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn fmaxhv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f16, scalar.fmaxh, a, b);
}

/// Computes fmax of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn fmaxfv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f32, scalar.fmaxf, a, b);
}

/// Computes fmax of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn fmaxdv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f64, scalar.fmaxd, a, b);
}

/// Computes fmin of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn fminv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return mapBinary(f16, scalar.fminh, a, b);
    if (T == f32) return mapBinary(f32, scalar.fminf, a, b);
    if (T == f64) return mapBinary(f64, scalar.fmind, a, b);
    @compileError("zmath.fminv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes fmin of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn fminhv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f16, scalar.fminh, a, b);
}

/// Computes fmin of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn fminfv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f32, scalar.fminf, a, b);
}

/// Computes fmin of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn fmindv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f64, scalar.fmind, a, b);
}

/// Computes fmaximum of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn fmaximumv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return mapBinary(f16, scalar.fmaximumh, a, b);
    if (T == f32) return mapBinary(f32, scalar.fmaximumf, a, b);
    if (T == f64) return mapBinary(f64, scalar.fmaximumd, a, b);
    @compileError("zmath.fmaximumv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes fmaximum of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn fmaximumhv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f16, scalar.fmaximumh, a, b);
}

/// Computes fmaximum of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn fmaximumfv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f32, scalar.fmaximumf, a, b);
}

/// Computes fmaximum of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn fmaximumdv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f64, scalar.fmaximumd, a, b);
}

/// Computes fminimum of corresponding lanes of `a` and `b` with the scalar kernel for their lane type.
pub fn fminimumv(a: anytype, b: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return mapBinary(f16, scalar.fminimumh, a, b);
    if (T == f32) return mapBinary(f32, scalar.fminimumf, a, b);
    if (T == f64) return mapBinary(f64, scalar.fminimumd, a, b);
    @compileError("zmath.fminimumv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes fminimum of corresponding lanes of `a` and `b` with the binary16 scalar kernel.
pub fn fminimumhv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f16, scalar.fminimumh, a, b);
}

/// Computes fminimum of corresponding lanes of `a` and `b` with the binary32 scalar kernel.
pub fn fminimumfv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f32, scalar.fminimumf, a, b);
}

/// Computes fminimum of corresponding lanes of `a` and `b` with the binary64 scalar kernel.
pub fn fminimumdv(a: anytype, b: anytype) @TypeOf(a) {
    return mapBinary(f64, scalar.fminimumd, a, b);
}
