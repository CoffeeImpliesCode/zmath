//! Guard over the `@Vector` surface of this library.
//!
//! The surface is machine-derived from the vector dispatch lines in `root.zig`,
//! and three defects are latent in that derivation. Each is pinned below.
//!
//! 1. Names collide under the `op + type + v` convention: `sinhv`, `coshv`,
//!    `tanhv`, `asinhv`, `acoshv`, and `atanhv` are the GENERIC entry points for
//!    `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, and `atanh`, so the binary16
//!    specializations of `sin`, `cos`, `tan`, `asin`, `acos`, and `atan` are
//!    intentionally absent.
//! 2. The binary16 kernel of a bare operation is not its hyperbolic twin:
//!    `root.sin` binds `scalar.sinh`, the binary16 SINE, and not `scalar.sinhh`,
//!    the binary16 HYPERBOLIC SINE. The same trap holds for `cos`, `tan`,
//!    `asin`, `acos`, and `atan`.
//! 3. A lane of the float-to-float surface maps to a lane of the SAME float
//!    type, so an operation returning a struct, an integer, a `bool`, or a
//!    narrowed float cannot be spelled `@Vector(N, T) -> @Vector(N, T)` and has
//!    no `OPv` name in `vector.zig`. Its lane-parallel entry point lives in one
//!    of the family modules `vector.zig` re-exports, and is named there.
//!
//! The first table pins the 55 float-to-float operations the surface started
//! with. The rest of the file pins the operations vectorized afterwards, in the
//! four family modules, and closes with a completeness check that ties the
//! pinned set to the vectorized `root.zig` surface in BOTH directions, so a new
//! dispatch line without a table row, or a table row for an operation that is
//! not vectorized, fails to compile.

const std = @import("std");
const testing = std.testing;
const root = @import("root.zig");
const scalar = @import("scalar.zig");
const vector = @import("vector.zig");

// Short module aliases keep every row of the surface tables on a single line.
const s = scalar;
const v = vector;

/// Returns the raw bits of a binary16 value, so comparisons see signed zero and
/// NaN payloads exactly.
fn bits16(x: f16) u16 {
    return @bitCast(x);
}

/// Asserts lane `i` of `vop(x)` equals `skernel(x[i])`, compared by raw bits.
fn expectUnaryLanes(
    comptime T: type,
    comptime Bits: type,
    comptime vop: anytype,
    comptime skernel: anytype,
    x: @Vector(4, T),
) !void {
    const got: @Vector(4, T) = vop(x);
    inline for (0..4) |i| {
        const want: T = skernel(x[i]);
        try testing.expectEqual(@as(Bits, @bitCast(want)), @as(Bits, @bitCast(got[i])));
    }
}

/// Asserts lane `i` of `vop(a, b)` equals `skernel(a[i], b[i])`, by raw bits.
fn expectBinaryLanes(
    comptime T: type,
    comptime Bits: type,
    comptime vop: anytype,
    comptime skernel: anytype,
    a: @Vector(4, T),
    b: @Vector(4, T),
) !void {
    const got: @Vector(4, T) = vop(a, b);
    inline for (0..4) |i| {
        const want: T = skernel(a[i], b[i]);
        try testing.expectEqual(@as(Bits, @bitCast(want)), @as(Bits, @bitCast(got[i])));
    }
}

/// Checks one unary operation against the binary16, binary32, and binary64
/// kernels that `root.zig` binds for it. The three vector slots take the
/// generic `opv` entry point; a specialized `*hv`, `*fv`, or `*dv` name drops
/// into the same slot.
fn checkUnary(
    comptime v16: anytype,
    comptime v32: anytype,
    comptime v64: anytype,
    comptime s16: anytype,
    comptime s32: anytype,
    comptime s64: anytype,
) !void {
    // Two four-lane vectors per format: signed zero and unity first, then the
    // fractional lanes plus the infinities and NaN that exercise each kernel's
    // special-case dispatch.
    const x16a: @Vector(4, f16) = .{ 0.0, -0.0, 1.0, -1.0 };
    const x16b: @Vector(4, f16) = .{ 0.5, -0.5, 1.5, -1.5 };
    const x32a: @Vector(4, f32) = .{ 0.0, -0.0, 1.0, -1.0 };
    const x32b: @Vector(4, f32) = .{ 0.5, -0.5, std.math.inf(f32), std.math.nan(f32) };
    const x64a: @Vector(4, f64) = .{ 0.0, -0.0, 1.0, -1.0 };
    const x64b: @Vector(4, f64) = .{ 0.5, -0.5, std.math.inf(f64), std.math.nan(f64) };

    try expectUnaryLanes(f16, u16, v16, s16, x16a);
    try expectUnaryLanes(f16, u16, v16, s16, x16b);
    try expectUnaryLanes(f32, u32, v32, s32, x32a);
    try expectUnaryLanes(f32, u32, v32, s32, x32b);
    try expectUnaryLanes(f64, u64, v64, s64, x64a);
    try expectUnaryLanes(f64, u64, v64, s64, x64b);
}

/// Checks one binary operation against the binary16, binary32, and binary64
/// kernels that `root.zig` binds for it.
fn checkBinary(
    comptime v16: anytype,
    comptime v32: anytype,
    comptime v64: anytype,
    comptime s16: anytype,
    comptime s32: anytype,
    comptime s64: anytype,
) !void {
    const a16a: @Vector(4, f16) = .{ 0.0, -0.0, 1.0, -1.0 };
    const b16a: @Vector(4, f16) = .{ 1.0, 0.5, 2.0, 2.0 };
    const a16b: @Vector(4, f16) = .{ 0.5, -0.5, 1.5, -1.5 };
    const b16b: @Vector(4, f16) = .{ 2.0, 1.5, 0.5, 0.5 };
    const a32a: @Vector(4, f32) = .{ 0.0, -0.0, 1.0, -1.0 };
    const b32a: @Vector(4, f32) = .{ 1.0, 0.5, 2.0, 2.0 };
    const a32b: @Vector(4, f32) = .{ 0.5, -0.5, 1.0, -1.0 };
    const b32b: @Vector(4, f32) = .{ std.math.inf(f32), 2.0, 2.0, std.math.nan(f32) };
    const a64a: @Vector(4, f64) = .{ 0.0, -0.0, 1.0, -1.0 };
    const b64a: @Vector(4, f64) = .{ 1.0, 0.5, 2.0, 2.0 };
    const a64b: @Vector(4, f64) = .{ 0.5, -0.5, 1.0, -1.0 };
    const b64b: @Vector(4, f64) = .{ std.math.inf(f64), 2.0, 2.0, std.math.nan(f64) };

    try expectBinaryLanes(f16, u16, v16, s16, a16a, b16a);
    try expectBinaryLanes(f16, u16, v16, s16, a16b, b16b);
    try expectBinaryLanes(f32, u32, v32, s32, a32a, b32a);
    try expectBinaryLanes(f32, u32, v32, s32, a32b, b32b);
    try expectBinaryLanes(f64, u64, v64, s64, a64a, b64a);
    try expectBinaryLanes(f64, u64, v64, s64, a64b, b64b);
}

test "every vectorized operation applies its root scalar kernel lane for lane" {
    // Each row is one operation, its generic vector entry point, and the three
    // scalar kernels `root.zig` binds for that operation. 45 unary, then 10
    // binary, in alphabetical order.
    try checkUnary(v.absv, v.absv, v.absv, s.absh, s.absf, s.absd);
    try checkUnary(v.acosv, v.acosv, v.acosv, s.acosh, s.acosf, s.acosd);
    try checkUnary(v.acoshv, v.acoshv, v.acoshv, s.acoshh, s.acoshf, s.acoshd);
    try checkUnary(v.acospiv, v.acospiv, v.acospiv, s.acospih, s.acospif, s.acospid);
    try checkUnary(v.asinv, v.asinv, v.asinv, s.asinh, s.asinf, s.asind);
    try checkUnary(v.asinhv, v.asinhv, v.asinhv, s.asinhh, s.asinhf, s.asinhd);
    try checkUnary(v.asinpiv, v.asinpiv, v.asinpiv, s.asinpih, s.asinpif, s.asinpid);
    try checkUnary(v.atanv, v.atanv, v.atanv, s.atanh, s.atanf, s.atand);
    try checkUnary(v.atanhv, v.atanhv, v.atanhv, s.atanhh, s.atanhf, s.atanhd);
    try checkUnary(v.atanpiv, v.atanpiv, v.atanpiv, s.atanpih, s.atanpif, s.atanpid);
    try checkUnary(v.cbrtv, v.cbrtv, v.cbrtv, s.cbrth, s.cbrtf, s.cbrtd);
    try checkUnary(v.ceilv, v.ceilv, v.ceilv, s.ceilh, s.ceilf, s.ceild);
    try checkUnary(v.cosv, v.cosv, v.cosv, s.cosh, s.cosf, s.cosd);
    try checkUnary(v.coshv, v.coshv, v.coshv, s.coshh, s.coshf, s.coshd);
    try checkUnary(v.cospiv, v.cospiv, v.cospiv, s.cospih, s.cospif, s.cospid);
    try checkUnary(v.erfv, v.erfv, v.erfv, s.erfh, s.erff, s.erfd);
    try checkUnary(v.erfcv, v.erfcv, v.erfcv, s.erfch, s.erfcf, s.erfcd);
    try checkUnary(v.expv, v.expv, v.expv, s.exph, s.expf, s.expd);
    try checkUnary(v.exp10v, v.exp10v, v.exp10v, s.exp10h, s.exp10f, s.exp10d);
    try checkUnary(v.exp10m1v, v.exp10m1v, v.exp10m1v, s.exp10m1h, s.exp10m1f, s.exp10m1d);
    try checkUnary(v.exp2v, v.exp2v, v.exp2v, s.exp2h, s.exp2f, s.exp2d);
    try checkUnary(v.exp2m1v, v.exp2m1v, v.exp2m1v, s.exp2m1h, s.exp2m1f, s.exp2m1d);
    try checkUnary(v.expm1v, v.expm1v, v.expm1v, s.expm1h, s.expm1f, s.expm1d);
    try checkUnary(v.floorv, v.floorv, v.floorv, s.floorh, s.floorf, s.floord);
    try checkUnary(v.logv, v.logv, v.logv, s.logh, s.logf, s.logd);
    try checkUnary(v.log10v, v.log10v, v.log10v, s.log10h, s.log10f, s.log10d);
    try checkUnary(v.log10p1v, v.log10p1v, v.log10p1v, s.log10p1h, s.log10p1f, s.log10p1d);
    try checkUnary(v.log1pv, v.log1pv, v.log1pv, s.log1ph, s.log1pf, s.log1pd);
    try checkUnary(v.log2v, v.log2v, v.log2v, s.log2h, s.log2f, s.log2d);
    try checkUnary(v.log2p1v, v.log2p1v, v.log2p1v, s.log2p1h, s.log2p1f, s.log2p1d);
    try checkUnary(v.logbv, v.logbv, v.logbv, s.logbh, s.logbf, s.logbd);
    try checkUnary(v.logp1v, v.logp1v, v.logp1v, s.logp1h, s.logp1f, s.logp1d);
    try checkUnary(
        v.nearbyintv,
        v.nearbyintv,
        v.nearbyintv,
        s.nearbyinth,
        s.nearbyintf,
        s.nearbyintd,
    );
    try checkUnary(v.rintv, v.rintv, v.rintv, s.rinth, s.rintf, s.rintd);
    try checkUnary(
        v.roundevenv,
        v.roundevenv,
        v.roundevenv,
        s.roundevenh,
        s.roundevenf,
        s.roundevend,
    );
    try checkUnary(v.rsqrtv, v.rsqrtv, v.rsqrtv, s.rsqrth, s.rsqrtf, s.rsqrtd);
    try checkUnary(v.sinv, v.sinv, v.sinv, s.sinh, s.sinf, s.sind);
    try checkUnary(v.sinhv, v.sinhv, v.sinhv, s.sinhh, s.sinhf, s.sinhd);
    try checkUnary(v.sinpiv, v.sinpiv, v.sinpiv, s.sinpih, s.sinpif, s.sinpid);
    try checkUnary(v.sqrtv, v.sqrtv, v.sqrtv, s.sqrth, s.sqrtf, s.sqrtd);
    try checkUnary(v.tanv, v.tanv, v.tanv, s.tanh, s.tanf, s.tand);
    try checkUnary(v.tanhv, v.tanhv, v.tanhv, s.tanhh, s.tanhf, s.tanhd);
    try checkUnary(v.tanpiv, v.tanpiv, v.tanpiv, s.tanpih, s.tanpif, s.tanpid);
    try checkUnary(v.tgammav, v.tgammav, v.tgammav, s.tgammah, s.tgammaf, s.tgammad);
    try checkUnary(v.truncv, v.truncv, v.truncv, s.trunch, s.truncf, s.truncd);

    try checkBinary(v.fdimv, v.fdimv, v.fdimv, s.fdimh, s.fdimf, s.fdimd);
    try checkBinary(v.fmaxv, v.fmaxv, v.fmaxv, s.fmaxh, s.fmaxf, s.fmaxd);
    try checkBinary(v.fmaximumv, v.fmaximumv, v.fmaximumv, s.fmaximumh, s.fmaximumf, s.fmaximumd);
    try checkBinary(v.fminv, v.fminv, v.fminv, s.fminh, s.fminf, s.fmind);
    try checkBinary(v.fminimumv, v.fminimumv, v.fminimumv, s.fminimumh, s.fminimumf, s.fminimumd);
    try checkBinary(v.fmodv, v.fmodv, v.fmodv, s.fmodh, s.fmodf, s.fmodd);
    try checkBinary(v.hypotv, v.hypotv, v.hypotv, s.hypoth, s.hypotf, s.hypotd);
    try checkBinary(v.powv, v.powv, v.powv, s.powh, s.powf, s.powd);
    try checkBinary(v.powrv, v.powrv, v.powrv, s.powrh, s.powrf, s.powrd);
    try checkBinary(
        v.remainderv,
        v.remainderv,
        v.remainderv,
        s.remainderh,
        s.remainderf,
        s.remainderd,
    );
}

test "defect 1: the six colliding names are the hyperbolic generic entry points" {
    // Under the `op + type + v` convention the binary16 specialization of `sin`
    // would be named `sinhv`, which is already the GENERIC entry point for
    // `sinh`. The same collision exists for `cos`/`cosh`, `tan`/`tanh`,
    // `asin`/`asinh`, `acos`/`acosh`, and `atan`/`atanh`, so the binary16
    // SPECIALIZATIONS of `sin`, `cos`, `tan`, `asin`, `acos`, and `atan` are
    // intentionally absent from `vector.zig`. Binary16 callers reach those six
    // kernels through the generic `sinv`, `cosv`, `tanv`, `asinv`, `acosv`, and
    // `atanv` instead.

    // Domain-clean lanes per kernel: `acosh` needs x >= 1, `atanh` needs
    // |x| < 1, and the rest take any real value including signed zero.
    const sym: @Vector(4, f16) = .{ 0.0, -0.0, 1.5, -1.5 };
    const unit: @Vector(4, f16) = .{ 0.0, -0.0, 0.5, -0.5 };
    const pos: @Vector(4, f16) = .{ 1.0, 1.5, 2.0, 2.5 };

    // The hyperbolic generics carry the doubled-`h` binary16 kernels.
    try expectUnaryLanes(f16, u16, v.sinhv, s.sinhh, sym);
    try expectUnaryLanes(f16, u16, v.coshv, s.coshh, sym);
    try expectUnaryLanes(f16, u16, v.tanhv, s.tanhh, sym);
    try expectUnaryLanes(f16, u16, v.asinhv, s.asinhh, sym);
    try expectUnaryLanes(f16, u16, v.atanhv, s.atanhh, unit);
    try expectUnaryLanes(f16, u16, v.acoshv, s.acoshh, pos);

    // The hyperbolic specializations do exist, and use the same kernels.
    try expectUnaryLanes(f16, u16, v.sinhhv, s.sinhh, sym);
    try expectUnaryLanes(f16, u16, v.coshhv, s.coshh, sym);
    try expectUnaryLanes(f16, u16, v.tanhhv, s.tanhh, sym);
    try expectUnaryLanes(f16, u16, v.asinhhv, s.asinhh, sym);
    try expectUnaryLanes(f16, u16, v.atanhhv, s.atanhh, unit);
    try expectUnaryLanes(f16, u16, v.acoshhv, s.acoshh, pos);

    // The hazard is concrete, not a naming nicety: at 0.5 every hyperbolic
    // result differs from the bare result, so a caller who reads `sinhv` as
    // "sine of a binary16 vector" silently gets hyperbolic sine.
    const half: @Vector(4, f16) = @splat(0.5);
    try testing.expect(bits16(v.sinhv(half)[0]) != bits16(v.sinv(half)[0]));
    try testing.expect(bits16(v.coshv(half)[0]) != bits16(v.cosv(half)[0]));
    try testing.expect(bits16(v.tanhv(half)[0]) != bits16(v.tanv(half)[0]));
    try testing.expect(bits16(v.asinhv(half)[0]) != bits16(v.asinv(half)[0]));
    try testing.expect(bits16(v.atanhv(half)[0]) != bits16(v.atanv(half)[0]));
    try testing.expect(bits16(v.acoshv(half)[0]) != bits16(v.acosv(half)[0]));
}

/// Asserts lane `i` of the vector overload `f` equals `f` on that one lane.
fn expectRootLanes(comptime f: anytype, x: @Vector(4, f16)) !void {
    const got: @Vector(4, f16) = f(x);
    inline for (0..4) |i| {
        const want: f16 = f(x[i]);
        try testing.expectEqual(bits16(want), bits16(got[i]));
    }
}

test "defect 2: the bare binary16 root overloads reach the non-hyperbolic kernels" {
    // The binary16 kernel of `sin` is `scalar.sinh`, the binary16 SINE, and not
    // `scalar.sinhh`, the binary16 HYPERBOLIC SINE. The same trap applies to
    // `cos` -> `scalar.cosh`, `tan` -> `scalar.tanh`, `asin` -> `scalar.asinh`,
    // `acos` -> `scalar.acosh`, and `atan` -> `scalar.atanh`. A generator that
    // assumes `op` maps to `op + h` binds the wrong kernel and silently
    // computes the wrong function. Drive the whole public path: the top-level
    // root overload on a binary16 vector must equal that overload applied to
    // one lane at a time. 1.0 is the last in-domain lane for `asin` and `acos`.
    const x: @Vector(4, f16) = .{ 0.0, -0.0, 0.5, 1.0 };
    try expectRootLanes(root.sin, x);
    try expectRootLanes(root.cos, x);
    try expectRootLanes(root.tan, x);
    try expectRootLanes(root.asin, x);
    try expectRootLanes(root.acos, x);
    try expectRootLanes(root.atan, x);
}

test "defect 3: struct- and integer-returning operations have no float vector entry point" {
    // The float vector surface in `vector.zig` is exactly
    // `@Vector(N, T) -> @Vector(N, T)`. An operation that returns a struct
    // (`Atan2Result`, `PowerResult`, `LgammaResult`, `RemquoResult`,
    // `FrexpResult`, `ModfResult`, `IntegerResult`, `PayloadResult`,
    // `CanonicalizeResult`), a bare integer (`ilogb`, `llogb`), a `bool`
    // (`totalOrder`, `totalOrderMag`), or a narrowed `f32` (`fadd` and its five
    // siblings) cannot be spelled that way, so their `v` names must stay absent
    // HERE. Their lane-parallel entry points live in the family modules that
    // `vector.zig` re-exports, and the tables further down pin those.
    comptime {
        if (!@hasDecl(vector, "absv")) @compileError("vector.absv must exist");
        if (!@hasDecl(vector, "powv")) @compileError("vector.powv must exist");
        if (!@hasDecl(vector, "sinv")) @compileError("vector.sinv must exist");
        if (@hasDecl(vector, "atan2v")) @compileError("atan2 returns Atan2Result");
        if (@hasDecl(vector, "atan2piv")) @compileError("atan2pi returns Atan2Result");
        if (@hasDecl(vector, "pownv")) @compileError("pown returns PowerResult");
        if (@hasDecl(vector, "compoundnv")) @compileError("compoundn returns PowerResult");
        if (@hasDecl(vector, "rootnv")) @compileError("rootn returns PowerResult");
        if (@hasDecl(vector, "lgammav")) @compileError("lgamma returns LgammaResult");
        if (@hasDecl(vector, "remquov")) @compileError("remquo returns RemquoResult");
        if (@hasDecl(vector, "frexpv")) @compileError("frexp returns FrexpResult");
        if (@hasDecl(vector, "modfv")) @compileError("modf returns ModfResult");
        if (@hasDecl(vector, "canonicalizev")) @compileError("canonicalize returns a struct");
        if (@hasDecl(vector, "setPayloadv")) @compileError("setPayload returns PayloadResult");
        if (@hasDecl(vector, "setPayloadSigv")) @compileError("setPayloadSig returns a struct");
        if (@hasDecl(vector, "lrintv")) @compileError("lrint returns IntegerResult");
        if (@hasDecl(vector, "llrintv")) @compileError("llrint returns IntegerResult");
        if (@hasDecl(vector, "lroundv")) @compileError("lround returns IntegerResult");
        if (@hasDecl(vector, "llroundv")) @compileError("llround returns IntegerResult");
        if (@hasDecl(vector, "ilogbv")) @compileError("ilogb returns an integer");
        if (@hasDecl(vector, "llogbv")) @compileError("llogb returns an integer");
        if (@hasDecl(vector, "totalOrderv")) @compileError("totalOrder returns a bool");
        if (@hasDecl(vector, "totalOrderMagv")) @compileError("totalOrderMag returns a bool");
    }

    // The scalar path for a struct-returning operation is unchanged: each
    // format still carries its own kernel's exact bits and resolution flag.
    const a16 = root.atan2(@as(f16, 1.0), @as(f16, 1.0));
    try expectResultBits(f16, u16, a16, s.atan2h(1.0, 1.0));
    const a32 = root.atan2(@as(f32, 1.0), @as(f32, 1.0));
    try expectResultBits(f32, u32, a32, s.atan2f(1.0, 1.0));
    const a64 = root.atan2(@as(f64, 1.0), @as(f64, 1.0));
    try expectResultBits(f64, u64, a64, s.atan2d(1.0, 1.0));
    const p32 = root.atan2pi(@as(f32, 1.0), @as(f32, 1.0));
    try expectResultBits(f32, u32, p32, s.atan2pif(1.0, 1.0));
}

/// Asserts an `Atan2Result` carries its kernel's exact bits and status flag.
fn expectResultBits(
    comptime T: type,
    comptime Bits: type,
    got: root.Atan2Result(T),
    want: root.Atan2Result(T),
) !void {
    try testing.expectEqual(want.resolved, got.resolved);
    try testing.expectEqual(@as(Bits, @bitCast(want.value)), @as(Bits, @bitCast(got.value)));
}

// ---------------------------------------------------------------------------
// The lane-parallel surface added after the original float-to-float operations
//
// The table test above pins 55 operations, all of them `@Vector(N, T) ->
// @Vector(N, T)`. The 47 public names vectorized afterwards live in the four
// family modules, where a lane can map to an integer, to a `bool`, to a result
// struct, to a narrowed float, or to a float under an extra non-lane argument.
// They are 41 `pub fn` operations and the six `pub const` extrema aliases,
// which `root.zig` exports under snake_case names while their entry points are
// camelCase. A future edit could re-bind one of them to the wrong kernel or to
// the wrong result type, and the table test would stay green.
//
// Every check below therefore compares a family entry point against
// `root.zig`'s OWN operation, lane for lane, never against a hardcoded kernel
// name. The binding between the two is what is under test. Floats are
// compared by raw bits, integers and `bool`s directly, and a result struct
// field by field.
// ---------------------------------------------------------------------------

/// The lane count of every table below. Four keeps one operation per line and
/// leaves room for a non-uniform input case.
const lane_count = 4;

/// Asserts lane `i` of `f(x)` equals `f` on that lane alone, by raw bits.
fn expectRootUnaryLanes(
    comptime T: type,
    comptime Bits: type,
    comptime f: anytype,
    x: @Vector(lane_count, T),
) !void {
    const got: @Vector(lane_count, T) = f(x);
    inline for (0..lane_count) |i| {
        const want: T = f(x[i]);
        try testing.expectEqual(@as(Bits, @bitCast(want)), @as(Bits, @bitCast(got[i])));
    }
}

/// Asserts lane `i` of `f(a, b)` equals `f` on that lane alone, by raw bits.
fn expectRootBinaryLanes(
    comptime T: type,
    comptime Bits: type,
    comptime f: anytype,
    a: @Vector(lane_count, T),
    b: @Vector(lane_count, T),
) !void {
    const got: @Vector(lane_count, T) = f(a, b);
    inline for (0..lane_count) |i| {
        const want: T = f(a[i], b[i]);
        try testing.expectEqual(@as(Bits, @bitCast(want)), @as(Bits, @bitCast(got[i])));
    }
}

/// Asserts lane `i` of `f(a, b, c)` equals `f` on that lane alone, by raw bits.
fn expectRootTernaryLanes(
    comptime T: type,
    comptime Bits: type,
    comptime f: anytype,
    a: @Vector(lane_count, T),
    b: @Vector(lane_count, T),
    c: @Vector(lane_count, T),
) !void {
    const got: @Vector(lane_count, T) = f(a, b, c);
    inline for (0..lane_count) |i| {
        const want: T = f(a[i], b[i], c[i]);
        try testing.expectEqual(@as(Bits, @bitCast(want)), @as(Bits, @bitCast(got[i])));
    }
}

/// Asserts lane `i` of the integer vector `f(x)` equals `f` on that lane.
fn expectRootIntUnaryLanes(
    comptime I: type,
    comptime T: type,
    comptime f: anytype,
    x: @Vector(lane_count, T),
) !void {
    const got: @Vector(lane_count, I) = f(x);
    inline for (0..lane_count) |i| try testing.expectEqual(f(x[i]), got[i]);
}

/// Asserts lane `i` of the `bool` mask `f(a, b)` equals `f` on that lane, and
/// that the mask holds BOTH a true and a false lane. A mask that came out
/// uniformly true or uniformly false would make the lane comparison vacuous.
fn expectMixedBoolLanes(
    comptime T: type,
    comptime f: anytype,
    a: @Vector(lane_count, T),
    b: @Vector(lane_count, T),
) !void {
    const got: @Vector(lane_count, bool) = f(a, b);
    try testing.expect(@TypeOf(got) == @Vector(lane_count, bool));
    inline for (0..lane_count) |i| try testing.expectEqual(f(a[i], b[i]), got[i]);
    var trues: usize = 0;
    inline for (0..lane_count) |i| {
        if (got[i]) trues += 1;
    }
    try testing.expect(trues > 0);
    try testing.expect(trues < lane_count);
}

/// Compares two lane results of the same result type `R`, field by field. A
/// float field is compared by raw bits, so a signed zero or a NaN payload
/// cannot slip through; every other field is compared directly.
fn expectLaneResult(comptime R: type, got: R, want: R) !void {
    // The field list is forced to comptime so the field names below resolve at
    // compile time while the field VALUES stay runtime.
    const fields = comptime std.meta.declarations(R);
    inline for (fields, 0..) |_, i| {
        const name = fields[i].name;
        const field_type = @field(R, name);
        const got_field = @field(got, name);
        const want_field = @field(want, name);
        switch (@typeInfo(field_type)) {
            .float => {
                const Bits = @Int(.unsigned, @bitSizeOf(field_type));
                try testing.expectEqual(
                    @as(Bits, @bitCast(want_field)),
                    @as(Bits, @bitCast(got_field)),
                );
            },
            else => try testing.expectEqual(want_field, got_field),
        }
    }
}

/// Asserts lane `i` of the lane-ordered `[N]Result` that `f(x)` returns
/// equals `f` on that lane alone, field by field.
fn expectRootUnaryStructLanes(
    comptime R: type,
    comptime T: type,
    comptime f: anytype,
    x: @Vector(lane_count, T),
) !void {
    const got: [lane_count]R = f(x);
    inline for (0..lane_count) |i| try expectLaneResult(R, got[i], f(x[i]));
}

/// Asserts lane `i` of the lane-ordered `[N]Result` that `f(a, b)` returns
/// equals `f` on that lane alone, field by field.
fn expectRootBinaryStructLanes(
    comptime R: type,
    comptime T: type,
    comptime f: anytype,
    a: @Vector(lane_count, T),
    b: @Vector(lane_count, T),
) !void {
    const got: [lane_count]R = f(a, b);
    inline for (0..lane_count) |i| try expectLaneResult(R, got[i], f(a[i], b[i]));
}

/// Runs one float-to-float root operation over the binary16, binary32, and
/// binary64 lanes it dispatches to.
fn checkRootUnary(
    comptime f: anytype,
    x16: @Vector(lane_count, f16),
    x32: @Vector(lane_count, f32),
    x64: @Vector(lane_count, f64),
) !void {
    try expectRootUnaryLanes(f16, u16, f, x16);
    try expectRootUnaryLanes(f32, u32, f, x32);
    try expectRootUnaryLanes(f64, u64, f, x64);
}

/// Runs one binary float-to-float root operation over the binary16, binary32,
/// and binary64 lanes it dispatches to.
fn checkRootBinary(
    comptime f: anytype,
    a16: @Vector(lane_count, f16),
    b16: @Vector(lane_count, f16),
    a32: @Vector(lane_count, f32),
    b32: @Vector(lane_count, f32),
    a64: @Vector(lane_count, f64),
    b64: @Vector(lane_count, f64),
) !void {
    try expectRootBinaryLanes(f16, u16, f, a16, b16);
    try expectRootBinaryLanes(f32, u32, f, a32, b32);
    try expectRootBinaryLanes(f64, u64, f, a64, b64);
}

/// Runs one ternary float-to-float root operation over the binary16,
/// binary32, and binary64 lanes it dispatches to.
fn checkRootTernary(
    comptime f: anytype,
    a16: @Vector(lane_count, f16),
    b16: @Vector(lane_count, f16),
    c16: @Vector(lane_count, f16),
    a32: @Vector(lane_count, f32),
    b32: @Vector(lane_count, f32),
    c32: @Vector(lane_count, f32),
    a64: @Vector(lane_count, f64),
    b64: @Vector(lane_count, f64),
    c64: @Vector(lane_count, f64),
) !void {
    try expectRootTernaryLanes(f16, u16, f, a16, b16, c16);
    try expectRootTernaryLanes(f32, u32, f, a32, b32, c32);
    try expectRootTernaryLanes(f64, u64, f, a64, b64, c64);
}

/// Runs one integer-returning root operation over the binary16, binary32,
/// and binary64 lanes it dispatches to.
fn checkRootIntUnary(
    comptime I: type,
    comptime f: anytype,
    x16: @Vector(lane_count, f16),
    x32: @Vector(lane_count, f32),
    x64: @Vector(lane_count, f64),
) !void {
    try expectRootIntUnaryLanes(I, f16, f, x16);
    try expectRootIntUnaryLanes(I, f32, f, x32);
    try expectRootIntUnaryLanes(I, f64, f, x64);
}

/// Runs one boolean root operation over the binary16, binary32, and binary64
/// lanes it dispatches to, rejecting a mask that is not mixed.
fn checkRootBoolBinary(
    comptime f: anytype,
    a16: @Vector(lane_count, f16),
    b16: @Vector(lane_count, f16),
    a32: @Vector(lane_count, f32),
    b32: @Vector(lane_count, f32),
    a64: @Vector(lane_count, f64),
    b64: @Vector(lane_count, f64),
) !void {
    try expectMixedBoolLanes(f16, f, a16, b16);
    try expectMixedBoolLanes(f32, f, a32, b32);
    try expectMixedBoolLanes(f64, f, a64, b64);
}

/// Runs one result-struct operation whose result type is `R(T)` for the lane
/// type `T`, over the binary16, binary32, and binary64 lanes.
fn checkRootUnaryStruct(
    comptime f: anytype,
    comptime R: anytype,
    x16: @Vector(lane_count, f16),
    x32: @Vector(lane_count, f32),
    x64: @Vector(lane_count, f64),
) !void {
    try expectRootUnaryStructLanes(R(f16), f16, f, x16);
    try expectRootUnaryStructLanes(R(f32), f32, f, x32);
    try expectRootUnaryStructLanes(R(f64), f64, f, x64);
}

/// The same as `checkRootUnaryStruct` for a result type that does not depend
/// on the lane type, such as `IntegerResult`.
fn checkRootUnaryFixedStruct(
    comptime R: type,
    comptime f: anytype,
    x16: @Vector(lane_count, f16),
    x32: @Vector(lane_count, f32),
    x64: @Vector(lane_count, f64),
) !void {
    try expectRootUnaryStructLanes(R, f16, f, x16);
    try expectRootUnaryStructLanes(R, f32, f, x32);
    try expectRootUnaryStructLanes(R, f64, f, x64);
}

/// Runs one binary result-struct operation whose result type is `R(T)` for the
/// lane type `T`, over the binary16, binary32, and binary64 lanes.
fn checkRootBinaryStruct(
    comptime f: anytype,
    comptime R: anytype,
    a16: @Vector(lane_count, f16),
    b16: @Vector(lane_count, f16),
    a32: @Vector(lane_count, f32),
    b32: @Vector(lane_count, f32),
    a64: @Vector(lane_count, f64),
    b64: @Vector(lane_count, f64),
) !void {
    try expectRootBinaryStructLanes(R(f16), f16, f, a16, b16);
    try expectRootBinaryStructLanes(R(f32), f32, f, a32, b32);
    try expectRootBinaryStructLanes(R(f64), f64, f, a64, b64);
}

/// Asserts the narrowing entry point `vop` turns `@Vector(N, f64)` lanes into
/// a `@Vector(N, f32)` result, and that lane `i` equals the root operation
/// `sop` on that lane, compared by raw bits.
fn expectNarrowUnary(
    comptime vop: anytype,
    comptime sop: anytype,
    x: @Vector(lane_count, f64),
) !void {
    const got = vop(x);
    try testing.expect(@TypeOf(x) == @Vector(lane_count, f64));
    try testing.expect(@TypeOf(got) == @Vector(lane_count, f32));
    inline for (0..lane_count) |i| {
        const want: f32 = sop(x[i]);
        try testing.expectEqual(@as(u32, @bitCast(want)), @as(u32, @bitCast(got[i])));
    }
}

/// The binary narrowing counterpart of `expectNarrowUnary`.
fn expectNarrowBinary(
    comptime vop: anytype,
    comptime sop: anytype,
    a: @Vector(lane_count, f64),
    b: @Vector(lane_count, f64),
) !void {
    const got = vop(a, b);
    try testing.expect(@TypeOf(a) == @Vector(lane_count, f64));
    try testing.expect(@TypeOf(got) == @Vector(lane_count, f32));
    inline for (0..lane_count) |i| {
        const want: f32 = sop(a[i], b[i]);
        try testing.expectEqual(@as(u32, @bitCast(want)), @as(u32, @bitCast(got[i])));
    }
}

/// The ternary narrowing counterpart of `expectNarrowUnary`.
fn expectNarrowTernary(
    comptime vop: anytype,
    comptime sop: anytype,
    a: @Vector(lane_count, f64),
    b: @Vector(lane_count, f64),
    c: @Vector(lane_count, f64),
) !void {
    const got = vop(a, b, c);
    try testing.expect(@TypeOf(a) == @Vector(lane_count, f64));
    try testing.expect(@TypeOf(got) == @Vector(lane_count, f32));
    inline for (0..lane_count) |i| {
        const want: f32 = sop(a[i], b[i], c[i]);
        try testing.expectEqual(@as(u32, @bitCast(want)), @as(u32, @bitCast(got[i])));
    }
}

/// Asserts lane `i` of the broadcast result `f(x, e)` equals the root
/// operation `f` on that one lane with the same extra argument.
fn expectBroadcastLanes(
    comptime I: type,
    comptime T: type,
    comptime Bits: type,
    comptime f: anytype,
    x: @Vector(lane_count, T),
    e: I,
) !void {
    const got: @Vector(lane_count, T) = f(x, e);
    inline for (0..lane_count) |i| {
        const want: T = f(x[i], e);
        try testing.expectEqual(@as(Bits, @bitCast(want)), @as(Bits, @bitCast(got[i])));
    }
}

/// Asserts lane `i` of the per-lane result `vop(x, e)` equals the root
/// operation `f` on that one lane with that one lane's own exponent. The
/// exponent vector is not uniform, so a form that broadcast a single exponent
/// across the whole vector cannot pass.
fn expectIntLanes(
    comptime I: type,
    comptime T: type,
    comptime Bits: type,
    comptime vop: anytype,
    comptime f: anytype,
    x: @Vector(lane_count, T),
    e: @Vector(lane_count, I),
) !void {
    const got: @Vector(lane_count, T) = vop(x, e);
    inline for (0..lane_count) |i| {
        const want: T = f(x[i], e[i]);
        try testing.expectEqual(@as(Bits, @bitCast(want)), @as(Bits, @bitCast(got[i])));
    }
}

/// Asserts lane `i` of the lane-ordered `PowerResult` array that `f(x, e)`
/// returns equals the root operation `f` on that one lane, field by field.
fn expectPowerBroadcast(
    comptime T: type,
    comptime f: anytype,
    x: @Vector(lane_count, T),
    e: i64,
) !void {
    const got: [lane_count]root.PowerResult(T) = f(x, e);
    inline for (0..lane_count) |i| {
        try expectLaneResult(root.PowerResult(T), got[i], f(x[i], e));
    }
}

/// Asserts lane `i` of the lane-ordered `PowerResult` array that the per-lane
/// entry point `vop(x, e)` returns equals the root operation `f` on that one
/// lane with that one lane's own exponent, field by field.
fn expectPowerLanes(
    comptime T: type,
    comptime vop: anytype,
    comptime f: anytype,
    x: @Vector(lane_count, T),
    e: @Vector(lane_count, i64),
) !void {
    const got: [lane_count]root.PowerResult(T) = vop(x, e);
    inline for (0..lane_count) |i| {
        try expectLaneResult(root.PowerResult(T), got[i], f(x[i], e[i]));
    }
}

/// Asserts lane `i` of `f(x, p)` equals the root operation `f` on that one
/// lane with the same extra argument, by raw bits.
fn expectExtraLanes(
    comptime T: type,
    comptime Bits: type,
    comptime f: anytype,
    x: @Vector(lane_count, T),
    p: anytype,
) !void {
    const got: @Vector(lane_count, T) = f(x, p);
    inline for (0..lane_count) |i| {
        const want: T = f(x[i], p);
        try testing.expectEqual(@as(Bits, @bitCast(want)), @as(Bits, @bitCast(got[i])));
    }
}

/// Asserts lane `i` of `f(x, p, q)` equals the root operation `f` on that one
/// lane with the same two extra arguments, by raw bits.
fn expectTwoExtrasLanes(
    comptime T: type,
    comptime Bits: type,
    comptime f: anytype,
    x: @Vector(lane_count, T),
    p: anytype,
    q: anytype,
) !void {
    const got: @Vector(lane_count, T) = f(x, p, q);
    inline for (0..lane_count) |i| {
        const want: T = f(x[i], p, q);
        try testing.expectEqual(@as(Bits, @bitCast(want)), @as(Bits, @bitCast(got[i])));
    }
}

test "the IEEE family float-to-float operations apply the root kernel lane for lane" {
    // Four unary float-to-float operations and the one ternary `fma`. Every
    // row carries a signed zero, a half-way tie, and a special value, so a lane
    // that dispatches into a kernel's special-case arm is compared there too.
    try checkRootUnary(
        root.round,
        .{ 0.5, -0.5, 2.5, -2.5 },
        .{ 0.5, -0.5, 2.5, std.math.inf(f32) },
        .{ 0.5, -0.5, 2.5, std.math.nan(f64) },
    );
    try checkRootUnary(
        root.nextUp,
        .{ 0.0, -0.0, 1.5, -1.5 },
        .{ 0.0, -0.0, 1.5, std.math.inf(f32) },
        .{ 0.0, -0.0, -1.5, std.math.nan(f64) },
    );
    try checkRootUnary(
        root.nextDown,
        .{ 0.0, -0.0, 1.5, -1.5 },
        .{ 0.0, -0.0, 1.5, -std.math.inf(f32) },
        .{ 0.0, -0.0, -1.5, std.math.nan(f64) },
    );
    try checkRootUnary(
        root.getPayload,
        .{ 1.5, 2.5, 3.5, 0.5 },
        .{ 1.5, 2.5, 3.5, 0.5 },
        .{ 1.5, 2.5, 3.5, 0.5 },
    );
    try checkRootTernary(
        root.fma,
        .{ 2.0, 3.0, -1.5, 1.0 },
        .{ 3.0, 0.5, 2.0, 1.0 },
        .{ 1.0, -1.0, 0.5, -1.0 },
        .{ 2.0, 3.0, -1.5, 1e30 },
        .{ 3.0, 0.5, 2.0, 1e30 },
        .{ 1.0, -1.0, 0.5, -1e30 },
        .{ 2.0, 3.0, -1.5, 1e300 },
        .{ 3.0, 0.5, 2.0, 1e300 },
        .{ 1.0, -1.0, 0.5, -1e300 },
    );
}

test "the IEEE family binary and extrema operations apply the root kernel lane for lane" {
    // `copySign` and `nextAfter`, then the six extrema aliases. Their root
    // names keep the underscores and their vector entry points do not, so
    // every row is pinned against the root name a caller actually writes.
    try checkRootBinary(
        root.copySign,
        .{ 3.0, -3.0, 0.0, 1.5 },
        .{ 1.0, -1.0, 0.0, -1.5 },
        .{ 3.0, -3.0, 0.0, 1.5 },
        .{ 1.0, -1.0, 0.0, -1.5 },
        .{ 3.0, -3.0, 0.0, std.math.nan(f64) },
        .{ 1.0, -1.0, 0.0, -std.math.inf(f64) },
    );
    try checkRootBinary(
        root.nextAfter,
        .{ 1.0, -1.0, 0.0, 2.0 },
        .{ 2.0, 2.0, -0.0, 1.0 },
        .{ 1.0, -1.0, 0.0, std.math.inf(f32) },
        .{ 2.0, 2.0, -0.0, 1.0 },
        .{ 1.0, -1.0, 0.0, std.math.nan(f64) },
        .{ 2.0, 2.0, -0.0, 1.0 },
    );
    try checkRootBinary(
        root.fmaximum_mag,
        .{ 1.0, -1.0, 2.0, 2.0 },
        .{ -1.0, 1.0, -2.0, 2.0 },
        .{ 1.0, -1.0, std.math.nan(f32), 2.0 },
        .{ -1.0, 1.0, 2.0, std.math.nan(f32) },
        .{ 1.0, -1.0, 2.0, std.math.inf(f64) },
        .{ -1.0, 1.0, 2.0, std.math.inf(f64) },
    );
    try checkRootBinary(
        root.fminimum_mag,
        .{ 1.0, -1.0, 2.0, 2.0 },
        .{ -1.0, 1.0, -2.0, 2.0 },
        .{ 1.0, -1.0, std.math.nan(f32), 2.0 },
        .{ -1.0, 1.0, 2.0, std.math.nan(f32) },
        .{ 1.0, -1.0, 2.0, std.math.inf(f64) },
        .{ -1.0, 1.0, 2.0, std.math.inf(f64) },
    );
    try checkRootBinary(
        root.fmaximum_num,
        .{ 1.0, -1.0, 2.0, 2.0 },
        .{ -1.0, 1.0, -2.0, 2.0 },
        .{ 1.0, -1.0, std.math.nan(f32), 2.0 },
        .{ -1.0, 1.0, 2.0, std.math.nan(f32) },
        .{ 1.0, -1.0, 2.0, std.math.inf(f64) },
        .{ -1.0, 1.0, 2.0, std.math.inf(f64) },
    );
    try checkRootBinary(
        root.fminimum_num,
        .{ 1.0, -1.0, 2.0, 2.0 },
        .{ -1.0, 1.0, -2.0, 2.0 },
        .{ 1.0, -1.0, std.math.nan(f32), 2.0 },
        .{ -1.0, 1.0, 2.0, std.math.nan(f32) },
        .{ 1.0, -1.0, 2.0, std.math.inf(f64) },
        .{ -1.0, 1.0, 2.0, std.math.inf(f64) },
    );
    try checkRootBinary(
        root.fmaximum_mag_num,
        .{ 1.0, -1.0, 2.0, 2.0 },
        .{ -1.0, 1.0, -2.0, 2.0 },
        .{ 1.0, -1.0, std.math.nan(f32), 2.0 },
        .{ -1.0, 1.0, 2.0, std.math.nan(f32) },
        .{ 1.0, -1.0, 2.0, std.math.inf(f64) },
        .{ -1.0, 1.0, 2.0, std.math.inf(f64) },
    );
    try checkRootBinary(
        root.fminimum_mag_num,
        .{ 1.0, -1.0, 2.0, 2.0 },
        .{ -1.0, 1.0, -2.0, 2.0 },
        .{ 1.0, -1.0, std.math.nan(f32), 2.0 },
        .{ -1.0, 1.0, 2.0, std.math.nan(f32) },
        .{ 1.0, -1.0, 2.0, std.math.inf(f64) },
        .{ -1.0, 1.0, 2.0, std.math.inf(f64) },
    );
}

test "the integer-returning family operations match the root exponents lane for lane" {
    const x16: @Vector(lane_count, f16) = .{ 1.0, 8.0, 0.5, std.math.inf(f16) };
    const x32: @Vector(lane_count, f32) = .{ 1.0, 8.0, 0.5, std.math.inf(f32) };
    const x64: @Vector(lane_count, f64) = .{ 1.0, 8.0, 0.5, std.math.inf(f64) };

    try checkRootIntUnary(i32, root.ilogb, x16, x32, x64);
    try checkRootIntUnary(i64, root.llogb, x16, x32, x64);

    // The exponents differ per lane, so the comparison above is not vacuous:
    // a kernel bound to the wrong format, or an entry point that broadcast
    // one lane's answer across the whole vector, would still be caught.
    const il: @Vector(lane_count, i32) = root.ilogb(x64);
    const ll: @Vector(lane_count, i64) = root.llogb(x64);
    try testing.expect(il[0] != il[1]);
    try testing.expect(ll[0] != ll[1]);
    try testing.expect(il[0] == 0);
    try testing.expect(il[1] == 3);
    try testing.expect(ll[0] == 0);
    try testing.expect(ll[1] == 3);
}

test "the boolean family operations return the root ordering mask lane for lane" {
    // Total order sorts the two zeros apart; totalOrderMag treats them as
    // equal, so each mask holds a true lane and a false lane. `expectMixedBoolLanes`
    // rejects a mask that is not mixed, which an always-true entry point could
    // never pass.
    try checkRootBoolBinary(
        root.totalOrder,
        .{ -0.0, 0.0, std.math.nan(f16), 1.0 },
        .{ 0.0, -0.0, 1.0, 1.0 },
        .{ -0.0, 0.0, std.math.nan(f32), 1.0 },
        .{ 0.0, -0.0, 1.0, 1.0 },
        .{ -0.0, 0.0, std.math.nan(f64), 1.0 },
        .{ 0.0, -0.0, 1.0, 1.0 },
    );
    try checkRootBoolBinary(
        root.totalOrderMag,
        .{ -0.0, 0.0, std.math.nan(f16), 1.0 },
        .{ 0.0, -0.0, 1.0, 1.0 },
        .{ -0.0, 0.0, std.math.nan(f32), 1.0 },
        .{ 0.0, -0.0, 1.0, 1.0 },
        .{ -0.0, 0.0, std.math.nan(f64), 1.0 },
        .{ 0.0, -0.0, 1.0, 1.0 },
    );
}

test "the result-struct family operations match the root fields lane for lane" {
    // The integer rounding family. Lane three is outside the binary64 integer
    // range, so the `ok` flag is false there and the flag is compared too.
    try checkRootUnaryFixedStruct(
        root.IntegerResult,
        root.lrint,
        .{ 2.5, -2.5, 0.0, 1.0 },
        .{ 2.5, -2.5, 0.0, 1.0 },
        .{ 2.5, -2.5, 0.0, 1e300 },
    );
    try checkRootUnaryFixedStruct(
        root.IntegerResult,
        root.llrint,
        .{ 2.5, -2.5, 0.0, 1.0 },
        .{ 2.5, -2.5, 0.0, 1.0 },
        .{ 2.5, -2.5, 0.0, 1e300 },
    );
    try checkRootUnaryFixedStruct(
        root.IntegerResult,
        root.lround,
        .{ 2.5, -2.5, 0.0, 1.0 },
        .{ 2.5, -2.5, 0.0, 1.0 },
        .{ 2.5, -2.5, 0.0, 1e300 },
    );
    try checkRootUnaryFixedStruct(
        root.IntegerResult,
        root.llround,
        .{ 2.5, -2.5, 0.0, 1.0 },
        .{ 2.5, -2.5, 0.0, 1.0 },
        .{ 2.5, -2.5, 0.0, 1e300 },
    );

    // A negative lane puts the `sign` of `lgamma` at both values, and lane
    // zero of `frexp` and `modf` puts both kernels on their signed zero arm.
    try checkRootUnaryStruct(
        root.lgamma,
        root.LgammaResult,
        .{ 0.5, 1.0, 6.0, -3.0 },
        .{ 0.5, 1.0, 6.0, -3.0 },
        .{ 0.5, 1.0, 6.0, -3.0 },
    );
    try checkRootUnaryStruct(
        root.frexp,
        root.FrexpResult,
        .{ 0.0, -0.0, 12.5, 0.75 },
        .{ 0.0, -0.0, 12.5, 0.75 },
        .{ 0.0, -0.0, 12.5, 0.75 },
    );
    try checkRootUnaryStruct(
        root.modf,
        root.ModfResult,
        .{ 3.75, -3.75, 0.0, -0.0 },
        .{ 3.75, -3.75, 0.0, -0.0 },
        .{ 3.75, -3.75, 0.0, -0.0 },
    );

    // The IEEE payload and canonicalization operations. `canonicalize` takes
    // a negative zero and an infinity, and `setPayload` and `setPayloadSig`
    // take payloads that are and are not admissible integers, so both the
    // quiet and the signalling `ok` outcomes are compared.
    try checkRootUnaryStruct(
        root.canonicalize,
        root.CanonicalizeResult,
        .{ 1.0, -0.0, 0.0, std.math.inf(f16) },
        .{ 1.0, -0.0, 0.0, std.math.inf(f32) },
        .{ 1.0, -0.0, 0.0, std.math.inf(f64) },
    );
    try checkRootUnaryStruct(
        root.setPayload,
        root.PayloadResult,
        .{ 1.0, 2.0, 1023.0, 2047.0 },
        .{ 1.0, 2.0, 1023.0, 2047.0 },
        .{ 1.0, 2.0, 1023.0, 2047.0 },
    );
    try checkRootUnaryStruct(
        root.setPayloadSig,
        root.PayloadResult,
        .{ 1.0, 2.0, 1023.0, 2047.0 },
        .{ 1.0, 2.0, 1023.0, 2047.0 },
        .{ 1.0, 2.0, 1023.0, 2047.0 },
    );

    // The binary result structs. The `atan2` lanes cover the signed zero
    // quadrants, and the `remquo` lanes cover a negative numerator and a
    // quotient that does not divide evenly.
    try checkRootBinaryStruct(
        root.atan2,
        root.Atan2Result,
        .{ 0.0, 1.0, -1.0, 2.0 },
        .{ -0.0, 0.0, 1.0, 2.0 },
        .{ 0.0, 1.0, -1.0, 2.0 },
        .{ -0.0, 0.0, 1.0, 2.0 },
        .{ 0.0, 1.0, -1.0, 2.0 },
        .{ -0.0, 0.0, 1.0, 2.0 },
    );
    try checkRootBinaryStruct(
        root.atan2pi,
        root.Atan2Result,
        .{ 0.0, 1.0, -1.0, 2.0 },
        .{ -0.0, 0.0, 1.0, 2.0 },
        .{ 0.0, 1.0, -1.0, 2.0 },
        .{ -0.0, 0.0, 1.0, 2.0 },
        .{ 0.0, 1.0, -1.0, 2.0 },
        .{ -0.0, 0.0, 1.0, 2.0 },
    );
    try checkRootBinaryStruct(
        root.remquo,
        root.RemquoResult,
        .{ 7.5, -7.5, 1.0, 10.25 },
        .{ 2.0, 2.0, 3.0, 4.0 },
        .{ 7.5, -7.5, 1.0, 10.25 },
        .{ 2.0, 2.0, 3.0, 4.0 },
        .{ 7.5, -7.5, 1.0, 10.25 },
        .{ 2.0, 2.0, 3.0, 4.0 },
    );
}

test "the narrowing family operations return binary32 lanes equal to the root operation" {
    // `root.fadd` and its five siblings declare their parameters as `f64`, so
    // the lane-parallel path is the family entry point the root dispatch line
    // names. The scalar reference is still the root operation itself, and the
    // result TYPE is asserted with `@TypeOf` rather than inferred from a
    // comparison that would also succeed for an `f64` result.
    const a: @Vector(lane_count, f64) = .{ 1.0, -2.5, 0.75, 1e300 };
    const b: @Vector(lane_count, f64) = .{ 3.0, 0.5, -1.25, 1e-300 };
    const c: @Vector(lane_count, f64) = .{ 0.25, 2.0, 4.0, -8.0 };

    const sq: @Vector(lane_count, f64) = .{ 4.0, 2.0, 0.25, 9.0 };
    try expectNarrowBinary(v.narrow.faddv, root.fadd, a, b);
    try expectNarrowBinary(v.narrow.fsubv, root.fsub, a, b);
    try expectNarrowBinary(v.narrow.fmulv, root.fmul, a, b);
    try expectNarrowBinary(v.narrow.fdivv, root.fdiv, a, b);
    try expectNarrowTernary(v.narrow.ffmav, root.ffma, a, b, c);
    try expectNarrowUnary(v.narrow.fsqrtv, root.fsqrt, sq);

    // The narrowing is real, not a widening that happens to compare equal: a
    // binary64 lane the binary32 result cannot hold exactly must round away.
    const ones: @Vector(lane_count, f64) = .{ 1.0, 1.0, 1.0, 1.0 };
    const threes: @Vector(lane_count, f64) = .{ 3.0, 3.0, 3.0, 3.0 };
    const thirds: @Vector(lane_count, f32) = v.narrow.fdivv(ones, threes);
    const exact: u32 = @bitCast(@as(f32, 1.0 / 3.0));
    try testing.expect(@as(u32, @bitCast(thirds[0])) == exact);
    try testing.expect(exact != @as(u32, @bitCast(@as(f32, 1.0 / 3.0 + 1.0 / 3.0 + 1.0 / 3.0))));
}

test "the integer-argument family operations pass their extra argument through" {
    // The broadcast forms. The exponent is one runtime scalar shared by every
    // lane, and the scalar reference is the root operation on one lane.
    const bx16: @Vector(lane_count, f16) = .{ 1.0, -2.5, 0.75, 2.0 };
    const bx32: @Vector(lane_count, f32) = .{ 1.0, -2.5, 0.75, 1e30 };
    const bx64: @Vector(lane_count, f64) = .{ 1.0, -2.5, 0.75, 1e30 };
    inline for (.{ f16, f32, f64 }) |T| {
        const Bits = @Int(.unsigned, @bitSizeOf(T));
        const x: @Vector(lane_count, T) = if (T == f16) bx16 else if (T == f32) bx32 else bx64;
        try expectBroadcastLanes(i32, T, Bits, root.ldexp, x, 3);
        try expectBroadcastLanes(i32, T, Bits, root.scalbn, x, 3);
        try expectBroadcastLanes(i64, T, Bits, root.scalbln, x, 3);
    }

    // The per-lane `Lanes` forms. The exponent vector is NOT uniform, so a
    // form that broadcast one lane's exponent across the whole vector, or that
    // read lane zero's exponent for every lane, cannot pass.
    const e32: @Vector(lane_count, i32) = .{ 0, 1, -2, 5 };
    const e64: @Vector(lane_count, i64) = .{ 0, 1, -2, 5 };
    inline for (.{ f16, f32, f64 }) |T| {
        const Bits = @Int(.unsigned, @bitSizeOf(T));
        const x: @Vector(lane_count, T) = if (T == f16) bx16 else if (T == f32) bx32 else bx64;
        try expectIntLanes(i32, T, Bits, v.intarg.ldexpvLanes, root.ldexp, x, e32);
        try expectIntLanes(i32, T, Bits, v.intarg.scalbnvLanes, root.scalbn, x, e32);
        try expectIntLanes(i64, T, Bits, v.intarg.scalblnvLanes, root.scalbln, x, e64);
    }

    // The three power operations, in both forms. `pown` uses 3 broadcast and
    // {3, 2, 0, 5} per lane, so lane two is a zeroth power and cannot match a
    // broadcast of any single exponent.
    const px16: @Vector(lane_count, f16) = .{ 2.0, -3.0, 0.5, 9.0 };
    const px32: @Vector(lane_count, f32) = .{ 2.0, -3.0, 0.5, 9.0 };
    const px64: @Vector(lane_count, f64) = .{ 2.0, -3.0, 0.5, 9.0 };
    const pe: @Vector(lane_count, i64) = .{ 3, 2, 0, 5 };
    const re: @Vector(lane_count, i64) = .{ 2, 3, 4, 1 };
    inline for (.{ f16, f32, f64 }) |T| {
        const x: @Vector(lane_count, T) = if (T == f16) px16 else if (T == f32) px32 else px64;
        try expectPowerBroadcast(T, root.pown, x, 3);
        try expectPowerBroadcast(T, root.compoundn, x, 3);
        try expectPowerBroadcast(T, root.rootn, x, 2);
        try expectPowerLanes(T, v.intarg.pownvLanes, root.pown, x, pe);
        try expectPowerLanes(T, v.intarg.compoundnvLanes, root.compoundn, x, pe);
        try expectPowerLanes(T, v.intarg.rootnvLanes, root.rootn, x, re);
    }

    // `nextToward` has a single `f64` target and no per-lane form, so it is
    // driven at two different targets: one that raises a lane and one that
    // lowers it, so a hardcoded target cannot pass.
    const nx16: @Vector(lane_count, f16) = .{ 1.0, -1.0, 0.0, 2.0 };
    const nx32: @Vector(lane_count, f32) = .{ 1.0, -1.0, 0.0, 2.0 };
    const nx64: @Vector(lane_count, f64) = .{ 1.0, -1.0, 0.0, 2.0 };
    inline for (.{ f16, f32, f64 }) |T| {
        const Bits = @Int(.unsigned, @bitSizeOf(T));
        const x: @Vector(lane_count, T) = if (T == f16) nx16 else if (T == f32) nx32 else nx64;
        try expectExtraLanes(T, Bits, root.nextToward, x, @as(f64, 2.0));
        try expectExtraLanes(T, Bits, root.nextToward, x, @as(f64, -2.0));
    }

    // The `fromfp` family, at two rounding directions and two widths, so a
    // hardcoded pair of extra arguments cannot pass.
    const sx16: @Vector(lane_count, f16) = .{ 2.5, -2.5, 7.75, -0.5 };
    const sx32: @Vector(lane_count, f32) = .{ 2.5, -2.5, 7.75, -0.5 };
    const sx64: @Vector(lane_count, f64) = .{ 2.5, -2.5, 7.75, -0.5 };
    // An unsigned conversion of a negative lane raises an invalid operation,
    // so the unsigned spellings only see non-negative lanes.
    const ux16: @Vector(lane_count, f16) = .{ 2.5, 7.75, 0.5, 0.0 };
    const ux32: @Vector(lane_count, f32) = .{ 2.5, 7.75, 0.5, 0.0 };
    const ux64: @Vector(lane_count, f64) = .{ 2.5, 7.75, 0.5, 0.0 };
    inline for (.{ f16, f32, f64 }) |T| {
        const Bits = @Int(.unsigned, @bitSizeOf(T));
        const signed_lanes: @Vector(lane_count, T) =
            if (T == f16) sx16 else if (T == f32) sx32 else sx64;
        const unsigned_lanes: @Vector(lane_count, T) =
            if (T == f16) ux16 else if (T == f32) ux32 else ux64;
        inline for (.{ root.FromFpRound.nearest_even, root.FromFpRound.toward_zero }) |direction| {
            inline for (.{ 32, 64 }) |width| {
                const w: u32 = width;
                try expectTwoExtrasLanes(T, Bits, root.fromfp, signed_lanes, direction, w);
                try expectTwoExtrasLanes(T, Bits, root.fromfpx, signed_lanes, direction, w);
                try expectTwoExtrasLanes(T, Bits, root.ufromfp, unsigned_lanes, direction, w);
                try expectTwoExtrasLanes(T, Bits, root.ufromfpx, unsigned_lanes, direction, w);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Completeness: this guard covers the WHOLE vectorized root surface
//
// The tables above name the operations they pin. Without a cross-check, a new
// `@Vector` dispatch line could be added to `root.zig` and no table would
// follow it, which is exactly the gap this file exists to close. The check
// below runs in both directions:
//
// - Every operation the tables name must be a real `root` declaration, and
//   the lane-parallel entry point the table claims must really be declared in
//   the family module the root dispatch line names.
// - Every `root` declaration whose value is a function must be either pinned
//   by a table or named in `scalar_only_ops` below.
// - The two sets must be the same SIZE, so a duplicate in the tables cannot
//   stand in for a real operation.
//
// The root declaration list comes from `std.meta.declarations`, so it cannot
// drift from the file. A declaration is an operation exactly when its value is
// a function; the module and result-type re-exports are values of a container
// type and are skipped by that rule alone, with no hand-maintained exclusion
// list to fall out of date.
// ---------------------------------------------------------------------------

/// The module that holds an operation's lane-parallel entry point.
const Family = enum {
    /// The float-to-float surface declared directly in `vector.zig`, whose
    /// entry point for `OP` is `OPv` by that module's own naming rule.
    surface,
    struct_family,
    intarg,
    narrow,
    ieee_family,
};

fn FamilyModule(comptime f: Family) type {
    return switch (f) {
        .surface => vector,
        .struct_family => vector.struct_family,
        .intarg => vector.intarg,
        .narrow => vector.narrow,
        .ieee_family => vector.ieee_family,
    };
}

/// The 55 float-to-float operations the first table pins, in the order that
/// table lists them: 45 unary, then 10 binary.
const surface_ops = [_][:0]const u8{
    "abs",       "acos",     "acosh",     "acospi", "asin",      "asinh",
    "asinpi",    "atan",     "atanh",     "atanpi", "cbrt",      "ceil",
    "cos",       "cosh",     "cospi",     "erf",    "erfc",      "exp",
    "exp10",     "exp10m1",  "exp2",      "exp2m1", "expm1",     "floor",
    "log",       "log10",    "log10p1",   "log1p",  "log2",      "log2p1",
    "logb",      "logp1",    "nearbyint", "rint",   "roundeven", "rsqrt",
    "sin",       "sinh",     "sinpi",     "sqrt",   "tan",       "tanh",
    "tanpi",     "tgamma",   "trunc",     "fdim",   "fmax",      "fmaximum",
    "fmin",      "fminimum", "fmod",      "hypot",  "pow",       "powr",
    "remainder",
};

/// One operation pinned by a table above, with the family module that holds
/// its entry point and the entry point's own name. The six extrema aliases
/// carry a snake_case root name and a camelCase entry point, so the entry point
/// is spelled out here instead of being derived from the root name.
const FamilyOp = struct {
    op: [:0]const u8,
    family: Family,
    entry: [:0]const u8,
};

/// The thirteen operations whose scalar kernel returns a result struct. They
/// cannot return a `@Vector`, so each lane is a lane-ordered array element.
const struct_family_ops = [_]FamilyOp{
    .{ .op = "lrint", .family = .struct_family, .entry = "lrintv" },
    .{ .op = "llrint", .family = .struct_family, .entry = "llrintv" },
    .{ .op = "lround", .family = .struct_family, .entry = "lroundv" },
    .{ .op = "llround", .family = .struct_family, .entry = "llroundv" },
    .{ .op = "lgamma", .family = .struct_family, .entry = "lgammav" },
    .{ .op = "frexp", .family = .struct_family, .entry = "frexpv" },
    .{ .op = "modf", .family = .struct_family, .entry = "modfv" },
    .{ .op = "canonicalize", .family = .struct_family, .entry = "canonicalizev" },
    .{ .op = "setPayload", .family = .struct_family, .entry = "setPayloadv" },
    .{ .op = "setPayloadSig", .family = .struct_family, .entry = "setPayloadSigv" },
    .{ .op = "atan2", .family = .struct_family, .entry = "atan2v" },
    .{ .op = "atan2pi", .family = .struct_family, .entry = "atan2piv" },
    .{ .op = "remquo", .family = .struct_family, .entry = "remquov" },
};

/// The eleven operations that take a non-float extra argument.
const intarg_ops = [_]FamilyOp{
    .{ .op = "ldexp", .family = .intarg, .entry = "ldexpv" },
    .{ .op = "scalbn", .family = .intarg, .entry = "scalbnv" },
    .{ .op = "scalbln", .family = .intarg, .entry = "scalblnv" },
    .{ .op = "pown", .family = .intarg, .entry = "pownv" },
    .{ .op = "compoundn", .family = .intarg, .entry = "compoundnv" },
    .{ .op = "rootn", .family = .intarg, .entry = "rootnv" },
    .{ .op = "nextToward", .family = .intarg, .entry = "nextTowardv" },
    .{ .op = "fromfp", .family = .intarg, .entry = "fromfpv" },
    .{ .op = "ufromfp", .family = .intarg, .entry = "ufromfpv" },
    .{ .op = "fromfpx", .family = .intarg, .entry = "fromfpxv" },
    .{ .op = "ufromfpx", .family = .intarg, .entry = "ufromfpxv" },
};

/// The ternary `fma` and the six narrowing operations, which take `f64` lanes
/// and, for the narrowing six, produce `f32` lanes.
const narrow_ops = [_]FamilyOp{
    .{ .op = "fma", .family = .narrow, .entry = "fmav" },
    .{ .op = "fadd", .family = .narrow, .entry = "faddv" },
    .{ .op = "fsub", .family = .narrow, .entry = "fsubv" },
    .{ .op = "fmul", .family = .narrow, .entry = "fmulv" },
    .{ .op = "fdiv", .family = .narrow, .entry = "fdivv" },
    .{ .op = "ffma", .family = .narrow, .entry = "ffmav" },
    .{ .op = "fsqrt", .family = .narrow, .entry = "fsqrtv" },
};

/// The remaining IEEE and miscellaneous operations, including the six extrema
/// aliases, the integer-returning pair, and the boolean-returning pair.
const ieee_family_ops = [_]FamilyOp{
    .{ .op = "round", .family = .ieee_family, .entry = "roundv" },
    .{ .op = "nextUp", .family = .ieee_family, .entry = "nextUpv" },
    .{ .op = "nextDown", .family = .ieee_family, .entry = "nextDownv" },
    .{ .op = "getPayload", .family = .ieee_family, .entry = "getPayloadv" },
    .{ .op = "copySign", .family = .ieee_family, .entry = "copySignv" },
    .{ .op = "nextAfter", .family = .ieee_family, .entry = "nextAfterv" },
    .{ .op = "ilogb", .family = .ieee_family, .entry = "ilogbv" },
    .{ .op = "llogb", .family = .ieee_family, .entry = "llogbv" },
    .{ .op = "totalOrder", .family = .ieee_family, .entry = "totalOrderv" },
    .{ .op = "totalOrderMag", .family = .ieee_family, .entry = "totalOrderMagv" },
    .{ .op = "fmaximum_mag", .family = .ieee_family, .entry = "fmaximumMagv" },
    .{ .op = "fminimum_mag", .family = .ieee_family, .entry = "fminimumMagv" },
    .{ .op = "fmaximum_num", .family = .ieee_family, .entry = "fmaximumNumv" },
    .{ .op = "fminimum_num", .family = .ieee_family, .entry = "fminimumNumv" },
    .{ .op = "fmaximum_mag_num", .family = .ieee_family, .entry = "fmaximumMagNumv" },
    .{ .op = "fminimum_mag_num", .family = .ieee_family, .entry = "fminimumMagNumv" },
};

/// Every operation the tables above pin outside the float-to-float surface.
const family_ops = struct_family_ops ++ intarg_ops ++ narrow_ops ++ ieee_family_ops;

/// The one root operation that deliberately carries no vector dispatch line.
/// `nan` takes a comptime type and a tag slice, so it has no lane to map over.
const scalar_only_ops = [_][:0]const u8{"nan"};

/// Reports whether `name` is in the comptime-known set `set`.
///
/// The set is walked by index, never through a `for` capture, so the walk has
/// no capture state to carry between iterations.
fn memberOf(comptime name: [:0]const u8, comptime set: anytype) bool {
    for (set, 0..) |_, i| {
        if (std.mem.eql(u8, set[i], name)) return true;
    }
    return false;
}

/// Reports whether the family tables pin the root operation `name`.
fn familyMemberOf(comptime name: [:0]const u8) bool {
    for (family_ops, 0..) |_, i| {
        if (std.mem.eql(u8, family_ops[i].op, name)) return true;
    }
    return false;
}

/// Reports whether the root declaration `name` is an operation, which is true
/// exactly when it is a function that does not return a `type`.
///
/// The rule is structural, so it cannot drift from `root.zig`. It excludes the
/// three container re-exports `ieee`, `scalar`, and `vector` and the ten
/// result-type re-exports, because those are values of a container type rather
/// than functions. It also excludes the eight result-type CONSTRUCTORS
/// `Atan2Result`, `PowerResult`, and the rest, because those are the only
/// functions in `root.zig` that return a `type`. A `null` return type means
/// the return type is derived from an `anytype` parameter, which is how almost
/// every operation is written, and counts as an operation.
fn isRootOperation(comptime name: [:0]const u8) bool {
    const info = @typeInfo(@TypeOf(@field(root, name)));
    if (info != .@"fn") return false;
    const return_type = info.@"fn".return_type orelse return true;
    return return_type != type;
}

/// The outcome of classifying every public `root` declaration against the
/// tables above. An empty name means that direction found nothing wrong.
const Audit = struct {
    /// How many public `root` declarations are operations.
    operations: usize,
    /// How many of those the tables above pin.
    covered: usize,
    /// The first table row naming a declaration `root.zig` does not export.
    bad_operation: [:0]const u8,
    /// The first table row whose lane-parallel entry point is not declared in
    /// the module its `root.zig` dispatch line names.
    missing_entry: [:0]const u8,
    /// The first `root` operation that carries a vector dispatch line and has
    /// no row in this guard.
    uncovered: [:0]const u8,
    /// The first declared scalar-only exception that a table pins anyway.
    misfiled: [:0]const u8,
};

/// Classifies every public `root` declaration against the tables above, in
/// both directions, entirely at comptime.
///
/// `isRootOperation` decides what counts as an operation, so the sweep has no
/// hand-maintained list of the names to skip and the classification cannot
/// drift from `root.zig`.
fn audit() Audit {
    // Walking 102 table rows, formatting their names, and then sweeping every
    // root declaration is a deliberately large comptime evaluation.
    @setEvalBranchQuota(1_000_000);

    var result: Audit = .{
        .operations = 0,
        .covered = 0,
        .bad_operation = "",
        .missing_entry = "",
        .uncovered = "",
        .misfiled = "",
    };

    // Direction one: every row names a real root operation, and the
    // lane-parallel entry point it claims is declared in the module its
    // dispatch line names. `vector.zig` derives a float-to-float entry point
    // from the root name, so that name is derived here too.
    for (surface_ops, 0..) |_, i| {
        const op = surface_ops[i];
        if (!@hasDecl(root, op)) {
            if (result.bad_operation.len == 0) result.bad_operation = op;
        } else if (!@hasDecl(vector, std.fmt.comptimePrint("{s}v", .{op}))) {
            if (result.missing_entry.len == 0) result.missing_entry = op;
        }
    }
    for (family_ops, 0..) |_, i| {
        const g = family_ops[i];
        if (!@hasDecl(root, g.op)) {
            if (result.bad_operation.len == 0) result.bad_operation = g.op;
        } else if (!@hasDecl(FamilyModule(g.family), g.entry)) {
            if (result.missing_entry.len == 0) result.missing_entry = g.op;
        }
    }
    for (scalar_only_ops, 0..) |_, i| {
        const name = scalar_only_ops[i];
        if (!@hasDecl(root, name)) {
            if (result.bad_operation.len == 0) result.bad_operation = name;
        }
        if (memberOf(name, surface_ops) or familyMemberOf(name)) {
            if (result.misfiled.len == 0) result.misfiled = name;
        }
    }

    // Direction two: nothing in `root.zig` escapes the tables. The declared
    // exceptions are subtracted from the sweep, so a root operation given a
    // `@Vector` dispatch line without a table row is reported here.
    for (std.meta.declarations(root), 0..) |_, i| {
        const name = std.meta.declarations(root)[i].name;
        if (!isRootOperation(name)) continue;
        result.operations += 1;
        if (memberOf(name, scalar_only_ops)) continue;
        if (memberOf(name, surface_ops) or familyMemberOf(name)) {
            result.covered += 1;
        } else if (result.uncovered.len == 0) {
            result.uncovered = name;
        }
    }
    return result;
}

test "the guard pins every vectorized root operation and nothing else" {
    const report = comptime audit();

    if (report.bad_operation.len != 0) {
        @compileError(std.fmt.comptimePrint(
            "a guard row names `{s}`, which `root.zig` does not export",
            .{report.bad_operation},
        ));
    }
    if (report.missing_entry.len != 0) {
        @compileError(std.fmt.comptimePrint(
            "the guard row for `{s}` names an entry point its family module lacks",
            .{report.missing_entry},
        ));
    }
    if (report.uncovered.len != 0) {
        @compileError(std.fmt.comptimePrint(
            "root operation `{s}` has a vector dispatch line and no row here",
            .{report.uncovered},
        ));
    }
    if (report.misfiled.len != 0) {
        @compileError(std.fmt.comptimePrint(
            "`{s}` is declared scalar-only but a table pins it",
            .{report.misfiled},
        ));
    }

    // The two directions meet: the tables name exactly as many distinct root
    // operations as the sweep found, so a duplicated row cannot stand in for a
    // missing one.
    try testing.expectEqual(report.covered, surface_ops.len + family_ops.len);
    try testing.expectEqual(report.operations - scalar_only_ops.len, report.covered);

    // The family split, so no table can be silently truncated.
    try testing.expectEqual(@as(usize, 55), surface_ops.len);
    try testing.expectEqual(@as(usize, 13), struct_family_ops.len);
    try testing.expectEqual(@as(usize, 11), intarg_ops.len);
    try testing.expectEqual(@as(usize, 7), narrow_ops.len);
    try testing.expectEqual(@as(usize, 16), ieee_family_ops.len);
    try testing.expectEqual(@as(usize, 1), scalar_only_ops.len);
}
