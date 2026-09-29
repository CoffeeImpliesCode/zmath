//! Guard over the `@Vector` surface in `vector.zig`.
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
//! 3. Operations that return a struct or an integer have no float vector entry
//!    point, because this surface is exactly `@Vector(N, T) -> @Vector(N, T)`.

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
    // The float vector surface is exactly `@Vector(N, T) -> @Vector(N, T)`.
    // Operations that return a struct (`Atan2Result`, `PowerResult`,
    // `LgammaResult`, `RemquoResult`, `FrexpResult`, `ModfResult`,
    // `IntegerResult`, `PayloadResult`, `CanonicalizeResult`), a bare integer
    // (`ilogb`, `llogb`), or a `bool` (`totalOrder`, `totalOrderMag`) would
    // need a `@Vector(N, Result)` return type, which this module does not
    // provide. Vector support for them is intentionally deferred, so their `v`
    // names must stay absent.
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
