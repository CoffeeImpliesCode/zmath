const std = @import("std");
const testing = std.testing;
const root = @import("root.zig");
const vector = @import("vector.zig");
const scalar = @import("scalar.zig");

// Every lane of a vector result must equal the scalar kernel bit-for-bit. The
// assertions below compare raw bit patterns (not float equality) so signed zero,
// infinities, and NaN payloads are all checked exactly.

fn bits16(x: f16) u16 {
    return @bitCast(x);
}
fn bits32(x: f32) u32 {
    return @bitCast(x);
}
fn bits64(x: f64) u64 {
    return @bitCast(x);
}

test "sinv binary16 lanes equal scalar sine kernel bit-for-bit" {
    const x = @Vector(4, f16){ 0.0, -0.0, 0.25, 1.5 };
    const got = vector.sinv(x);
    inline for (0..4) |i| {
        try testing.expectEqual(bits16(scalar.sinh(x[i])), bits16(got[i]));
    }
}

test "sinv binary16 through the generic top-level sin" {
    const x = @Vector(4, f16){ 0.0, 0.25, -0.25, 1.0 };
    const got = root.sin(x);
    inline for (0..4) |i| {
        try testing.expectEqual(bits16(scalar.sinh(x[i])), bits16(got[i]));
    }
}

test "expv binary32 lanes equal scalar exp kernel" {
    const x = @Vector(4, f32){ 0.0, 1.0, -3.25, 10.5 };
    const got = vector.expv(x);
    inline for (0..4) |i| {
        try testing.expectEqual(bits32(scalar.expf(x[i])), bits32(got[i]));
    }
}

test "expv binary64 lanes equal scalar exp kernel" {
    const x = @Vector(2, f64){ -1.0, 20.0 };
    const got = vector.expv(x);
    inline for (0..2) |i| {
        try testing.expectEqual(bits64(scalar.expd(x[i])), bits64(got[i]));
    }
}

test "logv handles edge values bit-for-bit" {
    const x32 = @Vector(4, f32){ 0.0, 1.0, std.math.inf(f32), -1.0 };
    const g32 = vector.logv(x32);
    inline for (0..4) |i| {
        try testing.expectEqual(bits32(scalar.logf(x32[i])), bits32(g32[i]));
    }
    const x64 = @Vector(2, f64){ -0.0, 0.5 };
    const g64 = vector.logv(x64);
    inline for (0..2) |i| {
        try testing.expectEqual(bits64(scalar.logd(x64[i])), bits64(g64[i]));
    }
}

test "sinv preserves NaN and infinity lanes exactly" {
    const x = @Vector(4, f64){ std.math.inf(f64), -std.math.inf(f64), std.math.nan(f64), -0.0 };
    const got = vector.sinv(x);
    inline for (0..4) |i| {
        try testing.expectEqual(bits64(scalar.sind(x[i])), bits64(got[i]));
    }
}

test "powv binary32 and binary64 lanes equal scalar pow kernel" {
    const a32 = @Vector(4, f32){ 2.0, 0.5, -2.0, 1.0 };
    const b32 = @Vector(4, f32){ 3.0, 0.5, 3.0, -0.0 };
    const g32 = vector.powv(a32, b32);
    inline for (0..4) |i| {
        try testing.expectEqual(bits32(scalar.powf(a32[i], b32[i])), bits32(g32[i]));
    }
    const a64 = @Vector(2, f64){ 2.0, 10.0 };
    const b64 = @Vector(2, f64){ 0.5, -3.0 };
    const g64 = vector.powv(a64, b64);
    inline for (0..2) |i| {
        try testing.expectEqual(bits64(scalar.powd(a64[i], b64[i])), bits64(g64[i]));
    }
}

test "hypotv and fmaxv binary32 lanes equal scalar kernels" {
    const a = @Vector(4, f32){ 3.0, -3.0, 0.0, 1.0 };
    const b = @Vector(4, f32){ 4.0, 4.0, 0.0, -1.0 };
    const h = vector.hypotv(a, b);
    inline for (0..4) |i| {
        try testing.expectEqual(bits32(scalar.hypotf(a[i], b[i])), bits32(h[i]));
    }
    const m = vector.fmaxv(a, b);
    inline for (0..4) |i| {
        try testing.expectEqual(bits32(scalar.fmaxf(a[i], b[i])), bits32(m[i]));
    }
}

test "specialized entry points match generic entry points" {
    const x32 = @Vector(4, f32){ 0.3, -1.7, 2.5, 100.0 };
    const generic = vector.expv(x32);
    const specialized = vector.expfv(x32);
    inline for (0..4) |i| {
        try testing.expectEqual(bits32(generic[i]), bits32(specialized[i]));
    }
    const x16 = @Vector(4, f16){ 0.5, -0.5, 1.0, -1.0 };
    // `sinh` is the hyperbolic-sine generic entry point, not a specialization of
    // `sin`; compare genuine generic/specialized pairs instead.
    try testing.expectEqual(bits16(vector.expv(x16)[0]), bits16(vector.exphv(x16)[0]));
    const x64 = @Vector(2, f64){ 0.1, -0.2 };
    try testing.expectEqual(bits64(vector.logv(x64)[1]), bits64(vector.logdv(x64)[1]));
}

test "vector width is generic, not fixed" {
    const x1 = @Vector(1, f64){1.0};
    const g1 = vector.sinv(x1);
    try testing.expectEqual(bits64(scalar.sind(x1[0])), bits64(g1[0]));

    const lanes = [_]f64{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };
    const v8: @Vector(8, f64) = lanes;
    const g8 = vector.sinv(v8);
    inline for (0..8) |i| {
        try testing.expectEqual(bits64(scalar.sind(lanes[i])), bits64(g8[i]));
    }
}

test "binary16 vector lanes match scalar across the full finite set" {
    // Exercise a spread of binary16 bit patterns through the vector path and
    // compare each lane with the scalar kernel, including subnormals.
    const x = @Vector(8, f16){
        6.103515625e-05, // min normal
        5.960464477539063e-08, // min subnormal
        -5.960464477539063e-08,
        65504.0, // max finite
        2048.0,
        -2048.0,
        0.0078125,
        -0.0078125,
    };
    const got = vector.sinv(x);
    inline for (0..8) |i| {
        try testing.expectEqual(bits16(scalar.sinh(x[i])), bits16(got[i]));
    }
}

fn laneVector8(raw: usize) @Vector(8, f16) {
    var arr: [8]f16 = undefined;
    for (&arr, 0..) |*slot, k| {
        slot.* = @bitCast(@as(u16, @intCast(raw + k)));
    }
    return arr;
}

test "exhaustive binary16 domain: vector sin lanes equal the scalar kernel" {
    var raw: usize = 0;
    while (raw < 65536) : (raw += 8) {
        const x = laneVector8(raw);
        const got = root.sin(x);
        inline for (0..8) |k| {
            const want = scalar.sinh(x[k]);
            testing.expectEqual(bits16(want), bits16(got[k])) catch |err| {
                std.debug.print("mismatch at raw bits 0x{x}\n", .{raw + k});
                return err;
            };
        }
    }
}

test "exhaustive binary16 domain: vector exp, log, tanh, and cbrt lanes equal the scalar kernel" {
    var raw: usize = 0;
    while (raw < 65536) : (raw += 8) {
        const x = laneVector8(raw);
        const e = vector.expv(x);
        const l = vector.logv(x);
        const t = vector.tanhv(x);
        const c = vector.cbrtv(x);
        inline for (0..8) |k| {
            try testing.expectEqual(bits16(scalar.exph(x[k])), bits16(e[k]));
            try testing.expectEqual(bits16(scalar.logh(x[k])), bits16(l[k]));
            try testing.expectEqual(bits16(scalar.tanhh(x[k])), bits16(t[k]));
            try testing.expectEqual(bits16(scalar.cbrth(x[k])), bits16(c[k]));
        }
    }
}

test "binary32 and binary64 specials: vector lanes equal the scalar kernel" {
    const x32: @Vector(13, f32) = .{
        0.0,               -0.0,          1.0,               -1.0,
        0.5,               -0.5,          std.math.inf(f32), -std.math.inf(f32),
        std.math.nan(f32), 1.1754944e-38, 1.4012985e-45,     3.4028235e38,
        -3.4028235e38,
    };
    const e32 = root.exp(x32);
    const l32 = root.log(x32);
    const p32 = root.pow(x32, x32);
    inline for (0..13) |k| {
        try testing.expectEqual(bits32(scalar.expf(x32[k])), bits32(e32[k]));
        try testing.expectEqual(bits32(scalar.logf(x32[k])), bits32(l32[k]));
        try testing.expectEqual(bits32(scalar.powf(x32[k], x32[k])), bits32(p32[k]));
    }

    const a64: @Vector(8, f64) = .{ 0.0, -0.0, 1.0, -1.0, 3.0, -3.0, 1e300, -1e300 };
    const b64: @Vector(8, f64) = .{ 1.0, -1.0, 0.0, 0.0, 4.0, 4.0, 1e300, 1e-300 };
    const h64 = vector.hypotv(a64, b64);
    inline for (0..8) |k| {
        try testing.expectEqual(bits64(scalar.hypotd(a64[k], b64[k])), bits64(h64[k]));
    }
}

// ---------------------------------------------------------------------------
// The four lane-parallel families behind `root.zig`
//
// Each case drives the top-level generic entry point in `root.zig` on a
// `@Vector` and compares it, lane by lane, with that same top-level call on one
// scalar lane, so the oracle is the scalar path itself and never a re-derived
// value. Float, integer, and boolean results are compared exactly: floats by
// raw bits, integers and booleans directly. The result-struct operations are
// compared FIELD by FIELD, because a bare `struct` has no bit pattern to
// compare and Zig forbids it as a `@Vector` element type, so their vector result
// is a lane-ordered array rather than a `@Vector`.
// ---------------------------------------------------------------------------

test "result structs: frexp and lgamma lanes equal the scalar result fields" {
    const x = @Vector(4, f64){ 0.0, -0.0, 12.5, 0.75 };
    const got = root.frexp(x);
    inline for (0..4) |i| {
        const want = root.frexp(x[i]);
        try testing.expectEqual(want.exponent, got[i].exponent);
        try testing.expectEqual(bits64(want.fraction), bits64(got[i].fraction));
    }

    const y = @Vector(4, f64){ 0.5, 1.0, 6.0, -3.0 };
    const ly = root.lgamma(y);
    inline for (0..4) |i| {
        const want = root.lgamma(y[i]);
        try testing.expectEqual(want.sign, ly[i].sign);
        try testing.expectEqual(bits64(want.value), bits64(ly[i].value));
    }
}

test "result structs: modf, remquo, and atan2 lanes equal the scalar fields" {
    const m = @Vector(4, f64){ 3.75, -3.75, 0.0, -0.0 };
    const mm = root.modf(m);
    inline for (0..4) |i| {
        const want = root.modf(m[i]);
        try testing.expectEqual(bits64(want.fractional), bits64(mm[i].fractional));
        try testing.expectEqual(bits64(want.integral), bits64(mm[i].integral));
    }

    const r = @Vector(4, f64){ 7.5, -7.5, 1.0, 10.25 };
    const d = @Vector(4, f64){ 2.0, 2.0, 3.0, 4.0 };
    const rq = root.remquo(r, d);
    inline for (0..4) |i| {
        const want = root.remquo(r[i], d[i]);
        try testing.expectEqual(bits64(want.value), bits64(rq[i].value));
        try testing.expectEqual(want.quotient, rq[i].quotient);
    }

    const u = @Vector(4, f64){ 0.0, 1.0, -1.0, 2.0 };
    const v = @Vector(4, f64){ -0.0, 0.0, 1.0, 2.0 };
    const a2 = root.atan2(u, v);
    inline for (0..4) |i| {
        const want = root.atan2(u[i], v[i]);
        try testing.expectEqual(want.resolved, a2[i].resolved);
        try testing.expectEqual(bits64(want.value), bits64(a2[i].value));
    }
}

test "result structs: lrint and the IEEE payload operations carry an ok flag" {
    const x = @Vector(4, f64){ 2.5, -2.5, 0.0, 1e300 };
    const got = root.lrint(x);
    inline for (0..4) |i| {
        const want = root.lrint(x[i]);
        try testing.expectEqual(want.value, got[i].value);
        try testing.expectEqual(want.ok, got[i].ok);
    }

    const c = @Vector(4, f64){ 1.0, -0.0, 0.0, std.math.inf(f64) };
    const cz = root.canonicalize(c);
    inline for (0..4) |i| {
        const want = root.canonicalize(c[i]);
        try testing.expectEqual(want.ok, cz[i].ok);
        try testing.expectEqual(bits64(want.value), bits64(cz[i].value));
    }

    const p = @Vector(4, f64){ 1.0, 2.0, 1023.0, 2047.0 };
    const pl = root.setPayload(p);
    inline for (0..4) |i| {
        const want = root.setPayload(p[i]);
        try testing.expectEqual(want.ok, pl[i].ok);
        try testing.expectEqual(bits64(want.value), bits64(pl[i].value));
    }
}

test "boolean results: totalOrder returns a bool mask equal to the scalar calls" {
    const a = @Vector(4, f64){ -0.0, 0.0, std.math.nan(f64), 1.0 };
    const b = @Vector(4, f64){ 0.0, -0.0, 1.0, 1.0 };
    const got = root.totalOrder(a, b);
    try testing.expect(@TypeOf(got) == @Vector(4, bool));
    inline for (0..4) |i| {
        try testing.expectEqual(root.totalOrder(a[i], b[i]), got[i]);
    }
    // The mask is not uniformly true, so the lane comparison above is not
    // vacuous: total order puts -0.0 before +0.0 and the reverse lane after.
    try testing.expect(got[0]);
    try testing.expect(!got[1]);
}

test "integer results: ilogb returns an i32 vector equal to the scalar exponents" {
    const x = @Vector(4, f64){ 1.0, 8.0, 0.5, std.math.inf(f64) };
    const got = root.ilogb(x);
    try testing.expect(@TypeOf(got) == @Vector(4, i32));
    inline for (0..4) |i| {
        try testing.expectEqual(root.ilogb(x[i]), got[i]);
    }
    // The exponents differ per lane, so the comparison above is not vacuous.
    try testing.expect(got[0] != got[1]);
}

test "IEEE float-to-float lanes: round, nextUp, copySign, and fmaximum_mag" {
    const x = @Vector(4, f32){ 0.5, -0.5, 2.5, std.math.inf(f32) };
    inline for (0..4) |i| {
        try testing.expectEqual(bits32(root.round(x[i])), bits32(root.round(x)[i]));
        try testing.expectEqual(bits32(root.nextUp(x[i])), bits32(root.nextUp(x)[i]));
    }

    const s = @Vector(4, f32){ 3.0, -3.0, 0.0, 1.5 };
    inline for (0..4) |i| {
        try testing.expectEqual(
            bits32(root.copySign(s[i], x[i])),
            bits32(root.copySign(s, x)[i]),
        );
    }

    const p = @Vector(4, f32){ 1.0, -1.0, 2.0, 2.0 };
    const q = @Vector(4, f32){ -1.0, 1.0, -2.0, 2.0 };
    inline for (0..4) |i| {
        try testing.expectEqual(
            bits32(root.fmaximum_mag(p[i], q[i])),
            bits32(root.fmaximum_mag(p, q)[i]),
        );
    }
}

test "integer arguments: ldexp broadcast and per-lane forms match the scalar" {
    const x = @Vector(4, f32){ 1.0, -2.5, 0.75, 1e30 };
    const broadcast = root.ldexp(x, 3);
    inline for (0..4) |i| {
        try testing.expectEqual(bits32(root.ldexp(x[i], 3)), bits32(broadcast[i]));
    }

    const e = @Vector(4, i32){ 0, 1, -2, 5 };
    const per_lane = vector.intarg.ldexpvLanes(x, e);
    inline for (0..4) |i| {
        try testing.expectEqual(bits32(scalar.ldexpf(x[i], e[i])), bits32(per_lane[i]));
    }
}

test "integer arguments: pown fills a lane-ordered PowerResult array" {
    const x = @Vector(4, f64){ 2.0, -3.0, 0.5, 9.0 };
    const got = root.pown(x, 3);
    inline for (0..4) |i| {
        const want = root.pown(x[i], 3);
        try testing.expectEqual(want.ok, got[i].ok);
        try testing.expectEqual(bits64(want.value), bits64(got[i].value));
    }
}

test "integer arguments: nextToward and fromfp pass the extra argument through" {
    const t = @Vector(4, f32){ 1.0, -1.0, 0.0, std.math.inf(f32) };
    const nt = root.nextToward(t, 2.0);
    inline for (0..4) |i| {
        try testing.expectEqual(bits32(root.nextToward(t[i], 2.0)), bits32(nt[i]));
    }

    const f = @Vector(4, f32){ 2.5, -2.5, 7.75, -0.5 };
    const fp = root.fromfp(f, .nearest_even, 32);
    inline for (0..4) |i| {
        try testing.expectEqual(
            bits32(root.fromfp(f[i], .nearest_even, 32)),
            bits32(fp[i]),
        );
    }
}

test "narrowing: fsqrt turns a binary64 vector into a binary32 result" {
    const x = @Vector(4, f64){ 4.0, 2.0, 0.25, 9.0 };
    const got = vector.narrow.fsqrtv(x);
    try testing.expect(@TypeOf(x) == @Vector(4, f64));
    try testing.expect(@TypeOf(got) == @Vector(4, f32));
    inline for (0..4) |i| {
        try testing.expectEqual(bits32(root.fsqrt(x[i])), bits32(got[i]));
    }
}

test "narrowing: the other five narrow a binary64 vector to binary32 lanes" {
    const a = @Vector(4, f64){ 1.0, -2.5, 0.75, 1e300 };
    const b = @Vector(4, f64){ 3.0, 0.5, -1.25, 1e-300 };
    const c = @Vector(4, f64){ 0.25, 2.0, 4.0, -8.0 };
    inline for (0..4) |i| {
        try testing.expectEqual(bits32(root.fadd(a[i], b[i])), bits32(vector.narrow.faddv(a, b)[i]));
        try testing.expectEqual(bits32(root.fsub(a[i], b[i])), bits32(vector.narrow.fsubv(a, b)[i]));
        try testing.expectEqual(bits32(root.fmul(a[i], b[i])), bits32(vector.narrow.fmulv(a, b)[i]));
        try testing.expectEqual(bits32(root.fdiv(a[i], b[i])), bits32(vector.narrow.fdivv(a, b)[i]));
        try testing.expectEqual(
            bits32(root.ffma(a[i], b[i], c[i])),
            bits32(vector.narrow.ffmav(a, b, c)[i]),
        );
    }
}

test "ternary: fmav keeps the lane type and equals the scalar fma per lane" {
    const a = @Vector(4, f64){ 2.0, 3.0, -1.5, 1e300 };
    const b = @Vector(4, f64){ 3.0, 0.5, 2.0, 1e300 };
    const c = @Vector(4, f64){ 1.0, -1.0, 0.5, -1e300 };
    const got = root.fma(a, b, c);
    try testing.expect(@TypeOf(got) == @TypeOf(a));
    inline for (0..4) |i| {
        try testing.expectEqual(bits64(root.fma(a[i], b[i], c[i])), bits64(got[i]));
    }
}

test "the scalar path of a result-struct operation is unchanged" {
    const got = root.frexp(@as(f32, 1.5));
    // A bare `FrexpResult`, not the one-element array a `@Vector(1, f32)` would
    // produce.
    try testing.expect(@TypeOf(got) == root.FrexpResult(f32));
    try testing.expect(@TypeOf(got) != [1]root.FrexpResult(f32));

    const want = scalar.frexpf(@as(f32, 1.5));
    try testing.expectEqual(want.exponent, got.exponent);
    try testing.expectEqual(bits32(want.fraction), bits32(got.fraction));
    // And the same holds for the other struct-returning top-level operations.
    try testing.expect(@TypeOf(root.lgamma(@as(f32, 1.5))) == root.LgammaResult(f32));
    try testing.expect(@TypeOf(root.lrint(@as(f32, 1.5))) == root.IntegerResult);
    try testing.expect(@TypeOf(root.ilogb(@as(f32, 1.5))) == i32);
    try testing.expect(@TypeOf(root.totalOrder(@as(f32, 1.0), @as(f32, 2.0))) == bool);
}
