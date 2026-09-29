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
