// Correctly rounded base-10 exponential function biased by one for binary64 values.
//
// Copyright (c) 2022-2025 Paul Zimmermann, Tom Hubrecht and Claude-Pierre Jeannerod
//
// This file is part of the CORE-MATH project
// (https://core-math.gitlabpages.inria.fr/).
//
// Ported from src/binary64/exp10m1/exp10m1.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
//
// Correct rounding requires the caller to use round-to-nearest, ties-to-even.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

const data = @import("scalar_exp10m1d_data.zig");

const magnitude_mask: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit: u64 = 0x0008_0000_0000_0000;
const ln10_high: f64 = 0x1.26bb1bbb55516p+1;
const ln10_low: f64 = -0x1.f48ad494ea3e9p-53;

const DoubleDouble = struct { high: f64, low: f64 };
const FastResult = struct { value: DoubleDouble, error_bound: f64 };

fn roundEvenFinite(value: f64) f64 {
    @setFloatMode(.strict);
    const truncated: i64 = @intFromFloat(value);
    const tv: f64 = @floatFromInt(truncated);
    const fraction = value - tv;
    const odd = (truncated & 1) != 0;
    if (fraction > 0.5 or (fraction == 0.5 and odd)) return @floatFromInt(truncated + 1);
    if (fraction < -0.5 or (fraction == -0.5 and odd)) return @floatFromInt(truncated - 1);
    return tv;
}

fn aMul(a: f64, b: f64) DoubleDouble {
    @setFloatMode(.strict);
    const high = a * b;
    return .{ .high = high, .low = @mulAdd(f64, a, b, -high) };
}
fn sMul(a: f64, bh: f64, bl: f64) DoubleDouble {
    @setFloatMode(.strict);
    var r = aMul(a, bh);
    r.low = @mulAdd(f64, a, bl, r.low);
    return r;
}
fn dMul(ah: f64, al: f64, bh: f64, bl: f64) DoubleDouble {
    @setFloatMode(.strict);
    var r = aMul(ah, bh);
    r.low = @mulAdd(f64, ah, bl, r.low);
    r.low = @mulAdd(f64, al, bh, r.low);
    return r;
}
fn fastTwoSum(a: f64, b: f64) DoubleDouble {
    @setFloatMode(.strict);
    const high = a + b;
    const e = high - a;
    return .{ .high = high, .low = b - e };
}
fn fastSum(a: f64, bh: f64, bl: f64) DoubleDouble {
    @setFloatMode(.strict);
    var r = fastTwoSum(a, bh);
    r.low += bl;
    return r;
}
fn mulAddDd(x: f64, value: DoubleDouble, ch: f64, cl: f64) DoubleDouble {
    @setFloatMode(.strict);
    var r = sMul(x, value.high, value.low);
    const sum = fastTwoSum(ch, r.high);
    r.high = sum.high;
    r.low += sum.low + cl;
    return r;
}
fn copySign(magnitude: f64, source: f64) f64 {
    const mb: u64 = @bitCast(magnitude);
    const sb = @as(u64, @bitCast(source)) & ~magnitude_mask;
    return @bitCast((mb & magnitude_mask) | sb);
}

fn q1(zh: f64, zl: f64) DoubleDouble {
    @setFloatMode(.strict);
    const z = zh + zl;
    var q = @mulAdd(f64, data.q1[4], zh, data.q1[3]);
    q = @mulAdd(f64, q, z, data.q1[2]);
    const sum = fastTwoSum(data.q1[1], q * z);
    const product = dMul(zh, zl, sum.high, sum.low);
    return fastSum(data.q1[0], product.high, product.low);
}
fn exp1(xh: f64, xl: f64) DoubleDouble {
    @setFloatMode(.strict);
    const k = roundEvenFinite(xh * 0x1.71547652b82fep+12);
    const kl = sMul(k, 0x1.62e42fefa39efp-13, 0x1.abc9e3b39803fp-68);
    var y = fastTwoSum(xh - kl.high, xl);
    y.low -= kl.low;
    const ki: i64 = @intFromFloat(k);
    var m = (ki >> 12) + 0x3ff;
    const table_index_2: usize = @intCast((ki >> 6) & 0x3f);
    const table_index_1: usize = @intCast(ki & 0x3f);
    var r = dMul(data.t2[table_index_1][0], data.t2[table_index_1][1], data.t1[table_index_2][0], data.t1[table_index_2][1]);
    const poly = q1(y.high, y.low);
    r = dMul(r.high, r.low, poly.high, poly.low);
    if (m == 0x7ff) {
        r.high *= 2.0;
        r.low *= 2.0;
        m -= 1;
    }
    const scale: f64 = @bitCast(@as(u64, @intCast(m)) << 52);
    r.high *= scale;
    r.low *= scale;
    return r;
}

fn q2(zh: f64, zl: f64) DoubleDouble {
    @setFloatMode(.strict);
    const z = zh + zl;
    var q = @mulAdd(f64, data.q2[8], zh, data.q2[7]);
    q = @mulAdd(f64, q, z, data.q2[6]);
    q = @mulAdd(f64, q, z, data.q2[5]);
    var r = aMul(q, z);
    var sum = fastTwoSum(data.q2[3], r.high);
    r.high = sum.high;
    r.low += sum.low + data.q2[4];
    r = dMul(r.high, r.low, zh, zl);
    sum = fastTwoSum(data.q2[2], r.high);
    r.high = sum.high;
    r.low += sum.low;
    r = dMul(r.high, r.low, zh, zl);
    sum = fastTwoSum(data.q2[1], r.high);
    r.high = sum.high;
    r.low += sum.low;
    r = dMul(r.high, r.low, zh, zl);
    sum = fastTwoSum(data.q2[0], r.high);
    r.high = sum.high;
    r.low += sum.low;
    return r;
}
fn exp2(value: f64) DoubleDouble {
    @setFloatMode(.strict);
    var k = roundEvenFinite(value * 0x1.a934f0979a371p+13);
    if (k == 4194304.0) k = 4194303.0;
    var yh = @mulAdd(f64, -k, 0x1.34413509f79ffp-14, value);
    var y = fastTwoSum(yh, -k * -0x1.9dc1da98p-71);
    yh = y.high;
    y.low = @mulAdd(f64, -k, -0x1.4fd20dba1f655p-103, y.low);
    const e = sMul(yh, ln10_high, ln10_low);
    y.high = e.high;
    y.low = @mulAdd(f64, y.low, ln10_high, e.low);
    const ki: i64 = @intFromFloat(k);
    var m = (ki >> 12) + 0x3ff;
    const table_index_2: usize = @intCast((ki >> 6) & 0x3f);
    const table_index_1: usize = @intCast(ki & 0x3f);
    var r = dMul(data.t2[table_index_1][0], data.t2[table_index_1][1], data.t1[table_index_2][0], data.t1[table_index_2][1]);
    const poly = q2(y.high, y.low);
    r = dMul(r.high, r.low, poly.high, poly.low);
    if (m == 0x7ff) {
        r.high *= 2.0;
        r.low *= 2.0;
        m -= 1;
    }
    const scale: f64 = @bitCast(@as(u64, @intCast(m)) << 52);
    r.high *= scale;
    r.low *= scale;
    return r;
}

fn exp10m1FastTiny(x: f64) FastResult {
    @setFloatMode(.strict);
    const x2 = x * x;
    const x4 = x2 * x2;
    var c9 = @mulAdd(f64, data.p_fast[12], x, data.p_fast[11]);
    const c7 = @mulAdd(f64, data.p_fast[10], x, data.p_fast[9]);
    var c5 = @mulAdd(f64, data.p_fast[8], x, data.p_fast[7]);
    c9 = @mulAdd(f64, data.p_fast[13], x2, c9);
    c5 = @mulAdd(f64, c7, x2, c5);
    c5 = @mulAdd(f64, c9, x4, c5);
    var r = aMul(c5, x);
    const sum = fastTwoSum(data.p_fast[6], r.high);
    r.high = sum.high;
    r.low += sum.low;
    r = mulAddDd(x, r, data.p_fast[4], data.p_fast[5]);
    r = mulAddDd(x, r, data.p_fast[2], data.p_fast[3]);
    r = mulAddDd(x, r, data.p_fast[0], data.p_fast[1]);
    r = sMul(x, r.high, r.low);
    return .{ .value = r, .error_bound = 0x1.0ap-68 * r.high };
}
fn exp10m1Fast(x: f64, tiny: bool) FastResult {
    @setFloatMode(.strict);
    if (tiny) return exp10m1FastTiny(x);
    const logarithm = sMul(x, ln10_high, ln10_low);
    var r = exp1(logarithm.high, logarithm.low);
    const sum = if (x >= 0.0) fastTwoSum(r.high, -1.0) else fastTwoSum(-1.0, r.high);
    r.high = sum.high;
    r.low += sum.low;
    return .{ .value = r, .error_bound = 0x1.7ap-72 * r.high };
}
fn exceptionLookup(comptime n: usize, x: f64, table: *const [n][3]f64) ?f64 {
    @setFloatMode(.strict);
    var a: usize = 0;
    var b: usize = n;
    while (a + 1 != b) {
        const c = (a + b) / 2;
        if (table[c][0] <= x) a = c else b = c;
    }
    if (x == table[a][0]) return table[a][1] + table[a][2];
    return null;
}
fn exp10m1AccurateTiny(x: f64) f64 {
    @setFloatMode(.strict);
    if (exceptionLookup(data.exceptions_accurate_tiny.len, x, &data.exceptions_accurate_tiny)) |r| return r;
    const x2 = x * x;
    const x4 = x2 * x2;
    var c15 = @mulAdd(f64, data.q_accurate[23], x, data.q_accurate[22]);
    const c13 = @mulAdd(f64, data.q_accurate[21], x, data.q_accurate[20]);
    c15 = @mulAdd(f64, data.q_accurate[24], x2, c15);
    var r = fastTwoSum(data.q_accurate[18], data.q_accurate[19] * x + c13 * x2 + c15 * x4);
    r = mulAddDd(x, r, data.q_accurate[17], 0.0);
    r = mulAddDd(x, r, data.q_accurate[16], 0.0);
    r = mulAddDd(x, r, data.q_accurate[14], data.q_accurate[15]);
    r = mulAddDd(x, r, data.q_accurate[12], data.q_accurate[13]);
    r = mulAddDd(x, r, data.q_accurate[10], data.q_accurate[11]);
    r = mulAddDd(x, r, data.q_accurate[8], data.q_accurate[9]);
    r = mulAddDd(x, r, data.q_accurate[6], data.q_accurate[7]);
    r = mulAddDd(x, r, data.q_accurate[4], data.q_accurate[5]);
    r = mulAddDd(x, r, data.q_accurate[2], data.q_accurate[3]);
    r = mulAddDd(x, r, data.q_accurate[0], data.q_accurate[1]);
    r = sMul(x, r.high, r.low);
    return r.high + r.low;
}
fn exp10m1Accurate(x: f64) f64 {
    @setFloatMode(.strict);
    const ax = @as(u64, @bitCast(x)) & magnitude_mask;
    if (ax <= 0x3fb0_0000_0000_0000) return exp10m1AccurateTiny(x);
    if (exceptionLookup(data.exceptions_accurate.len, x, &data.exceptions_accurate)) |r| return r;
    var v = exp2(x);
    const sum = if (x >= 0.0) fastTwoSum(v.high, -1.0) else fastTwoSum(-1.0, v.high);
    v.high = sum.high;
    v.low += sum.low;
    return v.high + v.low;
}

pub fn exp10m1d(value: f64) f64 {
    @setFloatMode(.strict);
    const ux: u64 = @bitCast(value);
    const ax = ux & magnitude_mask;
    if (ux >= 0xc030_4170_4c06_8ef0) {
        if ((ux >> 52) == 0xfff) {
            if (ux > 0xfff0_0000_0000_0000) return @bitCast(ux | quiet_nan_bit);
            return -1.0;
        }
        return -1.0 + 0x1p-54;
    } else if (ax > 0x4073_4413_509f_79fe) {
        if (ax > infinity_bits) return @bitCast(ux | quiet_nan_bit);
        if (ax == infinity_bits) return value;
        return 0x1.fffffffffffffp+1023 * value;
    } else if (ax <= 0x3c90_0000_0000_0000) {
        if (ax <= 0x3970_0000_0000_0000) {
            if (value == 0.0) return value;
            const absolute: f64 = @bitCast(ax);
            if (absolute == 0x0.086c73059343cp-1022)
                return @mulAdd(f64, -copySign(0x1.001p-537, value), 0x1p-538, copySign(0x1.36568740cb56p-1026, value));
            if (absolute == 0x0.13a7b70d0248cp-1022)
                return @mulAdd(f64, copySign(0x0.fffp-537, value), 0x1p-538, copySign(0x1.6a0f9dcb97e38p-1025, value));
            const sx = value * 0x1p106;
            var v = aMul(ln10_high, sx);
            v.low = @mulAdd(f64, ln10_low, sx, v.low);
            const h2 = (v.high + v.low) * 0x1p-106;
            v.high = @mulAdd(f64, -h2, 0x1p106, v.high);
            v.high += v.low;
            return @mulAdd(f64, v.high, 0x1p-106, h2);
        }
        if (exceptionLookup(data.exceptions_tiny.len, value, &data.exceptions_tiny)) |r| return r;
        var v = aMul(ln10_high, value);
        v.low = @mulAdd(f64, ln10_low, value, v.low);
        v.low = @mulAdd(f64, 0x1.53524c73cea69p+1 * value, value, v.low);
        return v.high + v.low;
    }
    if (ux << 15 == 0) {
        const integer: i32 = @intFromFloat(value);
        if (value == @as(f64, @floatFromInt(integer)) and integer >= 1 and integer <= 15)
            return data.exact_powers[@intCast(integer)];
    }
    const fast = exp10m1Fast(value, ax <= 0x3fb0_0000_0000_0000);
    const left = fast.value.high + (fast.value.low - fast.error_bound);
    const right = fast.value.high + (fast.value.low + fast.error_bound);
    if (left == right) return left;
    return exp10m1Accurate(value);
}
