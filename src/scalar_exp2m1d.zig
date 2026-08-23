// Correctly rounded base-2 exponential biased by one for binary64 values.
//
// Copyright (c) 2022-2025 Paul Zimmermann, Tom Hubrecht and Claude-Pierre Jeannerod
//
// This file is part of the CORE-MATH project
// (https://core-math.gitlabpages.inria.fr/).
//
// Ported from src/binary64/exp2m1/exp2m1.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
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
//
// Correct rounding requires the caller to use round-to-nearest, ties-to-even.

const data = @import("scalar_exp2m1d_data.zig");

const magnitude_mask: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit: u64 = 0x0008_0000_0000_0000;
const ln2_high: f64 = 0x1.62e42fefa39efp-1;
const ln2_low: f64 = 0x1.abc9e3b39803fp-56;

const DoubleDouble = struct { high: f64, low: f64 };
const FastResult = struct { value: DoubleDouble, error_bound: f64 };

fn roundEvenFinite(value: f64) f64 {
    @setFloatMode(.strict);
    const truncated: i64 = @intFromFloat(value);
    const truncated_float: f64 = @floatFromInt(truncated);
    const fraction = value - truncated_float;
    const odd = (truncated & 1) != 0;
    if (fraction > 0.5 or (fraction == 0.5 and odd)) return @floatFromInt(truncated + 1);
    if (fraction < -0.5 or (fraction == -0.5 and odd)) return @floatFromInt(truncated - 1);
    return truncated_float;
}

fn aMul(a: f64, b: f64) DoubleDouble {
    @setFloatMode(.strict);
    const high = a * b;
    return .{ .high = high, .low = @mulAdd(f64, a, b, -high) };
}

fn sMul(a: f64, bh: f64, bl: f64) DoubleDouble {
    @setFloatMode(.strict);
    var result = aMul(a, bh);
    result.low = @mulAdd(f64, a, bl, result.low);
    return result;
}

fn dMul(ah: f64, al: f64, bh: f64, bl: f64) DoubleDouble {
    @setFloatMode(.strict);
    var result = aMul(ah, bh);
    result.low = @mulAdd(f64, ah, bl, result.low);
    result.low = @mulAdd(f64, al, bh, result.low);
    return result;
}

fn fastTwoSum(a: f64, b: f64) DoubleDouble {
    @setFloatMode(.strict);
    const high = a + b;
    const rounded_a = high - a;
    return .{ .high = high, .low = b - rounded_a };
}

fn fastSum(a: f64, bh: f64, bl: f64) DoubleDouble {
    @setFloatMode(.strict);
    var result = fastTwoSum(a, bh);
    result.low += bl;
    return result;
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
    const logarithm = sMul(k, 0x1.62e42fefa39efp-13, 0x1.abc9e3b39803fp-68);
    var y = fastTwoSum(xh - logarithm.high, xl);
    y.low -= logarithm.low;

    const k_integer: i64 = @intFromFloat(k);
    var exponent = (k_integer >> 12) + 0x3ff;
    const table_2_index: usize = @intCast((k_integer >> 6) & 0x3f);
    const table_1_index: usize = @intCast(k_integer & 0x3f);
    var result = dMul(
        data.t2[table_1_index][0],
        data.t2[table_1_index][1],
        data.t1[table_2_index][0],
        data.t1[table_2_index][1],
    );
    const polynomial = q1(y.high, y.low);
    result = dMul(result.high, result.low, polynomial.high, polynomial.low);
    if (exponent == 0x7ff) {
        result.high *= 2.0;
        result.low *= 2.0;
        exponent -= 1;
    }
    const scale: f64 = @bitCast(@as(u64, @intCast(exponent)) << 52);
    result.high *= scale;
    result.low *= scale;
    return result;
}

fn q2(zh: f64, zl: f64) DoubleDouble {
    @setFloatMode(.strict);
    const z = zh + zl;
    var q = @mulAdd(f64, data.q2[8], zh, data.q2[7]);
    q = @mulAdd(f64, q, z, data.q2[6]);
    q = @mulAdd(f64, q, z, data.q2[5]);
    var result = aMul(q, z);
    var sum = fastTwoSum(data.q2[3], result.high);
    result.high = sum.high;
    result.low += sum.low + data.q2[4];
    result = dMul(result.high, result.low, zh, zl);
    sum = fastTwoSum(data.q2[2], result.high);
    result.high = sum.high;
    result.low += sum.low;
    result = dMul(result.high, result.low, zh, zl);
    sum = fastTwoSum(data.q2[1], result.high);
    result.high = sum.high;
    result.low += sum.low;
    result = dMul(result.high, result.low, zh, zl);
    sum = fastTwoSum(data.q2[0], result.high);
    result.high = sum.high;
    result.low += sum.low;
    return result;
}

fn exp2Accurate(value: f64) DoubleDouble {
    @setFloatMode(.strict);
    const k = roundEvenFinite(value * 0x1p12);
    var y_high = @mulAdd(f64, -k, 0x1p-12, value);
    const y = sMul(y_high, ln2_high, ln2_low);
    y_high = y.high;

    const k_integer: i64 = @intFromFloat(k);
    var exponent = (k_integer >> 12) + 0x3ff;
    const table_2_index: usize = @intCast((k_integer >> 6) & 0x3f);
    const table_1_index: usize = @intCast(k_integer & 0x3f);
    var result = dMul(
        data.t2[table_1_index][0],
        data.t2[table_1_index][1],
        data.t1[table_2_index][0],
        data.t1[table_2_index][1],
    );
    const polynomial = q2(y_high, y.low);
    result = dMul(result.high, result.low, polynomial.high, polynomial.low);
    if (exponent == 0x7ff) {
        result.high *= 2.0;
        result.low *= 2.0;
        exponent -= 1;
    }
    const scale: f64 = @bitCast(@as(u64, @intCast(exponent)) << 52);
    result.high *= scale;
    result.low *= scale;
    return result;
}

fn exp2m1FastTiny(value: f64) FastResult {
    @setFloatMode(.strict);
    const value_squared = value * value;
    const value_fourth = value_squared * value_squared;
    var coefficient_8 = @mulAdd(f64, data.p_fast[10], value, data.p_fast[9]);
    const coefficient_6 = @mulAdd(f64, data.p_fast[8], value, data.p_fast[7]);
    var coefficient_4 = @mulAdd(f64, data.p_fast[6], value, data.p_fast[5]);
    coefficient_8 = @mulAdd(f64, data.p_fast[11], value_squared, coefficient_8);
    coefficient_4 = @mulAdd(f64, coefficient_6, value_squared, coefficient_4);
    coefficient_4 = @mulAdd(f64, coefficient_8, value_fourth, coefficient_4);

    var result = aMul(coefficient_4, value);
    var sum = fastTwoSum(data.p_fast[4], result.high);
    result.high = sum.high;
    result.low += sum.low;
    result = sMul(value, result.high, result.low);
    sum = fastTwoSum(data.p_fast[2], result.high);
    result.high = sum.high;
    result.low += sum.low + data.p_fast[3];
    result = sMul(value, result.high, result.low);
    sum = fastTwoSum(data.p_fast[0], result.high);
    result.high = sum.high;
    result.low += sum.low + data.p_fast[1];
    result = sMul(value, result.high, result.low);
    return .{ .value = result, .error_bound = 0x1.4ep-66 * result.high };
}

fn exp2m1Fast(value: f64, tiny: bool) FastResult {
    @setFloatMode(.strict);
    if (tiny) return exp2m1FastTiny(value);
    const logarithm = sMul(value, ln2_high, ln2_low);
    var result = exp1(logarithm.high, logarithm.low);
    const sum = if (value >= 0.0)
        fastTwoSum(result.high, -1.0)
    else
        fastTwoSum(-1.0, result.high);
    result.high = sum.high;
    result.low += sum.low;
    return .{ .value = result, .error_bound = 0x1.42p-71 * result.high };
}

fn exceptionLookup(comptime count: usize, value: f64, table: *const [count][3]f64) ?f64 {
    @setFloatMode(.strict);
    var lower: usize = 0;
    var upper: usize = count;
    while (lower + 1 != upper) {
        const middle = (lower + upper) / 2;
        if (table[middle][0] <= value) lower = middle else upper = middle;
    }
    if (value == table[lower][0]) return table[lower][1] + table[lower][2];
    return null;
}

fn exp2m1AccurateTiny(value: f64) f64 {
    @setFloatMode(.strict);
    if (exceptionLookup(
        data.exceptions_accurate_tiny.len,
        value,
        &data.exceptions_accurate_tiny,
    )) |result| return result;

    const value_squared = value * value;
    const value_fourth = value_squared * value_squared;
    var coefficient_13 = @mulAdd(f64, data.q_accurate[20], value, data.q_accurate[19]);
    const coefficient_11 = @mulAdd(f64, data.q_accurate[18], value, data.q_accurate[17]);
    coefficient_13 = @mulAdd(f64, data.q_accurate[21], value_squared, coefficient_13);
    var result = fastTwoSum(
        data.q_accurate[15],
        data.q_accurate[16] * value +
            coefficient_11 * value_squared + coefficient_13 * value_fourth,
    );
    result = sMul(value, result.high, result.low);
    var sum = fastTwoSum(data.q_accurate[14], result.high);
    result.high = sum.high;
    result.low += sum.low;
    result = sMul(value, result.high, result.low);
    sum = fastTwoSum(data.q_accurate[12], result.high);
    result.high = sum.high;
    result.low += sum.low + data.q_accurate[13];
    result = sMul(value, result.high, result.low);
    sum = fastTwoSum(data.q_accurate[10], result.high);
    result.high = sum.high;
    result.low += sum.low + data.q_accurate[11];
    result = sMul(value, result.high, result.low);
    sum = fastTwoSum(data.q_accurate[8], result.high);
    result.high = sum.high;
    result.low += sum.low + data.q_accurate[9];
    result = sMul(value, result.high, result.low);
    sum = fastTwoSum(data.q_accurate[6], result.high);
    result.high = sum.high;
    result.low += sum.low + data.q_accurate[7];
    result = sMul(value, result.high, result.low);
    sum = fastTwoSum(data.q_accurate[4], result.high);
    result.high = sum.high;
    result.low += sum.low + data.q_accurate[5];
    result = sMul(value, result.high, result.low);
    sum = fastTwoSum(data.q_accurate[2], result.high);
    result.high = sum.high;
    result.low += sum.low + data.q_accurate[3];
    result = sMul(value, result.high, result.low);
    sum = fastTwoSum(data.q_accurate[0], result.high);
    result.high = sum.high;
    result.low += sum.low + data.q_accurate[1];
    result = sMul(value, result.high, result.low);
    return result.high + result.low;
}

fn exp2m1Accurate(value: f64) f64 {
    @setFloatMode(.strict);
    const absolute_bits = @as(u64, @bitCast(value)) & magnitude_mask;
    if (absolute_bits <= 0x3fc0_0000_0000_0000) return exp2m1AccurateTiny(value);
    if (exceptionLookup(data.exceptions_accurate.len, value, &data.exceptions_accurate)) |result|
        return result;
    var result = exp2Accurate(value);
    const sum = if (value >= 0.0)
        fastTwoSum(result.high, -1.0)
    else
        fastTwoSum(-1.0, result.high);
    result.high = sum.high;
    result.low += sum.low;
    return result.high + result.low;
}

/// Correctly-rounded base-2 exponential biased by one for binary64 values.
/// Correct rounding requires the caller to use round-to-nearest, ties-to-even.
pub fn exp2m1d(value: f64) f64 {
    @setFloatMode(.strict);
    const input_bits: u64 = @bitCast(value);
    const absolute_bits = input_bits & magnitude_mask;

    if (input_bits >= 0xc04b_0000_0000_0000) {
        if ((input_bits >> 52) == 0xfff) {
            if (input_bits > 0xfff0_0000_0000_0000)
                return @bitCast(input_bits | quiet_nan_bit);
            return -1.0;
        }
        return -1.0 + 0x1p-54;
    }
    if (absolute_bits >= 0x4090_0000_0000_0000) {
        if (absolute_bits > infinity_bits) return @bitCast(input_bits | quiet_nan_bit);
        if (absolute_bits == infinity_bits) return value;
        return 0x1.fffffffffffffp+1023 + value * 0x1.fffffffffffffp+960;
    }
    if (absolute_bits <= 0x3cc0_527d_bd87_e24d) {
        if (absolute_bits <= 0x3970_0000_0000_0000) {
            if (value == 0.0) return value;
            const scaled_value = value * 0x1p106;
            var result = aMul(ln2_high, scaled_value);
            result.low = @mulAdd(f64, ln2_low, scaled_value, result.low);
            const rounded = (result.high + result.low) * 0x1p-106;
            result.high = @mulAdd(f64, -rounded, 0x1p106, result.high);
            result.high += result.low;
            return @mulAdd(f64, result.high, 0x1p-106, rounded);
        }
        if (exceptionLookup(data.exceptions_tiny.len, value, &data.exceptions_tiny)) |result|
            return result;
        var result = aMul(ln2_high, value);
        result.low = @mulAdd(f64, ln2_low, value, result.low);
        result.low += 0x1.ebfbdff82c58fp-3 * value * value;
        return result.high + result.low;
    }

    if (input_bits << 17 == 0) {
        const integer: i64 = @intFromFloat(@floor(value));
        if (value == @as(f64, @floatFromInt(integer)) and integer >= -53 and integer <= 53) {
            if (integer >= 0) {
                const shift: u6 = @intCast(integer);
                return @floatFromInt((@as(u64, 1) << shift) - 1);
            }
            const power_bits = @as(u64, @intCast(integer + 1023)) << 52;
            return -1.0 + @as(f64, @bitCast(power_bits));
        }
    }

    const fast = exp2m1Fast(value, absolute_bits <= 0x3fc0_0000_0000_0000);
    const left = fast.value.high + (fast.value.low - fast.error_bound);
    const right = fast.value.high + (fast.value.low + fast.error_bound);
    if (left == right) return left;
    return exp2m1Accurate(value);
}
