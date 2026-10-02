// Correctly rounded 10^x function for binary64 values.
//
// Copyright (c) 2023-2025 Alexei Sibidanov.
//
// This file is part of the CORE-MATH project
// (https://core-math.gitlabpages.inria.fr/).
//
// Ported from src/binary64/exp10/exp10.c at revision
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

const fma_impl = @import("scalar_fma.zig");

const data = @import("scalar_exp10d_data.zig");

const sign_mask: u64 = 0x8000_0000_0000_0000;
const magnitude_mask: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit: u64 = 0x0008_0000_0000_0000;

const DoubleDouble = data.DoubleDouble;
const Rounded = struct { value: f64, integer: i64 };

fn roundEvenFinite(value: f64) Rounded {
    @setFloatMode(.strict);

    // All callers bound |value| below 2^24. The biased encoding gives both
    // the rounded value and its exact signed integer without a checked cast.
    const rounding_constant: f64 = 0x1.8p52;
    const biased = value + rounding_constant;
    const biased_bits: u64 = @bitCast(biased);
    return .{
        .value = biased - rounding_constant,
        .integer = @bitCast(biased_bits -% 0x4338_0000_0000_0000),
    };
}

fn fastTwoSum(first: f64, second: f64) DoubleDouble {
    @setFloatMode(.strict);

    const sum = first + second;
    const rounded_first = sum - first;
    return .{ .high = sum, .low = second - rounded_first };
}

fn fastSum(first_high: f64, first_low: f64, second_high: f64, second_low: f64) DoubleDouble {
    @setFloatMode(.strict);

    const sum = fastTwoSum(first_high, second_high);
    return .{ .high = sum.high, .low = (first_low + second_low) + sum.low };
}

fn multiplyDoubleDouble(
    first_high: f64,
    first_low: f64,
    second_high: f64,
    second_low: f64,
) DoubleDouble {
    @setFloatMode(.strict);

    const high = second_high * first_high;
    const low = second_low * first_high + second_high * first_low +
        fma_impl.fma64(second_high, first_high, -high);
    return .{ .high = high, .low = low };
}

fn polynomialDoubleDouble(value_high: f64, value_low: f64) DoubleDouble {
    @setFloatMode(.strict);

    const coefficients = [_]DoubleDouble{
        .{ .high = 0x1.26bb1bbb55516p+1, .low = -0x1.f48ad494ea102p-53 },
        .{ .high = 0x1.53524c73cea69p+1, .low = -0x1.e2bfab318d399p-53 },
        .{ .high = 0x1.0470591de2ca4p+1, .low = 0x1.81f50779e162bp-53 },
        .{ .high = 0x1.2bd7609fd98c4p+0, .low = 0x1.31a5cc5d3d313p-54 },
        .{ .high = 0x1.1429ffd336aa3p-1, .low = 0x1.910de8c68a0c2p-55 },
        .{ .high = 0x1.a7ed7086882b4p-3, .low = -0x1.05e703d496537p-57 },
    };
    var index = coefficients.len - 1;
    var high = coefficients[index].high;
    var low = coefficients[index].low;
    while (index != 0) {
        index -= 1;
        const product = multiplyDoubleDouble(value_high, value_low, high, low);
        high = product.high;
        low = product.low;
        const sum_high = high + coefficients[index].high;
        const sum_low = (coefficients[index].high - sum_high) + high;
        high = sum_high;
        low += sum_low + coefficients[index].low;
    }
    return .{ .high = high, .low = low };
}

fn scaleNormal(value: f64, exponent: i64) f64 {
    const value_bits: u64 = @bitCast(value);
    const exponent_bits = @as(u64, @bitCast(exponent)) << 52;
    return @bitCast(value_bits +% exponent_bits);
}

fn toDenormal(value: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    return @bitCast(value_bits & 0x000f_ffff_ffff_ffff);
}

fn databaseLookup(value: f64, fallback: f64) f64 {
    @setFloatMode(.strict);

    const value_bits: u64 = @bitCast(value);
    var lower: isize = 0;
    var upper: isize = @intCast(data.database.len - 1);
    while (lower <= upper) {
        const middle: isize = @divTrunc(lower + upper, 2);
        const middle_index: usize = @intCast(middle);
        const database_bits: u64 = @bitCast(data.database[middle_index]);
        if (database_bits < value_bits) {
            lower = middle + 1;
        } else if (database_bits > value_bits) {
            upper = middle - 1;
        } else {
            const rounding_patterns = [2]u64{ 0x7eb3_7ef5_ac3f_e7c6, 0x0000_0003_781b_19e1 };
            const sign_pattern: u64 = 371470981966157;
            const middle_shift: u6 = @intCast(middle);
            const pattern = (rounding_patterns[middle_index >> 5] >> @intCast((middle * 2) & 63)) & 3;
            const delta_bits = @as(u64, 0x3c90 | ((sign_pattern >> middle_shift) << 15)) << 48;
            const delta: f64 = @bitCast(delta_bits);
            const fallback_bits: u64 = @bitCast(fallback);
            var offset: i32 = -1;
            while (offset <= 1) : (offset += 1) {
                const candidate_bits = if (offset < 0)
                    fallback_bits -% 1
                else if (offset > 0)
                    fallback_bits +% 1
                else
                    fallback_bits;
                if ((candidate_bits & 3) == pattern) return @as(f64, @bitCast(candidate_bits)) + delta;
            }
            break;
        }
    }
    return fallback;
}

noinline fn accurate(value: f64) f64 {
    @setFloatMode(.strict);

    const value_bits: u64 = @bitCast(value);
    const reduction = roundEvenFinite(0x1.a934f0979a371p+13 * value);
    const scaled = reduction.value;
    const rounded_integer = reduction.integer;
    const table_1_index: usize = @intCast(rounded_integer & 0x3f);
    const table_0_index: usize = @intCast((rounded_integer >> 6) & 0x3f);
    const exponent = rounded_integer >> 12;
    const table_0 = data.table_0[table_0_index];
    const table_1 = data.table_1[table_1_index];
    const table_product = multiplyDoubleDouble(table_0.high, table_0.low, table_1.high, table_1.low);
    const low_0: f64 = 0x1.34413508p-14;
    const low_1: f64 = -0x1.f79fef311f12bp-46;
    const low_2: f64 = -0x1.ac0b7c917826bp-101;
    const high_delta = value - low_0 * scaled;
    const low_delta = low_1 * scaled;
    const tail_delta = low_2 * scaled + fma_impl.fma64(low_1, scaled, -low_delta);
    const reduced_sum = high_delta + low_delta;
    const reduced_high = reduced_sum;
    const reduced_low = ((high_delta - reduced_sum) + low_delta) + tail_delta;
    var polynomial = polynomialDoubleDouble(reduced_high, reduced_low);
    polynomial = multiplyDoubleDouble(reduced_high, reduced_low, polynomial.high, polynomial.low);

    if (value_bits < 0xc073_3a71_46f7_2a42) {
        if ((rounded_integer & 0xfff) == 0) {
            const first_sum = fastTwoSum(polynomial.high, polynomial.low);
            const second_sum = fastTwoSum(table_product.high, first_sum.high);
            const result = fastTwoSum(second_sum.low, first_sum.low);
            var result_high = result.high;
            if ((@as(u64, @bitCast(result_high)) << 12) == 0) {
                const result_low_bits: u64 = @bitCast(result.low);
                const same_sign = ((@as(u64, @bitCast(result_high)) ^ result_low_bits) & sign_mask) == 0;
                const delta: u64 = if (same_sign)
                    @as(u64, 1) << 51
                else
                    ~(@as(u64, 1) << 51);
                result_high = @bitCast(@as(u64, @bitCast(result_high)) +% delta);
            }
            const combined = second_sum.high + result_high;
            return scaleNormal(combined, exponent);
        }

        const product = multiplyDoubleDouble(
            polynomial.high,
            polynomial.low,
            table_product.high,
            table_product.low,
        );
        var result = fastSum(
            table_product.high,
            table_product.low,
            product.high,
            product.low,
        );
        result = fastTwoSum(result.high, result.low);
        const residual_bits: u64 = @bitCast(result.low);
        const residual_exponent = (residual_bits >> 52) & 0x7ff;
        if (((residual_bits +% 4) & 0x000f_ffff_ffff_ffff) <= 4 or residual_exponent < 918) {
            result.high = databaseLookup(value, result.high);
        }
        return scaleNormal(result.high, exponent);
    }

    const scale_bits = @as(u64, @intCast(1 - exponent)) << 52;
    const scale: f64 = @bitCast(scale_bits);
    const product = multiplyDoubleDouble(
        polynomial.high,
        polynomial.low,
        table_product.high,
        table_product.low,
    );
    var result = fastSum(
        table_product.high,
        table_product.low,
        product.high,
        product.low,
    );
    const scaled_sum = fastTwoSum(scale, result.high);
    result.high = scaled_sum.high;
    result.low += scaled_sum.low;
    return toDenormal(result.high + result.low);
}

/// Computes 10^value in binary64 with caller-supplied nearest-even rounding.
pub fn exp10d(value: f64) f64 {
    @setFloatMode(.strict);

    const bits: u64 = @bitCast(value);
    const absolute_bits = bits & magnitude_mask;
    if (absolute_bits > 0x4073_4413_509f_79fe) {
        if (absolute_bits > infinity_bits) return @bitCast(bits | quiet_nan_bit);
        if (absolute_bits == infinity_bits) {
            if ((bits & sign_mask) != 0) return 0.0;
            return value;
        }
        if ((bits & sign_mask) == 0) return @bitCast(infinity_bits);
        if (absolute_bits > 0x4074_39b7_46e3_6b52) {
            return @as(f64, 0x1.8p-1022) * 0x1p-55;
        }
    }

    if (absolute_bits <= 0x3c7b_cb7b_1526_e50e) return 1.0 + value;

    const reduction = roundEvenFinite(0x1.a934f0979a371p+13 * value);
    const scaled = reduction.value;
    const rounded_integer = reduction.integer;
    const table_1_index: usize = @intCast(rounded_integer & 0x3f);
    const table_0_index: usize = @intCast((rounded_integer >> 6) & 0x3f);
    const exponent = rounded_integer >> 12;
    const table_0 = data.table_0[table_0_index];
    const table_1 = data.table_1[table_1_index];
    const table_product = multiplyDoubleDouble(table_0.high, table_0.low, table_1.high, table_1.low);
    const low_0: f64 = 0x1.34413508p-14;
    const low_1: f64 = 0x1.f79fef311f12bp-46;
    const reduced = (value - low_0 * scaled) - low_1 * scaled;
    const reduced_squared = reduced * reduced;
    const coefficients = [_]f64{
        0x1.26bb1bbb55516p+1,
        0x1.53524c73cea69p+1,
        0x1.0470591fd74e1p+1,
        0x1.2bd760a1f32a5p+0,
    };
    const polynomial = (coefficients[0] + reduced * coefficients[1]) +
        reduced_squared * (coefficients[2] + reduced * coefficients[3]);
    const high = table_product.high;
    const high_times_reduced = high * reduced;
    var low = table_product.low + high_times_reduced * polynomial;
    const error_bound: f64 = 1.63e-19;

    if (bits < 0xc073_3a71_46f7_2a42) {
        const upper = high + (low + error_bound);
        const lower = high + (low - error_bound);
        if (upper != lower) return accurate(value);
        return scaleNormal(high + low, exponent);
    }

    const scale_bits = @as(u64, @intCast(1 - exponent)) << 52;
    const scale: f64 = @bitCast(scale_bits);
    const scaled_sum = fastTwoSum(scale, high);
    const scaled_high = scaled_sum.high;
    low += scaled_sum.low;
    const upper = scaled_high + (low + error_bound);
    const lower = scaled_high + (low - error_bound);
    if (upper != lower) return accurate(value);
    return toDenormal(lower);
}
