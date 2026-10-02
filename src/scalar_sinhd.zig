const fma_impl = @import("scalar_fma.zig");

const data = @import("scalar_sinhd_data.zig");

const sign_mask: u64 = 0x8000_0000_0000_0000;
const magnitude_mask: u64 = 0x7fff_ffff_ffff_ffff;
const quiet_nan_bit: u64 = 0x0008_0000_0000_0000;
const infinity_bits: u64 = 0x7ff0_0000_0000_0000;
const fraction_mask: u64 = 0x000f_ffff_ffff_ffff;

fn fastTwoSum(x: f64, y: f64, error_value: *f64) f64 {
    @setFloatMode(.strict);
    const sum = x + y;
    const z = sum - x;
    error_value.* = y - z;
    return sum;
}

fn mulDd(x_high: f64, x_low: f64, c_high: f64, c_low: f64, error_value: *f64) f64 {
    @setFloatMode(.strict);
    const high = x_high * c_high;
    error_value.* = fma_impl.fma64(x_high, c_high, -high) +
        x_high * c_low + x_low * c_high;
    return high;
}

fn mulDdd(x_high: f64, x_low: f64, c: f64, error_value: *f64) f64 {
    @setFloatMode(.strict);
    const high = x_high * c;
    error_value.* = fma_impl.fma64(x_high, c, -high) + x_low * c;
    return high;
}

fn polyDd3(
    x_high: f64,
    x_low: f64,
    coefficients: *const [3][2]f64,
    error_value: *f64,
) f64 {
    @setFloatMode(.strict);
    var coefficient_low: f64 = undefined;
    var high = fastTwoSum(coefficients[2][0], error_value.*, &coefficient_low);
    coefficient_low += coefficients[2][1];
    high = mulDd(x_high, x_low, high, coefficient_low, &coefficient_low);
    var sum_error: f64 = undefined;
    high = fastTwoSum(coefficients[1][0], high, &sum_error);
    coefficient_low = (coefficient_low + coefficients[1][1]) + sum_error;
    high = mulDd(x_high, x_low, high, coefficient_low, &coefficient_low);
    high = fastTwoSum(coefficients[0][0], high, &sum_error);
    coefficient_low = (coefficient_low + coefficients[0][1]) + sum_error;
    error_value.* = coefficient_low;
    return high;
}
fn polyDd5(
    x_high: f64,
    x_low: f64,
    coefficients: *const [5][2]f64,
    error_value: *f64,
) f64 {
    @setFloatMode(.strict);
    var coefficient_low: f64 = undefined;
    var sum_error: f64 = undefined;
    var high = fastTwoSum(coefficients[4][0], error_value.*, &coefficient_low);
    coefficient_low += coefficients[4][1];
    var coefficient_index: i32 = 4;
    while (coefficient_index > 0) : (coefficient_index -= 1) {
        high = mulDd(x_high, x_low, high, coefficient_low, &coefficient_low);
        const index: usize = @intCast(coefficient_index - 1);
        high = fastTwoSum(coefficients[index][0], high, &sum_error);
        coefficient_low = (coefficient_low + coefficients[index][1]) + sum_error;
    }
    error_value.* = coefficient_low;
    return high;
}

fn asExpAccurate(
    x: f64,
    index: f64,
    table_high: f64,
    table_low: f64,
    error_value: *f64,
) f64 {
    @setFloatMode(.strict);
    const coefficients = [3][2]f64{
        .{ 0x1p+0, 0x1.6c16bd194535dp-94 },
        .{ 0x1p-1, -0x1.8259d904fd34fp-93 },
        .{ 0x1.5555555555555p-3, 0x1.53e93e9f26e62p-57 },
    };
    const log2_high = 0x1.62e42ffp-13;
    const log2_low = 0x1.718432a1b0e26p-47;
    const log2_tail = 0x1.9ff0342542fc3p-102;
    const reduced = x - log2_high * index;
    const reduced_low = log2_low * index;
    const reduced_tail = log2_tail * index +
        fma_impl.fma64(log2_low, index, -reduced_low);
    const reduced_high = reduced + reduced_low;
    const reduced_error = (reduced - reduced_high) + reduced_low + reduced_tail;
    var polynomial_low = reduced_high *
        (0x1.5555555555555p-5 + reduced_high *
            (0x1.11111113e93e9p-7 + reduced_high * 0x1.6c16c169400a7p-10));
    var polynomial_high = polyDd3(reduced_high, reduced_error, &coefficients, &polynomial_low);
    polynomial_high = mulDd(
        reduced_high,
        reduced_error,
        polynomial_high,
        polynomial_low,
        &polynomial_low,
    );
    polynomial_high = mulDd(
        table_high,
        table_low,
        polynomial_high,
        polynomial_low,
        &polynomial_low,
    );
    const sum_high = table_high + polynomial_high;
    const sum_low = (table_high - sum_high) + polynomial_high;
    const second_high = sum_high + table_low;
    const second_low = ((sum_high - second_high) + table_low) + sum_low;
    const result_high = second_high + polynomial_low;
    error_value.* = ((second_high - result_high) + polynomial_low) + second_low;
    return result_high;
}

fn asSinhZero(value: f64) f64 {
    @setFloatMode(.strict);
    const coefficients = [5][2]f64{
        .{ 0x1.5555555555555p-3, 0x1.555555555552fp-57 },
        .{ 0x1.1111111111111p-7, 0x1.11111115cf00dp-63 },
        .{ 0x1.a01a01a01a01ap-13, 0x1.a0011c925b85cp-73 },
        .{ 0x1.71de3a556c734p-19, -0x1.b4e2835532bcdp-73 },
        .{ 0x1.ae64567f54482p-26, -0x1.defcf17a6ab79p-81 },
    };
    const square_high = value * value;
    const square_low = fma_impl.fma64(value, value, -square_high);
    var polynomial_low = square_high *
        (0x1.6124613aef206p-33 + square_high *
            (0x1.ae7f36beea815p-41 + square_high * 0x1.95785063cd974p-49));
    var polynomial_high = polyDd5(square_high, square_low, &coefficients, &polynomial_low);
    polynomial_high = mulDdd(polynomial_high, polynomial_low, value, &polynomial_low);
    polynomial_high = mulDd(
        polynomial_high,
        polynomial_low,
        square_high,
        square_low,
        &polynomial_low,
    );
    var first_low: f64 = undefined;
    const first_high = fastTwoSum(value, polynomial_high, &first_low);
    polynomial_high = fastTwoSum(first_low, polynomial_low, &polynomial_low);
    const first_bits: u64 = @bitCast(polynomial_high);
    if ((first_bits & fraction_mask) == 0) {
        var adjusted_bits: u64 = first_bits;
        const low_bits: u64 = @bitCast(polynomial_low);
        if (((low_bits ^ first_bits) >> 63) != 0) {
            adjusted_bits -%= 1;
        } else {
            adjusted_bits +%= 1;
        }
        polynomial_high = @bitCast(adjusted_bits);
    }
    return first_high + polynomial_high;
}

fn asSinhDatabase(value: f64, fallback: f64) f64 {
    @setFloatMode(.strict);
    const absolute_bits = @as(u64, @bitCast(value)) & magnitude_mask;
    const absolute_value: f64 = @bitCast(absolute_bits);
    var lower: isize = 0;
    var upper: isize = @intCast(data.database.len - 1);
    while (lower <= upper) {
        const middle = @divTrunc(lower + upper, 2);
        const middle_index: usize = @intCast(middle);
        const entry = data.database[middle_index];
        if (entry[0] < absolute_value) {
            lower = middle + 1;
        } else if (entry[0] == absolute_value) {
            const sign: f64 = if ((@as(u64, @bitCast(value)) & sign_mask) == 0)
                1.0
            else
                -1.0;
            return sign * entry[1] + sign * entry[2];
        } else {
            upper = middle - 1;
        }
    }
    return fallback;
}

/// Correctly-rounded hyperbolic sine function for binary64 values.
///
/// Copyright (c) 2023-2026 Alexei Sibidanov and Paul Zimmermann.
///
/// This file is part of the CORE-MATH project
/// (https://core-math.gitlabpages.inria.fr/).
///
/// Ported from src/binary64/sinh/sinh.c at revision
/// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in all
/// copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
/// SOFTWARE.
///
/// The caller must keep IEEE round-to-nearest, ties-to-even active.
pub fn sinhd(value: f64) f64 {
    @setFloatMode(.strict);
    const input_bits: u64 = @bitCast(value);
    const absolute_bits = input_bits & magnitude_mask;

    if (absolute_bits < 0x3fd0_0000_0000_0000) {
        if (absolute_bits < 0x3e57_1374_4912_3ef7) {
            // The exact result rounds to the input encoding in this region.
            return value;
        }
        const coefficients = [5]f64{
            0x1.5555555555555p-3,
            0x1.111111111151ep-7,
            0x1.a01a019d0c767p-13,
            0x1.71de444a96e11p-19,
            0x1.ae8465375242p-26,
        };
        const square = value * value;
        const cube = square * value;
        const fourth = square * square;
        const polynomial = cube *
            ((coefficients[0] + square * coefficients[1]) +
                fourth * ((coefficients[2] + square * coefficients[3]) +
                    fourth * coefficients[4]));
        const error_bound = cube * 0x1.cp-53;
        const lower = value + (polynomial - error_bound);
        const upper = value + (polynomial + error_bound);
        if (lower == upper) return lower;
        return asSinhZero(value);
    }
    if (absolute_bits > 0x4086_33ce_8fb9_f87d) {
        if (absolute_bits >= infinity_bits) {
            if (absolute_bits == infinity_bits) return value;
            return @bitCast(input_bits | quiet_nan_bit);
        }
        return @bitCast((input_bits & sign_mask) | infinity_bits);
    }

    const scale = 0x1.71547652b82fep+12;
    const absolute_value: f64 = @bitCast(absolute_bits);
    const rounded_product = fma_impl.fma64(absolute_value, scale, 0x1.8000002p+26);
    const rounded_bits: u64 = @bitCast(rounded_product);
    const reduction_mask: u64 = ~@as(u64, (1 << 26) - 1);
    const reduction = @as(f64, @bitCast(rounded_bits & reduction_mask)) - 0x1.8p+26;
    const integer_index: i64 = @intCast((rounded_bits << 14) >> 40);
    const positive_table_1_index: usize = @intCast(integer_index & 0x3f);
    const positive_table_0_index: usize = @intCast((integer_index >> 6) & 0x3f);
    const positive_exponent: i64 = integer_index >> 12;
    const positive_scale: f64 = @bitCast(
        @as(u64, @bitCast(@as(i64, 1022 + positive_exponent))) << 52,
    );
    const table_0_high = data.t0[positive_table_0_index][1];
    const table_0_low = data.t0[positive_table_0_index][0];
    const table_1_high = data.t1[positive_table_1_index][1];
    const table_1_low = data.t1[positive_table_1_index][0];
    const table_high = table_0_high * table_1_high;
    const table_low = table_0_high * table_1_low + table_1_high * table_0_low +
        fma_impl.fma64(table_0_high, table_1_high, -table_high);
    const log2_high = 0x1.62e42ffp-13;
    const log2_low = 0x1.718432a1b0e26p-47;
    const reduced = (absolute_value - log2_high * reduction) + log2_low * reduction;
    const reduced_squared = reduced * reduced;
    const coefficients = [4]f64{ 0x1p+0, 0x1p-1, 0x1.5555555aaaaaep-3, 0x1.55555551c98cp-5 };
    const positive_high = table_high * positive_scale;
    const positive_low = table_low * positive_scale;
    const positive_polynomial = reduced *
        ((coefficients[0] + reduced * coefficients[1]) +
            reduced_squared * (coefficients[2] + reduced * coefficients[3]));

    if (absolute_bits > 0x4014_0000_0000_0000) {
        if (absolute_bits > 0x4042_5e4f_7b27_37fa) {
            const overflow_scale: f64 = @bitCast(
                @as(u64, @bitCast(@as(i64, 1021 + positive_exponent))) << 52,
            );
            var high = table_high;
            var low = table_low + table_high * positive_polynomial;
            const sign: f64 = if ((input_bits & sign_mask) == 0) 1.0 else -1.0;
            high *= sign;
            low *= sign;
            const error_bound = 0x1.1b6p-63 * table_high;
            const lower = high + (low - error_bound);
            const upper = high + (low + error_bound);
            if (lower == upper) return (lower * overflow_scale) * 2.0;
            var accurate_low: f64 = undefined;
            var accurate_high = asExpAccurate(absolute_value, reduction, table_high, table_low, &accurate_low);
            accurate_high = fastTwoSum(accurate_high, accurate_low, &accurate_low);
            accurate_high *= sign;
            accurate_low *= sign;
            const high_bits: u64 = @bitCast(accurate_high);
            const low_bits: u64 = @bitCast(accurate_low);
            const high_exponent: i64 = @intCast((high_bits >> 52) & 0x7ff);
            const low_exponent: i64 = @intCast((low_bits >> 52) & 0x7ff);
            const low_magnitude = (low_bits + 8) & fraction_mask;
            accurate_high += accurate_low;
            accurate_high *= 2.0;
            accurate_high *= overflow_scale;
            if (low_magnitude <= 16 or high_exponent - low_exponent > 103) {
                return asSinhDatabase(value, accurate_high);
            }
            return accurate_high;
        }

        const negative_integer_index = -integer_index;
        const negative_table_1_index: usize = @intCast(negative_integer_index & 0x3f);
        const negative_table_0_index: usize = @intCast((negative_integer_index >> 6) & 0x3f);
        const negative_exponent: i64 = negative_integer_index >> 12;
        const negative_scale: f64 = @bitCast(
            @as(u64, @bitCast(@as(i64, 1022 + negative_exponent))) << 52,
        );
        const negative_table_0_high = data.t0[negative_table_0_index][1];
        const negative_table_0_low = data.t0[negative_table_0_index][0];
        const negative_table_1_high = data.t1[negative_table_1_index][1];
        const negative_table_1_low = data.t1[negative_table_1_index][0];
        const negative_table_high = negative_table_0_high * negative_table_1_high;
        const negative_table_product_high = negative_table_high * negative_scale;
        const negative_table_low_unscaled = negative_table_0_high * negative_table_1_low +
            negative_table_1_high * negative_table_0_low +
            fma_impl.fma64(negative_table_0_high, negative_table_1_high, -negative_table_high);
        const negative_table_product_low = negative_table_low_unscaled * negative_scale;
        const negative_reduced = -reduced;
        const negative_even = 0x1p+0 + reduced_squared * 0x1.5555555aaaaaep-3;
        const negative_odd = 0x1p-1 + reduced_squared * 0x1.55555551c98cp-5;
        const negative_exponential = negative_table_product_high +
            negative_table_product_high * (negative_reduced *
                (negative_even + negative_reduced * negative_odd));
        var high = positive_high;
        var low = (positive_low - negative_exponential) + positive_high * positive_polynomial;
        const sign: f64 = if ((input_bits & sign_mask) == 0) 1.0 else -1.0;
        high *= sign;
        low *= sign;
        const error_bound = 0x1.202p-63 * high;
        const lower = high + (low - error_bound);
        const upper = high + (low + error_bound);
        if (lower == upper) return lower;
        var accurate_low: f64 = undefined;
        const accurate_high = asExpAccurate(absolute_value, reduction, positive_high, positive_low, &accurate_low);
        var result_high: f64 = undefined;
        var result_low: f64 = undefined;
        if (absolute_bits > 0x403f_6666_6666_6666) {
            result_high = accurate_high - negative_table_product_high;
            result_low = ((accurate_high - result_high) - negative_table_product_high) + accurate_low;
        } else {
            var accurate_negative_low: f64 = undefined;
            const accurate_negative_high = asExpAccurate(
                -absolute_value,
                -reduction,
                negative_table_product_high,
                negative_table_product_low,
                &accurate_negative_low,
            );
            result_high = accurate_high - accurate_negative_high;
            result_low = (((accurate_high - result_high) - accurate_negative_high) -
                accurate_negative_low) + accurate_low;
        }
        high = fastTwoSum(result_high, result_low, &result_low);
        const high_bits: u64 = @bitCast(high);
        const low_bits: u64 = @bitCast(result_low);
        const high_exponent: i64 = @intCast((high_bits >> 52) & 0x7ff);
        const low_exponent: i64 = @intCast((low_bits >> 52) & 0x7ff);
        const low_magnitude = (low_bits + 8) & fraction_mask;
        const sign_final: f64 = if ((input_bits & sign_mask) == 0) 1.0 else -1.0;
        high *= sign_final;
        result_low *= sign_final;
        high += result_low;
        if (low_magnitude <= 16 or high_exponent - low_exponent > 103) {
            return asSinhDatabase(value, high);
        }
        return high;
    }

    const negative_integer_index = -integer_index;
    const negative_table_1_index: usize = @intCast(negative_integer_index & 0x3f);
    const negative_table_0_index: usize = @intCast((negative_integer_index >> 6) & 0x3f);
    const negative_exponent: i64 = negative_integer_index >> 12;
    const negative_scale: f64 = @bitCast(
        @as(u64, @bitCast(@as(i64, 1022 + negative_exponent))) << 52,
    );
    const negative_table_0_high = data.t0[negative_table_0_index][1];
    const negative_table_0_low = data.t0[negative_table_0_index][0];
    const negative_table_1_high = data.t1[negative_table_1_index][1];
    const negative_table_1_low = data.t1[negative_table_1_index][0];
    const negative_high_unscaled = negative_table_0_high * negative_table_1_high;
    const negative_low_unscaled = negative_table_0_high * negative_table_1_low +
        negative_table_1_high * negative_table_0_low +
        fma_impl.fma64(negative_table_0_high, negative_table_1_high, -negative_high_unscaled);
    const negative_high_scaled = negative_high_unscaled * negative_scale;
    const negative_low_scaled = negative_low_unscaled * negative_scale;
    const negative_reduced = -reduced;
    const negative_polynomial = negative_reduced *
        ((coefficients[0] + negative_reduced * coefficients[1]) +
            reduced_squared * (coefficients[2] + negative_reduced * coefficients[3]));
    const positive_high_result = positive_high;
    const positive_low_result = positive_low + positive_high * positive_polynomial;
    const negative_high_result = negative_high_scaled;
    const negative_low_result = negative_low_scaled + negative_high_scaled * negative_polynomial;
    var high = positive_high_result - negative_high_result;
    var low = ((positive_high_result - high) - negative_high_result) -
        negative_low_result + positive_low_result;
    const sign: f64 = if ((input_bits & sign_mask) == 0) 1.0 else -1.0;
    high *= sign;
    low *= sign;
    const error_bound = 0x1.c0ap-62 * high;
    const lower = high + (low - error_bound);
    const upper = high + (low + error_bound);
    if (lower == upper) return lower;
    var accurate_low: f64 = undefined;
    const accurate_positive_high = asExpAccurate(
        absolute_value,
        reduction,
        positive_high,
        positive_low,
        &accurate_low,
    );
    var accurate_negative_low: f64 = undefined;
    const accurate_negative_high = asExpAccurate(
        -absolute_value,
        -reduction,
        negative_high_scaled,
        negative_low_scaled,
        &accurate_negative_low,
    );
    high = accurate_positive_high - accurate_negative_high;
    low = ((accurate_positive_high - high) - accurate_negative_high) -
        accurate_negative_low + accurate_low;
    high = fastTwoSum(high, low, &low);
    const high_bits: u64 = @bitCast(high);
    const low_bits: u64 = @bitCast(low);
    const high_exponent: i64 = @intCast((high_bits >> 52) & 0x7ff);
    const low_exponent: i64 = @intCast((low_bits >> 52) & 0x7ff);
    const low_magnitude = (low_bits + 8) & fraction_mask;
    high *= sign;
    low *= sign;
    high += low;
    if (low_magnitude <= 16 or high_exponent - low_exponent > 103) {
        return asSinhDatabase(value, high);
    }
    return high;
}
