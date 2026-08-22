// Correctly rounded inverse hyperbolic sine for binary64 values.
//
// Ported from CORE-MATH src/binary64/asinh/asinh.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa (Copyright (c) 2023-2025
// Alexei Sibidanov) and used under the MIT License.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// The caller must keep IEEE round-to-nearest, ties-to-even active.

const data = @import("scalar_asinhd_data.zig");

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

fn addDd(x_high: f64, x_low: f64, c_high: f64, c_low: f64, error_value: *f64) f64 {
    @setFloatMode(.strict);
    const sum = x_high + c_high;
    const delta = sum - x_high;
    error_value.* = ((c_high - delta) + (x_high + (delta - sum))) + (x_low + c_low);
    return sum;
}

fn mulDdAcc(x_high: f64, x_low: f64, c_high: f64, c_low: f64, error_value: *f64) f64 {
    @setFloatMode(.strict);
    const high_low = c_high * x_low;
    const low_high = c_low * x_high;
    const high_high = c_high * x_high;
    var high_low_error = @mulAdd(f64, c_high, x_high, -high_high);
    high_low_error += low_high + high_low;
    const result_high = high_high + high_low_error;
    error_value.* = (high_high - result_high) + high_low_error;
    return result_high;
}

fn mulDdD(x_high: f64, x_low: f64, c: f64, error_value: *f64) f64 {
    @setFloatMode(.strict);
    const high = c * x_high;
    const low = @mulAdd(f64, c, x_high, -high) + c * x_low;
    const result_high = high + low;
    error_value.* = (high - result_high) + low;
    return result_high;
}

fn polyDd(comptime length: usize, x_high: f64, x_low: f64, coefficients: *const [length][2]f64, error_value: *f64) f64 {
    @setFloatMode(.strict);
    var index: usize = length - 1;
    var coefficient_high = coefficients[index][0] + error_value.*;
    var coefficient_low = ((coefficients[index][0] - coefficient_high) + error_value.*) + coefficients[index][1];
    while (index > 0) : (index -= 1) {
        coefficient_high = mulDdAcc(x_high, x_low, coefficient_high, coefficient_low, &coefficient_low);
        const sum_high = coefficient_high + coefficients[index - 1][0];
        const sum_low = (coefficients[index - 1][0] - sum_high) + coefficient_high;
        coefficient_high = sum_high;
        coefficient_low += sum_low + coefficients[index - 1][1];
    }
    error_value.* = coefficient_low;
    return coefficient_high;
}

fn asAsinhZero(x: f64, x_squared_high: f64, x_squared_low: f64) f64 {
    @setFloatMode(.strict);
    const coefficients = [12][2]f64{
        .{ -0x1.5555555555555p-3, -0x1.5555555555555p-57 }, .{ 0x1.3333333333333p-4, 0x1.99999999949dfp-59 },
        .{ -0x1.6db6db6db6db7p-5, 0x1.2492496091b0cp-60 },  .{ 0x1.f1c71c71c71c7p-6, 0x1.c71a35cfa0671p-62 },
        .{ -0x1.6e8ba2e8ba2e9p-6, 0x1.17f937248cf81p-60 },  .{ 0x1.1c4ec4ec4ec4fp-6, -0x1.74e3c1dfd4c3dp-60 },
        .{ -0x1.c999999999977p-7, -0x1.38e7a467ecc55p-61 }, .{ 0x1.7a87878786c7ep-7, 0x1.a83c7bace55ebp-61 },
        .{ -0x1.3fde50d764083p-7, -0x1.d024df7fa0542p-61 }, .{ 0x1.12ef3ceae4d12p-7, -0x1.ba9c13deb261fp-61 },
        .{ -0x1.df3bd104aa267p-8, -0x1.546da9bc5b32ap-62 }, .{ 0x1.a685fc5de7a04p-8, 0x1.40d284a1d67f9p-62 },
    };
    const low_coefficients = [5]f64{
        -0x1.7828d553ec8p-8, 0x1.51712f7bee368p-8, -0x1.2e6d98527bcc6p-8, 0x1.0095da47b392cp-8, -0x1.3b92d6368192cp-9,
    };
    var polynomial_low = x_squared_high * (low_coefficients[0] + x_squared_high *
        (low_coefficients[1] + x_squared_high *
            (low_coefficients[2] + x_squared_high *
                (low_coefficients[3] + x_squared_high * low_coefficients[4]))));
    var polynomial_high = polyDd(12, x_squared_high, x_squared_low, &coefficients, &polynomial_low);
    polynomial_high = mulDdAcc(polynomial_high, polynomial_low, x_squared_high, x_squared_low, &polynomial_low);
    polynomial_high = mulDdD(polynomial_high, polynomial_low, x, &polynomial_low);
    var first_low: f64 = undefined;
    const first_high = fastTwoSum(x, polynomial_high, &first_low);
    polynomial_high = fastTwoSum(first_low, polynomial_low, &polynomial_low);

    const first_bits: u64 = @bitCast(polynomial_high);
    if ((first_bits & fraction_mask) == 0) {
        var adjusted_bits = first_bits;
        const low_bits: u64 = @bitCast(polynomial_low);
        if (((low_bits ^ first_bits) >> 63) != 0) adjusted_bits -%= 1 else adjusted_bits +%= 1;
        polynomial_high = @bitCast(adjusted_bits);
    }
    return first_high + polynomial_high;
}

fn asAsinhDatabase(x: f64, fallback: f64) f64 {
    @setFloatMode(.strict);
    const absolute_value = if (x < 0.0) -x else x;
    var lower: i32 = 0;
    var upper: i32 = @intCast(data.db.len - 1);
    while (lower <= upper) {
        const middle = @divTrunc(lower + upper, 2);
        const middle_index: usize = @intCast(middle);
        const entry = data.db[middle_index];
        if (entry[0] < absolute_value) {
            lower = middle + 1;
        } else if (entry[0] == absolute_value) {
            const sign: f64 = if (x < 0.0) -1.0 else 1.0;
            return sign * entry[1] + sign * entry[2];
        } else {
            upper = middle - 1;
        }
    }
    return fallback;
}

fn asAsinhRefine(x: f64, z_high: f64, z_low: f64, a: f64) f64 {
    @setFloatMode(.strict);
    const t1 = data.t1;
    const t2 = data.t2;
    const t3 = data.t3;
    const t4 = data.t4;
    const ll = data.LL;
    const high_coefficients = [3][2]f64{
        .{ 0x1p-1, 0x1.24b67ee516e3bp-111 },              .{ -0x1p-2, -0x1.932ce43199a8dp-110 },
        .{ 0x1.5555555555555p-3, 0x1.55540c15cf91fp-57 },
    };
    const low_coefficients = [3]f64{ -0x1p-3, 0x1.9999999a0754fp-4, -0x1.55555555c3157p-4 };
    const z_bits: u64 = @bitCast(z_high);
    const exponent: i64 = @as(i64, @intCast(z_bits >> 52)) - 0x3ff + @intFromBool(z_low == 0.0);
    var normalized_bits = z_bits & fraction_mask;
    normalized_bits |= @as(u64, 0x3ff) << 52;
    const normalized: f64 = @bitCast(normalized_bits);
    const exponent_float: f64 = @floatFromInt(exponent);
    const v_bits: u64 = @bitCast(a - exponent_float + 0x1.00008p+0);
    const table_index: u64 = (v_bits -% (@as(u64, 0x3ff) << 52)) >> (52 - 16);
    const index_1: usize = @intCast((table_index >> 12) & 0x1f);
    const index_2: usize = @intCast((table_index >> 8) & 0xf);
    const index_3: usize = @intCast((table_index >> 4) & 0xf);
    const index_4: usize = @intCast(table_index & 0xf);
    const log2_center = 0x1.62e42fefa38p-2;
    const log2_low = 0x1.ef35793c768p-46;
    const log2_tail = -0x1.9ff0342542fc3p-91;
    const log2_error_2 = log2_tail * exponent_float;
    const log2_error_1 = log2_low * exponent_float;
    const log2_error_0 = log2_center * exponent_float;
    const logarithm = [3]f64{
        ll[0][index_1][0] + ll[1][index_2][0] + (ll[2][index_3][0] + ll[3][index_4][0]) + log2_error_0,
        ll[0][index_1][1] + ll[1][index_2][1] + (ll[2][index_3][1] + ll[3][index_4][1]),
        ll[0][index_1][2] + ll[1][index_2][2] + (ll[2][index_3][2] + ll[3][index_4][2]),
    };
    const product_12 = t1[index_1] * t2[index_2];
    const product_34 = t3[index_3] * t4[index_4];
    const product_high = product_12 * product_34;
    const product_low = @mulAdd(f64, product_12, product_34, -product_high);
    const delta_high = product_high * normalized;
    const delta_low = @mulAdd(f64, product_high, normalized, -delta_high);
    const scaled_high = product_low * normalized;
    const scaled_low = @mulAdd(f64, product_low, normalized, -scaled_high);
    var x_low: f64 = undefined;
    var x_high = fastTwoSum(delta_high - 1.0, delta_low, &x_low);
    if (z_low != 0.0) {
        const z_low_bits: u64 = @bitCast(z_low);
        const adjusted_low_bits = z_low_bits -% @as(u64, @bitCast(exponent << 52));
        x_low += product_high * @as(f64, @bitCast(adjusted_low_bits));
    }
    x_high = addDd(x_high, x_low, scaled_high, scaled_low, &x_low);
    var low = x_high * (low_coefficients[0] + x_high *
        (low_coefficients[1] + x_high * low_coefficients[2]));
    var high = polyDd(3, x_high, x_low, &high_coefficients, &low);
    high = mulDdAcc(x_high, x_low, high, low, &low);
    high = addDd(high, low, log2_error_1, log2_error_2, &low);
    high = addDd(high, low, logarithm[1], logarithm[2], &low);
    var second: f64 = undefined;
    var first = fastTwoSum(logarithm[0], high, &second);
    var third = fastTwoSum(second, low, &second);
    first = fastTwoSum(first, third, &third);
    third = fastTwoSum(third, second, &second);
    const sign: f64 = if (x < 0.0) -2.0 else 2.0;
    first *= sign;
    third *= sign;
    second *= sign;
    var correction_bits: u64 = @bitCast(third);
    if ((correction_bits & fraction_mask) == 0) {
        const low_bits: u64 = @bitCast(second);
        if (((low_bits ^ correction_bits) >> 63) != 0) correction_bits -%= 1 else correction_bits +%= 1;
        third = @bitCast(correction_bits);
    }
    const first_bits: u64 = @bitCast(first);
    const error_region = (correction_bits +% 41) & fraction_mask;
    const exponent_difference = ((first_bits >> 52) & 0x7ff) - ((correction_bits >> 52) & 0x7ff);
    const result = first + third;
    if (exponent_difference > 99 or error_region < 80) return asAsinhDatabase(x, result);
    return result;
}

/// Correctly rounded inverse hyperbolic sine for binary64 values.
pub fn asinhd(value: f64) f64 {
    @setFloatMode(.strict);
    const input_bits: u64 = @bitCast(value);
    const magnitude = input_bits & magnitude_mask;
    if (magnitude > infinity_bits) return @bitCast(input_bits | quiet_nan_bit);
    if (magnitude == 0 or magnitude == infinity_bits) return value;

    const absolute_value: f64 = @bitCast(magnitude);
    if (magnitude < 0x3fbb_0000_0000_0000) {
        if (magnitude < 0x3e57_1374_4912_3ef7) {
            const result = @mulAdd(f64, -0x1p-60, value, value);
            return result;
        }
        const square_high = value * value;
        const square_low = @mulAdd(f64, value, value, -square_high);
        const cube_high = square_high * value;
        var correction: f64 = undefined;
        if (magnitude < 0x3f93_0000_0000_0000) {
            if (magnitude < 0x3f30_0000_0000_0000) {
                if (magnitude < 0x3e5a_0000_0000_0000) {
                    correction = cube_high * -0x1.5555555555555p-3;
                } else {
                    correction = cube_high * (-0x1.5555555555555p-3 + square_high *
                        0x1.3333327c57c6p-4);
                }
            } else {
                correction = cube_high * (-0x1.5555555555555p-3 + square_high *
                    (0x1.333333332f2ffp-4 + square_high *
                        (-0x1.6db6d9a665159p-5 + square_high * 0x1.f186866d775fp-6)));
            }
        } else {
            correction = cube_high * (-0x1.5555555555555p-3 + square_high *
                (0x1.333333333331p-4 + square_high *
                    (-0x1.6db6db6da466cp-5 + square_high *
                        (0x1.f1c71c2ea7be4p-6 + square_high *
                            (-0x1.6e8b651b09d72p-6 + square_high *
                                (0x1.1c309fc0e69c2p-6 + square_high * -0x1.bab7833c1ep-7))))));
        }
        const epsilon = 0x1.79p-53 * cube_high;
        const lower = value + (correction - epsilon);
        const upper = value + (correction + epsilon);
        if (lower == upper) return lower;
        return asAsinhZero(value, square_high, square_low);
    }

    var square_high: f64 = 0.0;
    var square_low: f64 = 0.0;
    var absolute_high: f64 = undefined;
    var absolute_low: f64 = undefined;
    var offset: i64 = 0x3ff;
    if (magnitude < 0x4190_0000_0000_0000) {
        square_high = value * value;
        square_low = @mulAdd(f64, value, value, -square_high);
        var sum_high: f64 = undefined;
        var sum_low: f64 = undefined;
        if (magnitude < 0x3ff0_0000_0000_0000) {
            sum_high = fastTwoSum(1.0, square_high, &sum_low);
        } else {
            sum_high = fastTwoSum(square_high, 1.0, &sum_low);
        }
        sum_low += square_low;
        absolute_high = @sqrt(sum_high);
        const reciprocal = 0.5 / sum_high;
        absolute_low = (sum_low - @mulAdd(f64, absolute_high, absolute_high, -sum_high)) * (reciprocal * absolute_high);
        absolute_high = fastTwoSum(absolute_high, absolute_value, &sum_low);
        absolute_low += sum_low;
    } else if (magnitude < 0x4330_0000_0000_0000) {
        absolute_high = 2.0 * absolute_value;
        absolute_low = 0.5 / absolute_value;
    } else {
        offset = 0x3fe;
        absolute_high = absolute_value;
        absolute_low = 0.0;
    }

    const absolute_high_bits: u64 = @bitCast(absolute_high);
    const exponent: i64 = @as(i64, @intCast(absolute_high_bits >> 52)) - offset;
    var normalized_bits = absolute_high_bits & fraction_mask;

    const exponent_index: usize = @intCast(normalized_bits >> (52 - 5));
    const residual = normalized_bits & 0x0000_7fff_ffff_ffff;
    const index_numerator: i128 = @as(i128, normalized_bits) +
        (@as(i128, data.B[exponent_index].c0) << 33) +
        @as(i128, data.B[exponent_index].c1) * @as(i128, @intCast(residual >> 16));
    const table_index: u64 = @intCast(index_numerator >> (52 - 10));
    normalized_bits |= @as(u64, 0x3ff) << 52;
    const normalized_value: f64 = @bitCast(normalized_bits);
    const index_1: usize = @intCast(table_index >> 5);
    const index_2: usize = @intCast(table_index & 0x1f);
    const reciprocal = data.r1[index_1] * data.r2[index_2];
    const delta = @mulAdd(f64, reciprocal, normalized_value, -1.0);
    const delta_squared = delta * delta;
    const polynomial = delta_squared * ((data.c[0] + delta * data.c[1]) +
        delta_squared * ((data.c[2] + delta * data.c[3]) + delta_squared * data.c[4]));
    const log2_high = 0x1.62e42fefa38p-1;
    const log2_low = 0x1.ef35793c7673p-45;
    var log_high = log2_high * @as(f64, @floatFromInt(exponent)) + (data.l1[index_1][1] + data.l2[index_2][1]);
    var log_low = log2_low * @as(f64, @floatFromInt(exponent)) + data.l1[index_1][0] + data.l2[index_2][0] + absolute_low / absolute_high + polynomial;
    log_low += delta;
    const sign: f64 = if (value < 0.0) -1.0 else 1.0;
    log_high *= sign;
    log_low *= sign;
    const epsilon = 1.63e-19;
    const lower = log_high + (log_low - epsilon);
    const upper = log_high + (log_low + epsilon);
    if (lower == upper) return lower;
    if (absolute_value < 0x1p-2) return asAsinhZero(value, square_high, square_low);
    return asAsinhRefine(value, absolute_high, absolute_low, 0x1.71547652b82fep+0 *
        (if (lower < 0.0) -lower else lower));
}
