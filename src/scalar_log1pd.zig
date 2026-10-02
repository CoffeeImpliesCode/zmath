// Correctly rounded log(1+x) for binary64 values.
//
// Copyright (c) 2024-2025 Alexei Sibidanov.
//
// This file is part of the CORE-MATH project
// (https://core-math.gitlabpages.inria.fr/).
//
// This code is derived from CORE-MATH src/binary64/log1p/log1p.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT License:
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

const fma_impl = @import("scalar_fma.zig");

const data = @import("scalar_log1pd_data.zig");
const std = @import("std");

const DoubleDouble = struct {
    hi: f64,
    lo: f64,
};

const canonical_nan_bits: u64 = 0x7ff8_0000_0000_0000;
const infinity_bits: u64 = 0x7ff0_0000_0000_0000;
const negative_one_bits: u64 = 0xbff0_0000_0000_0000;
const quiet_nan_bit: u64 = 0x0008_0000_0000_0000;

pub fn log1pd(x: f64) f64 {
    const x_bits: u64 = @bitCast(x);
    const magnitude_bits = x_bits & 0x7fff_ffff_ffff_ffff;
    const absolute_shifted = x_bits << 1;
    var ln_high: f64 = undefined;
    var ln_low: f64 = undefined;
    var error_bound: f64 = undefined;

    if (absolute_shifted < 0x7f60_0000_0000_0000) {
        if (absolute_shifted < 0x7940_0000_0000_0000) {
            if (absolute_shifted == 0) return x;
            return x;
        }

        const x_squared = x * x;
        if (absolute_shifted < 0x7e60_0000_0000_0000) {
            ln_high = x;
            error_bound = 0x1.6p-64 * x;
            if (absolute_shifted < 0x7d43_3600_0000_0000) {
                const coefficients = [_]f64{
                    -0x1.00000000001d1p-1,
                    0x1.55555555558f7p-2,
                };
                ln_low = x_squared * (coefficients[0] + x * coefficients[1]);
            } else {
                const coefficients = [_]f64{
                    -0x1.ffffffffffffdp-2,
                    0x1.5555555555551p-2,
                    -0x1.000000d5555e1p-2,
                    0x1.99999b442f73fp-3,
                };
                ln_low = x_squared *
                    ((coefficients[0] + x * coefficients[1]) +
                        x_squared * (coefficients[2] + x * coefficients[3]));
            }
        } else {
            const coefficients = [_]f64{
                0x1.5555555555555p-2,
                -0x1p-2,
                0x1.9999999999b41p-3,
                -0x1.555555555583bp-3,
                0x1.24924923f39ep-3,
                -0x1.fffffffe42e43p-4,
                0x1.c71c75511d70bp-4,
                -0x1.99999de10510fp-4,
                0x1.7457e81b175f6p-4,
                -0x1.554fb43e54e0fp-4,
                0x1.3ed68744f3d18p-4,
                -0x1.28558ad5a7ac4p-4,
            };
            const x_cubed = x_squared * x;
            const x_fourth = x_squared * x_squared;
            const half_x = -0.5 * x;
            ln_high = fma_impl.fma64(half_x, x, x);
            ln_low = fma_impl.fma64(half_x, x, x - ln_high);
            const polynomial =
                ((coefficients[0] + x * coefficients[1]) +
                    x_squared * (coefficients[2] + x * coefficients[3])) +
                x_fourth *
                    (((coefficients[4] + x * coefficients[5]) +
                        x_squared * (coefficients[6] + x * coefficients[7])) +
                        x_fourth *
                            ((coefficients[8] + x * coefficients[9]) +
                                x_squared * (coefficients[10] + x * coefficients[11])));
            ln_low += x_cubed * polynomial;
            error_bound = x_cubed * 0x1.94p-52;
        }
    } else {
        if (magnitude_bits >= infinity_bits) {
            @branchHint(.unlikely);
            const fraction_bits = magnitude_bits & 0x000f_ffff_ffff_ffff;
            if (fraction_bits != 0) return fromBits(x_bits | quiet_nan_bit);
            if ((x_bits >> 63) == 0) return x;
            return fromBits(canonical_nan_bits);
        }
        if (x_bits == negative_one_bits) {
            @branchHint(.unlikely);
            return fromBits(0xfff0_0000_0000_0000);
        }
        if ((x_bits >> 63) != 0 and x_bits > negative_one_bits) {
            @branchHint(.unlikely);
            return fromBits(canonical_nan_bits);
        }
        const coefficients = [_]f64{
            -0x1.000000000003dp-1,
            0x1.5555555554cf5p-2,
            -0x1.ffffffeca2939p-3,
            0x1.99999a3661724p-3,
            -0x1.555d345bfe6fdp-3,
            0x1.247b887a6e5edp-3,
        };

        var t: f64 = undefined;
        var delta_t: f64 = undefined;
        const signed_x_bits: i64 = @bitCast(x_bits);
        if (signed_x_bits < 0x4340_0000_0000_0000 and x_bits < negative_one_bits) {
            const t_sum = fastTwoSum(1.0, x);
            t = t_sum.hi;
            delta_t = t_sum.lo;
        } else if (x_bits < 0x4690_0000_0000_0000) {
            t = x;
            delta_t = 1.0;
        } else {
            t = x;
            delta_t = 0.0;
        }

        const t_bits_initial: u64 = @bitCast(t);
        const reduction: i64 = @bitCast(t_bits_initial -% 0x3fe6_a000_0000_0000);
        const table_index: usize = @intCast((reduction >> 46) & 0x3f);
        const exponent = reduction >> 52;
        var reciprocal_bits: u64 = @bitCast(data.rf[table_index]);
        if (exponent < 1022) {
            reciprocal_bits -%= exponentOffset(exponent);
        } else {
            reciprocal_bits -%= @as(u64, 1021) << 52;
            const scales = [_]f64{ 0x1p-1, 0x1p-2, 0x1p-3 };
            const scale = scales[@intCast(exponent - 1022)];
            t *= scale;
            delta_t *= scale;
        }

        const reciprocal = fromBits(reciprocal_bits);
        const product_high = reciprocal * t;
        const product_low = fma_impl.fma64(reciprocal, t, -product_high) + reciprocal * delta_t;
        const reduced_sum = fastTwoSum(product_high - 1.0, product_low);
        const reduced_squared = reduced_sum.hi * reduced_sum.hi;
        const reduced_low = reduced_sum.lo + reduced_squared *
            ((coefficients[0] + reduced_sum.hi * coefficients[1]) +
                reduced_squared *
                    ((coefficients[2] + reduced_sum.hi * coefficients[3]) +
                        reduced_squared *
                            (coefficients[4] + reduced_sum.hi * coefficients[5])));

        const exponent_f64: f64 = @floatFromInt(exponent);
        const exponent_high = 0x1.62e42fefa4p-1 * exponent_f64;
        const exponent_low = -0x1.8432a1b0e2634p-43 * exponent_f64;
        ln_high = data.lf[table_index][1] + exponent_high;
        ln_low = data.lf[table_index][0] + exponent_low;
        const combined = fastSum(ln_high, ln_low, reduced_sum.hi, reduced_low);
        ln_high = combined.hi;
        ln_low = combined.lo;
        error_bound = 0x1.ap-65;
    }

    const lower_bound = ln_high + (ln_low - error_bound);
    const upper_bound = ln_high + (ln_low + error_bound);
    if (lower_bound != upper_bound) {
        @branchHint(.unlikely);
        return log1pRefine(x, lower_bound);
    }
    return lower_bound;
}

noinline fn log1pRefine(x: f64, approximation: f64) f64 {
    const near_zero_coefficients = [5][2]f64{
        .{ 0x1.5555555555555p-2, 0x1.5555555555556p-56 },
        .{ -0x1p-2, 0x1.25558eff3c1efp-86 },
        .{ 0x1.999999999999ap-3, -0x1.999999a91d6cap-57 },
        .{ -0x1.5555555555555p-3, -0x1.588aab185593ap-57 },
        .{ 0x1.2492492492492p-3, 0x1.28c490abc528fp-57 },
    };
    const near_zero_low_coefficients = [_]f64{
        -0x1.fffffffffc555p-4,
        0x1.c71c71c7185aap-4,
        -0x1.9999d44449a31p-4,
        0x1.745d51f1817c7p-4,
    };
    const reduced_coefficients = [5][2]f64{
        .{ 1.0, 0.0 },
        .{ -0.5, 0.0 },
        .{ 0x1.5555555555555p-2, 0x1.5555555555555p-56 },
        .{ -0x1p-2, -0x1.80007ed9858b7p-107 },
        .{ 0x1.999999999999ap-3, -0x1.9999999999991p-57 },
    };
    const reduced_low_coefficients = [_]f64{
        -0x1.5555555555555p-3,
        0x1.2492492492492p-3,
        -0x1.0000000073334p-3,
        0x1.c71c71c802f68p-4,
    };

    const x_bits: u64 = @bitCast(x);
    const absolute_shifted = x_bits << 1;
    var result_high: f64 = undefined;
    var result_middle: f64 = undefined;
    var result_low: f64 = undefined;

    if (absolute_shifted < 0x7ea0_0000_0000_0000) {
        if (absolute_shifted < 0x7940_0000_0000_0000) {
            if (absolute_shifted == 0) return x;
            return fma_impl.fma64(fromBits(x_bits & 0x7fff_ffff_ffff_ffff), -0x1p-54, x);
        }

        const x_squared_high = x * x;
        const x_squared_low = fma_impl.fma64(x, x, -x_squared_high);
        const x_cubed = mulDdd(x, x_squared_high, x_squared_low);
        var polynomial_low = x *
            ((near_zero_low_coefficients[0] + x * near_zero_low_coefficients[1]) +
                x_squared_high *
                    (near_zero_low_coefficients[2] + x * near_zero_low_coefficients[3]));
        const polynomial_high = polyDdd(x, &near_zero_coefficients, &polynomial_low);
        const polynomial_product =
            mulDd(polynomial_high, polynomial_low, x_cubed.hi, x_cubed.lo);
        const scaled_x_squared_high = x_squared_high * -0.5;
        const scaled_x_squared_low = x_squared_low * -0.5;
        result_high = x;
        const middle_sum = fastSum(
            scaled_x_squared_high,
            scaled_x_squared_low,
            polynomial_product.hi,
            polynomial_product.lo,
        );
        result_middle = middle_sum.hi;
        result_low = middle_sum.lo;
    } else {
        const rounded = roundEvenFinite(approximation * 0x1.71547652b82fep+16);
        const index = rounded + 34952;
        const index_1: usize = @intCast((index >> 12) & 0xf);
        const index_2: usize = @intCast((index >> 8) & 0xf);
        const index_3: usize = @intCast((index >> 4) & 0xf);
        const index_4: usize = @intCast(index & 0xf);
        const exponent = index >> 16;

        const logarithm = [3]f64{
            (data.ln[0][index_1][0] + data.ln[1][index_2][0]) +
                (data.ln[2][index_3][0] + data.ln[3][index_4][0]),
            (data.ln[0][index_1][1] + data.ln[1][index_2][1]) +
                (data.ln[2][index_3][1] + data.ln[3][index_4][1]),
            (data.ln[0][index_1][2] + data.ln[1][index_2][2]) +
                (data.ln[2][index_3][2] + data.ln[3][index_4][2]),
        };

        var t: f64 = undefined;
        var delta_t: f64 = undefined;
        const signed_x_bits: i64 = @bitCast(x_bits);
        if (signed_x_bits < 0x4690_0000_0000_0000 and x_bits < 0xbfe0_0000_0000_0000) {
            const t_sum = twoSum(1.0, x);
            t = t_sum.hi;
            delta_t = t_sum.lo;
            const delta_bits: u64 = @bitCast(delta_t);
            if ((delta_bits << 1) == 0) delta_t = 0.0;
        } else {
            if (signed_x_bits >= 0x4690_0000_0000_0000) {
                t = x;
            } else {
                t = 1.0 + x;
            }
            delta_t = 0.0;
        }

        var t_bits: u64 = @bitCast(t);
        t_bits -%= exponentOffset(exponent);
        t = fromBits(t_bits);

        const table_product_12 = data.rt[0][index_1] * data.rt[1][index_2];
        const table_product_34 = data.rt[2][index_3] * data.rt[3][index_4];
        const table_high = table_product_12 * table_product_34;
        const table_low = fma_impl.fma64(table_product_12, table_product_34, -table_high);
        const product_high = table_high * t;
        const product_low = fma_impl.fma64(table_high, t, -product_high);
        const low_product_high = table_low * t;
        const low_product_low = fma_impl.fma64(table_low, t, -low_product_high);
        var reduced = fastTwoSum(product_high - 1.0, product_low);
        reduced = fastSum(reduced.hi, reduced.lo, low_product_high, low_product_low);

        if (delta_t != 0.0) {
            var delta_bits: u64 = @bitCast(delta_t);
            delta_bits -%= exponentOffset(exponent);
            delta_t = fromBits(delta_bits);
            const delta_high = table_high * delta_t;
            const delta_low =
                fma_impl.fma64(table_high, delta_t, -delta_high) + table_low * delta_t;
            reduced = fastSum(reduced.hi, reduced.lo, delta_high, delta_low);
        }

        var polynomial_low = reduced.hi *
            ((reduced_low_coefficients[0] + reduced.hi * reduced_low_coefficients[1]) +
                (reduced.hi * reduced.hi) *
                    (reduced_low_coefficients[2] +
                        reduced.hi * reduced_low_coefficients[3]));
        const polynomial_high =
            polyDd(reduced.hi, reduced.lo, &reduced_coefficients, &polynomial_low);
        const polynomial_product =
            mulDd(reduced.hi, reduced.lo, polynomial_high, polynomial_low);

        if (exponent != 0) {
            const exponent_f64: f64 = @floatFromInt(exponent);
            result_high = 0x1.62e42fefa4p-1 * exponent_f64;
            result_middle = -0x1.8432a1b0e28p-43 * exponent_f64;
            result_low = 0x1.cc01f97b57a08p-87 * exponent_f64;
            result_high += logarithm[2];
            const middle_sum = fastSum(
                result_middle,
                result_low,
                logarithm[1],
                logarithm[0],
            );
            result_middle = middle_sum.hi;
            result_low = middle_sum.lo;
        } else {
            result_high = logarithm[2];
            result_middle = logarithm[1];
            result_low = logarithm[0];
        }
        const final_sum = addDoubleDouble(
            result_middle,
            result_low,
            polynomial_product.hi,
            polynomial_product.lo,
        );
        result_middle = final_sum.hi;
        result_low = final_sum.lo;
    }

    var normalized = fastTwoSum(result_high, result_middle);
    result_high = normalized.hi;
    result_middle = normalized.lo;
    normalized = fastTwoSum(result_middle, result_low);
    result_middle = normalized.hi;
    result_low = normalized.lo;

    var middle_bits: u64 = @bitCast(result_middle);
    if ((middle_bits & 0x000f_ffff_ffff_ffff) == 0) {
        const low_bits: u64 = @bitCast(result_low);
        if (((low_bits ^ middle_bits) >> 63) != 0) {
            middle_bits -%= 1;
        } else {
            middle_bits +%= 1;
        }
        result_middle = fromBits(middle_bits);
    }
    return result_high + result_middle;
}

fn fastTwoSum(x: f64, y: f64) DoubleDouble {
    const sum = x + y;
    const rounded_y = sum - x;
    return .{ .hi = sum, .lo = y - rounded_y };
}

fn twoSum(x: f64, y: f64) DoubleDouble {
    const sum = x + y;
    const rounded_y = sum - x;
    return .{
        .hi = sum,
        .lo = (y - rounded_y) + (x + (rounded_y - sum)),
    };
}

fn fastSum(x_high: f64, x_low: f64, y_high: f64, y_low: f64) DoubleDouble {
    const high_sum = fastTwoSum(x_high, y_high);
    return .{ .hi = high_sum.hi, .lo = (x_low + y_low) + high_sum.lo };
}

fn addDoubleDouble(
    x_high: f64,
    x_low: f64,
    y_high: f64,
    y_low: f64,
) DoubleDouble {
    const high_sum = twoSum(x_high, y_high);
    return .{ .hi = high_sum.hi, .lo = (x_low + y_low) + high_sum.lo };
}

fn mulDd(x_high: f64, x_low: f64, y_high: f64, y_low: f64) DoubleDouble {
    const high_product = y_high * x_high;
    return .{
        .hi = high_product,
        .lo = (y_low * x_high + y_high * x_low) +
            fma_impl.fma64(y_high, x_high, -high_product),
    };
}

fn mulDdd(x: f64, y_high: f64, y_low: f64) DoubleDouble {
    const high_product = y_high * x;
    return .{
        .hi = high_product,
        .lo = y_low * x + fma_impl.fma64(y_high, x, -high_product),
    };
}

fn polyDd(
    x_high: f64,
    x_low: f64,
    coefficients: *const [5][2]f64,
    low: *f64,
) f64 {
    var index: usize = 4;
    const initial = fastTwoSum(coefficients[index][0], low.*);
    var high = initial.hi;
    var low_accumulator = coefficients[index][1] + initial.lo;
    while (index != 0) {
        index -= 1;
        const product = mulDd(x_high, x_low, high, low_accumulator);
        const total = fastSum(
            coefficients[index][0],
            coefficients[index][1],
            product.hi,
            product.lo,
        );
        high = total.hi;
        low_accumulator = total.lo;
    }
    low.* = low_accumulator;
    return high;
}

fn polyDdd(x: f64, coefficients: *const [5][2]f64, low: *f64) f64 {
    var index: usize = 4;
    const initial = fastTwoSum(coefficients[index][0], low.*);
    var high = initial.hi;
    var low_accumulator = coefficients[index][1] + initial.lo;
    while (index != 0) {
        index -= 1;
        const product = mulDdd(x, high, low_accumulator);
        const total = fastSum(
            coefficients[index][0],
            coefficients[index][1],
            product.hi,
            product.lo,
        );
        high = total.hi;
        low_accumulator = total.lo;
    }
    low.* = low_accumulator;
    return high;
}

fn roundEvenFinite(value: f64) i64 {
    const value_bits: u64 = @bitCast(value);
    const magnitude_bits = value_bits & 0x7fff_ffff_ffff_ffff;
    const exponent_field: u16 = @intCast((magnitude_bits >> 52) & 0x7ff);
    const exponent = @as(i32, exponent_field) - 1023;

    if (exponent < -1) return 0;
    if (exponent == -1) {
        if (magnitude_bits <= 0x3fe0_0000_0000_0000) return 0;
        return if ((value_bits >> 63) != 0) -1 else 1;
    }
    // Only log1pRefine calls this function. Its scaled approximation has
    // magnitude below 2^27 for every finite input greater than -1.
    std.debug.assert(exponent < 52);

    const significand = (magnitude_bits & 0x000f_ffff_ffff_ffff) |
        0x0010_0000_0000_0000;
    const shift: u6 = @intCast(52 - exponent);
    var integer = significand >> shift;
    const remainder_mask = (@as(u64, 1) << shift) - 1;
    const remainder = significand & remainder_mask;
    const halfway = @as(u64, 1) << (shift - 1);
    if (remainder > halfway or (remainder == halfway and (integer & 1) != 0)) {
        integer += 1;
    }

    const signed_integer: i64 = @intCast(integer);
    return if ((value_bits >> 63) != 0) -signed_integer else signed_integer;
}

fn exponentOffset(exponent: i64) u64 {
    return @as(u64, @bitCast(exponent)) << 52;
}

fn fromBits(value: u64) f64 {
    return @bitCast(value);
}
