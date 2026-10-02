const fma_impl = @import("scalar_fma.zig");

const data = @import("scalar_exp2_data.zig");

const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const quiet_nan_bit_h: u16 = 0x0200;

const sign_mask_f: u32 = 0x8000_0000;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const infinity_bits_f: u32 = 0x7f80_0000;
const quiet_nan_bit_f: u32 = 0x0040_0000;

const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;

// Correctly rounded binary16 exp2, ported from CORE-MATH
// src/binary16/exp2/exp2f16.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
// Copyright (c) 2025 Maxence Ponsardin and Paul Zimmermann.
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

pub fn exp2h(value: f16) f16 {
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;

    if (magnitude > infinity_bits_h) return @bitCast(bits | quiet_nan_bit_h);
    if (magnitude == infinity_bits_h) {
        if ((bits & sign_mask_h) != 0) return @bitCast(@as(u16, 0));
        return value;
    }

    const table_1_index: usize = @intCast(bits >> 5);
    const table_2_index: usize = @intCast(((bits >> 10) << 5) | (bits & 0x001f));
    const table_1: f32 = @bitCast(data.t1[table_1_index]);
    const table_2: f32 = @bitCast(data.t2[table_2_index]);
    return @floatCast(table_1 * table_2);
}

noinline fn exp2fSpecial(value: f32, bits: u32) f32 {
    const magnitude = bits & magnitude_mask_f;
    if (magnitude > infinity_bits_f) return @bitCast(bits | quiet_nan_bit_f);
    if (magnitude == infinity_bits_f) {
        if ((bits & sign_mask_f) != 0) return 0.0;
        return @bitCast(infinity_bits_f);
    }
    if ((bits & sign_mask_f) == 0) return @bitCast(infinity_bits_f);

    const value_wide: f64 = value;
    var result = 0x1p-149 + (value_wide + 149.0) * 0x1p-150;
    if (result < 0x1p-151) result = 0x1p-151;
    return @floatCast(result);
}

noinline fn exp2fAccurate(scale: f64, rounded_value: f64, rounded_squared: f64) f32 {
    @setFloatMode(.strict);
    const accurate_coefficients = [_]f64{
        0x1.62e42fefa39efp-1,
        0x1.ebfbdff82c58fp-3,
        0x1.c6b08d702e0edp-5,
        0x1.3b2ab6fb92e5ep-7,
        0x1.5d886e6d54203p-10,
        0x1.430976b8ce6efp-13,
    };
    const result = scale + (scale * rounded_value) *
        ((accurate_coefficients[0] + rounded_value * accurate_coefficients[1]) +
            rounded_squared * ((accurate_coefficients[2] + rounded_value * accurate_coefficients[3]) +
                rounded_squared * (accurate_coefficients[4] + rounded_value * accurate_coefficients[5])));
    return @floatCast(result);
}

fn exp2fFinite(value: f32, bits: u32) f32 {
    @setFloatMode(.strict);

    const unsigned_magnitude = bits << 1;
    if (unsigned_magnitude >= 0x8600_0000 or unsigned_magnitude < 0x6500_0000) {
        if (unsigned_magnitude < 0x6500_0000) return 1.0;
        if (!(bits >= 0xc300_0000 and bits < 0xc315_0000)) return exp2fSpecial(value, bits);
    }

    // The binary32 bias supplies the same 1/64 grid as a binary64 reduction.
    const scale_biased: f32 = value + @as(f32, 0x1.8p17);
    const scale_biased_bits: u32 = @bitCast(scale_biased);
    const table_index: u6 = @truncate(scale_biased_bits);
    const scale_bits = data.exp2f_scale_table[table_index] +%
        (@as(u64, scale_biased_bits) << 46);
    const scale: f64 = @bitCast(scale_bits);

    // Both subtractions are exact in binary32. A nonzero nearest
    // grid point and value satisfy Sterbenz; a zero grid leaves value.
    const rounded_grid: f32 = scale_biased - @as(f32, 0x1.8p17);
    const rounded_value_narrow: f32 = value - rounded_grid;
    const rounded_value: f64 = rounded_value_narrow;
    const fast_coefficients = [_]f64{
        1.0,
        0x1.62e42fef4c4e7p-1,
        0x1.ebfd1b232f475p-3,
        0x1.c6b19384ecd93p-5,
    };
    var polynomial = fma_impl.fma64(fast_coefficients[3], rounded_value, fast_coefficients[2]);
    polynomial = fma_impl.fma64(polynomial, rounded_value, fast_coefficients[1]);
    polynomial = fma_impl.fma64(polynomial, rounded_value, fast_coefficients[0]);
    const result = scale * polynomial;
    // Horner changes the original evaluation by less than 16u, u=2^-53.
    // Its total relative error remains below 2^-32. Every binary32
    // midpoint has 28 low binary64 bits clear; 2^21 adjacent encodings
    // enclose the error on either side, including binade transitions.
    const result_bits: u64 = @bitCast(result);
    if (((result_bits +% 0x0020_0000) & 0x0fff_ffff) > 0x0040_0000) return @floatCast(result);
    if (unsigned_magnitude <= 0x79e7_526e) {
        if (bits == 0x3b42_9d37) return @as(f32, 0x1.00870ap+0) - @as(f32, 0x1p-25);
        if (bits == 0xbcf3_a937) return @as(f32, 0x1.f58d62p-1) - @as(f32, 0x1p-26);
        if (bits == 0xb8d3_d026) return @as(f32, 0x1.fff6d2p-1) + @as(f32, 0x1p-26);
    }
    return exp2fAccurate(scale, rounded_value, rounded_value * rounded_value);
}

// Correctly rounded binary32 exp2, ported from CORE-MATH
// src/binary32/exp2/exp2f.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
// Copyright (c) 2023-2025 Alexei Sibidanov.
// This code is derived from CORE-MATH src/binary32/exp2/exp2f.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT License.

pub fn exp2f(value: f32) f32 {
    const bits: u32 = @bitCast(value);
    return exp2fFinite(value, bits);
}

// Correctly rounded binary64 exp2, ported from CORE-MATH
// src/binary64/exp2/exp2.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
// Copyright (c) 2023-2025 Alexei Sibidanov. Used under the MIT License.
const Exp2dDoubleDouble = data.DoubleDouble;

fn exp2dRoundEvenFinite(value: f64) i64 {
    @setFloatMode(.strict);

    // All callers prove |value| < 2^23. Under RNE, the fixed-exponent
    // addition rounds to the nearest even integer and the bit delta is its
    // signed two's-complement representation.
    const rounding_constant: f64 = 0x1.8p52;
    const biased = value + rounding_constant;
    const biased_bits: u64 = @bitCast(biased);
    const rounding_bits: u64 = @bitCast(rounding_constant);
    return @bitCast(biased_bits -% rounding_bits);
}

fn exp2dFastTwoSum(x: f64, y: f64) Exp2dDoubleDouble {
    const sum = x + y;
    const rounded_x = sum - x;
    return .{ .high = sum, .low = y - rounded_x };
}

fn exp2dFastSum(x_high: f64, x_low: f64, y_high: f64, y_low: f64) Exp2dDoubleDouble {
    const sum = exp2dFastTwoSum(x_high, y_high);
    return .{ .high = sum.high, .low = (x_low + y_low) + sum.low };
}

fn exp2dMultiplyDoubleDouble(
    x_high: f64,
    x_low: f64,
    c_high: f64,
    c_low: f64,
) Exp2dDoubleDouble {
    const high = c_high * x_high;
    const low = (c_low * x_high + c_high * x_low) + fma_impl.fma64(c_high, x_high, -high);
    return .{ .high = high, .low = low };
}

fn exp2dMultiplyDouble(
    x_high: f64,
    c_high: f64,
    c_low: f64,
) Exp2dDoubleDouble {
    const high = c_high * x_high;
    const low = c_low * x_high + fma_impl.fma64(c_high, x_high, -high);
    return .{ .high = high, .low = low };
}

fn exp2dPolynomial(value: f64, coefficients: []const Exp2dDoubleDouble) Exp2dDoubleDouble {
    var index = coefficients.len - 1;
    var high = coefficients[index].high;
    var low = coefficients[index].low;
    while (index != 0) {
        index -= 1;
        const product = exp2dMultiplyDouble(value, high, low);
        high = product.high;
        low = product.low;
        const sum_high = high + coefficients[index].high;
        const sum_low = (coefficients[index].high - sum_high) + high;
        high = sum_high;
        low += sum_low + coefficients[index].low;
    }
    return .{ .high = high, .low = low };
}

fn exp2dScaleNormal(value: f64, exponent: i64) f64 {
    const value_bits: u64 = @bitCast(value);
    const exponent_bits: u64 = @as(u64, @bitCast(exponent)) << 52;
    return @bitCast(value_bits +% exponent_bits);
}

fn exp2dToDenormal(value: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    return @bitCast(value_bits & 0x000f_ffff_ffff_ffff);
}

fn exp2dDatabaseLookup(value: f64, fallback: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    var lower: isize = 0;
    var upper: isize = @intCast(data.exp2d_database.len - 1);
    while (lower <= upper) {
        const middle: isize = @divTrunc(lower + upper, 2);
        const middle_index: usize = @intCast(middle);
        const database_bits: u64 = @bitCast(data.exp2d_database[middle_index]);
        if (database_bits < value_bits) {
            lower = middle + 1;
        } else if (database_bits > value_bits) {
            upper = middle - 1;
        } else {
            const rounding_patterns = [2]u64{ 0x3b216fbd5fd7665f, 0x0034c797 };
            const sign_pattern: u64 = 8677191773140;
            const middle_shift: u6 = @intCast(middle);
            const pattern = (rounding_patterns[middle_index >> 5] >> @intCast((middle * 2) & 63)) & 3;
            const delta_bits: u64 = @as(u64, 0x3c90 | ((sign_pattern >> middle_shift) << 15)) << 48;
            const delta: f64 = @bitCast(delta_bits);
            const fallback_bits: u64 = @bitCast(fallback);
            var offset: i64 = -1;
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

noinline fn exp2dAccurate(value: f64) f64 {
    @setFloatMode(.strict);

    const value_bits: u64 = @bitCast(value);
    const scaled_value = 4096.0 * value;
    const rounded_integer = exp2dRoundEvenFinite(scaled_value);
    const reduced = scaled_value - @as(f64, @floatFromInt(rounded_integer));
    const table_0_index: usize = @intCast((rounded_integer >> 6) & 0x3f);
    const table_1_index: usize = @intCast(rounded_integer & 0x3f);
    const exponent = rounded_integer >> 12;
    const table_0 = data.exp2d_table_0[table_0_index];
    const table_1 = data.exp2d_table_1[table_1_index];
    const table_product = exp2dMultiplyDoubleDouble(table_0.high, table_0.low, table_1.high, table_1.low);
    const coefficients = [_]Exp2dDoubleDouble{
        .{ .high = 0x1.62e42fefa39efp-13, .low = 0x1.abc9e3b39873ep-68 },
        .{ .high = 0x1.ebfbdff82c58fp-27, .low = -0x1.5e43a53e4495p-81 },
        .{ .high = 0x1.c6b08d704a0cp-41, .low = -0x1.d3a15710d3d83p-95 },
        .{ .high = 0x1.3b2ab6fba4e77p-55, .low = 0x1.4dd5d2a5e025ap-110 },
        .{ .high = 0x1.5d87fe7a66459p-70, .low = -0x1.dc47e47beb9ddp-124 },
        .{ .high = 0x1.430912f9fb79dp-85, .low = -0x1.4fcd51fcb764p-139 },
    };
    var polynomial = exp2dPolynomial(reduced, &coefficients);
    polynomial = exp2dMultiplyDouble(reduced, polynomial.high, polynomial.low);

    if (value_bits <= 0xc08f_f000_0000_0000) {
        if (-0x1.71547652b82fep-54 <= value and value <= 0x1.71547652b82fdp-53) {
            return fma_impl.fma64(value, 0.5, 1.0);
        }
        if ((rounded_integer & 0xfff) == 0) {
            const sum = exp2dFastTwoSum(table_product.high, polynomial.high);
            const low_sum = exp2dFastTwoSum(sum.low, polynomial.low);
            var result: Exp2dDoubleDouble = .{ .high = sum.high, .low = low_sum.high };
            var low_bits: u64 = @bitCast(result.low);
            if ((low_bits & 0x000f_ffff_ffff_ffff) == 0 and ((low_bits >> 52) & 0x7ff) != 0) {
                const error_bits: u64 = @bitCast(low_sum.low);
                const delta: u64 = if ((low_bits ^ error_bits) & sign_mask_d == 0)
                    1
                else
                    ~@as(u64, 0);
                low_bits +%= delta;
                result.low = @bitCast(low_bits);
            }
            result = exp2dFastTwoSum(result.high, result.low);
            const result_low_bits: u64 = @bitCast(result.low);
            if (((result_low_bits +% 2) & 0x000f_ffff_ffff_ffff) <= 2) {
                result.high = exp2dDatabaseLookup(value, result.high);
            }
            return exp2dScaleNormal(result.high, exponent);
        }

        const product = exp2dMultiplyDoubleDouble(
            polynomial.high,
            polynomial.low,
            table_product.high,
            table_product.low,
        );
        var result = exp2dFastSum(
            table_product.high,
            table_product.low,
            product.high,
            product.low,
        );
        result = exp2dFastTwoSum(result.high, result.low);
        const result_low_bits: u64 = @bitCast(result.low);
        if (((result_low_bits +% 2) & 0x000f_ffff_ffff_ffff) <= 2) {
            result.high = exp2dDatabaseLookup(value, result.high);
        }
        return exp2dScaleNormal(result.high, exponent);
    }

    const scale_bits = @as(u64, @bitCast(@as(i64, 1) - exponent)) << 52;
    const scale: f64 = @bitCast(scale_bits);
    const product = exp2dMultiplyDoubleDouble(
        polynomial.high,
        polynomial.low,
        table_product.high,
        table_product.low,
    );
    var result = exp2dFastSum(
        table_product.high,
        table_product.low,
        product.high,
        product.low,
    );
    const scaled_sum = exp2dFastTwoSum(scale, result.high);
    result.high = scaled_sum.high;
    result.low += scaled_sum.low;
    return exp2dToDenormal(result.high + result.low);
}

pub fn exp2d(value: f64) f64 {
    @setFloatMode(.strict);

    const bits: u64 = @bitCast(value);
    const unsigned_magnitude = bits << 1;
    if (unsigned_magnitude == 0) return 1.0;
    if (unsigned_magnitude >= 0x8120_0000_0000_0000) {
        if (unsigned_magnitude > 0xffe0_0000_0000_0000) return @bitCast(bits | quiet_nan_bit_d);
        if (unsigned_magnitude == 0xffe0_0000_0000_0000) {
            if ((bits & sign_mask_d) != 0) return 0.0;
            return value;
        }
        if ((bits & sign_mask_d) != 0) {
            if (bits >= 0xc090_cc00_0000_0000) return 0.0;
        } else {
            return @bitCast(infinity_bits_d);
        }
    }

    if (unsigned_magnitude <= 0x792e_2a8e_ca57_05fc) return 1.0;

    const scaled_value = 4096.0 * value;
    const biased = scaled_value + 0x1.8p52;
    const biased_bits: u64 = @bitCast(biased);
    const rounded_integer: i64 = @bitCast(biased_bits -% 0x4338_0000_0000_0000);
    const reduced = scaled_value - (biased - 0x1.8p52);
    const table_0_index: usize = @intCast((rounded_integer >> 6) & 0x3f);
    const table_1_index: usize = @intCast(rounded_integer & 0x3f);
    const exponent = rounded_integer >> 12;
    const table_0 = data.exp2d_table_0[table_0_index];
    const table_1 = data.exp2d_table_1[table_1_index];
    const table_high = table_0.high * table_1.high;
    const table_error = fma_impl.fma64(table_0.high, table_1.high, -table_high);
    const table_low = fma_impl.fma64(
        table_0.low,
        table_1.high,
        fma_impl.fma64(table_0.high, table_1.low, table_error),
    );
    const table_product: Exp2dDoubleDouble = .{ .high = table_high, .low = table_low };
    const reduced_squared = reduced * reduced;
    const linear = fma_impl.fma64(reduced, 0x1.ebfbdff82c58fp-27, 0x1.62e42fefa39efp-13);
    const cubic = fma_impl.fma64(reduced, 0x1.3b2ab6fdda001p-55, 0x1.c6b08d73b3e01p-41);
    const polynomial = fma_impl.fma64(reduced_squared, cubic, linear);
    var high = table_product.high;
    const table_times_reduced = table_product.high * reduced;
    var low = table_product.low + table_times_reduced * polynomial;
    const error_bound: f64 = 1.64e-19;

    if (bits <= 0xc08f_f000_0000_0000) {
        // Integer inputs reduce to 1+0 and pass this same enclosure.
        const upper = high + (low + error_bound);
        const lower = high + (low - error_bound);
        if (upper != lower) return exp2dAccurate(value);
        high = lower;
        return exp2dScaleNormal(high, exponent);
    }

    const scale_bits = @as(u64, @bitCast(@as(i64, 1) - exponent)) << 52;
    const scale: f64 = @bitCast(scale_bits);
    const sum = exp2dFastTwoSum(scale, high);
    high = sum.high;
    low += sum.low;
    // Cover low-part addition and endpoint rounding after the coarse
    // subnormal scale has enlarged the low part.
    const denormal_error = fma_impl.fma64(@abs(low), 0x1p-51, error_bound);
    const upper = high + (low + denormal_error);
    const lower = high + (low - denormal_error);
    if (upper != lower) return exp2dAccurate(value);
    high = lower;
    return exp2dToDenormal(high);
}
