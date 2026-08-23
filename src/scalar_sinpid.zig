// Correctly rounded half-revolution sine for binary64 nearest-even rounding.
//
// Copyright (c) 2023-2025 Alexei Sibidanov.
//
// This code is derived from CORE-MATH src/binary64/sinpi/sinpi.c at revision
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

const data = @import("scalar_sinpid_data.zig");

const sign_mask: u64 = 0x8000_0000_0000_0000;
const magnitude_mask: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits: u64 = infinity_bits | quiet_nan_bit;

const SinCos = struct {
    sin_high: f64,
    sin_low: f64,
    cos_high: f64,
    cos_low: f64,
};

fn copySign(value: f64, bits: u64) f64 {
    return @bitCast((@as(u64, @bitCast(value)) & magnitude_mask) | (bits & sign_mask));
}

fn fastTwoSum(x: f64, y: f64, error_value: *f64) f64 {
    const sum = x + y;
    const rounded_x = sum - x;
    error_value.* = y - rounded_x;
    return sum;
}

fn addDoubleDouble(x_high: f64, x_low: f64, c_high: f64, c_low: f64, low: *f64) f64 {
    const sum = x_high + c_high;
    const difference = sum - x_high;
    low.* = ((c_high - difference) + (x_high + (difference - sum))) + (x_low + c_low);
    return sum;
}

fn mulDoubleDoubleAcc(
    x_high: f64,
    x_low: f64,
    c_high: f64,
    c_low: f64,
    low: *f64,
) f64 {
    const high_low = c_high * x_low;
    const low_high = c_low * x_high;
    const high_high = c_high * x_high;
    var high_error = @mulAdd(f64, c_high, x_high, -high_high);
    high_error += low_high + high_low;
    return fastTwoSum(high_high, high_error, low);
}

fn mulDoubleDouble(x_high: f64, x_low: f64, c_high: f64, low: *f64) f64 {
    const high_low = c_high * x_low;
    const high_high = c_high * x_high;
    const high_error = @mulAdd(f64, c_high, x_high, -high_high);
    const middle = high_error + high_low;
    const result = high_high + middle;
    low.* = (high_high - result) + middle;
    return result;
}

fn polyDoubleDouble(
    x_high: f64,
    x_low: f64,
    coefficients: []const [2]f64,
    low: *f64,
) f64 {
    var index: usize = coefficients.len - 1;
    var high = coefficients[index][0] + low.*;
    var low_part = ((coefficients[index][0] - high) + low.*) + coefficients[index][1];
    while (index != 0) {
        index -= 1;
        high = mulDoubleDoubleAcc(x_high, x_low, high, low_part, &low_part);
        const temporary = high + coefficients[index][0];
        const temporary_low = (coefficients[index][0] - temporary) + high;
        high = temporary;
        low_part += temporary_low + coefficients[index][1];
    }
    low.* = low_part;
    return high;
}

fn sinpiZero(value: f64) f64 {
    const square = value * value;
    const square_error = @mulAdd(f64, value, value, -square);
    const cube = square * value;
    const cube_error = @mulAdd(f64, square, value, -cube) + square_error * value;
    const coefficients = [_][2]f64{
        .{ -0x1.4abbce625be53p+2, 0x1.05511c68477bep-52 },
        .{ 0x1.466bc6775aae2p+1, -0x1.6dc0cbefae1dap-54 },
        .{ -0x1.32d2cce62bd86p-1, 0x1.066bd54973829p-55 },
        .{ 0x1.50783487ee781p-4, 0x1.832989f39a743p-58 },
    };
    const low_coefficients = [_]f64{
        -0x1.e3074fde861fp-8,
        0x1.e8f4344534da6p-12,
        -0x1.6f9cd7b8cb9dbp-16,
    };
    var low_part = square * (low_coefficients[0] + square *
        (low_coefficients[1] + square * low_coefficients[2]));
    var high = polyDoubleDouble(square, square_error, &coefficients, &low_part);
    high = mulDoubleDoubleAcc(high, low_part, cube, cube_error, &low_part);

    const pi_high = 0x1.92p+1;
    const pi_middle = 0x1.fb54442d1846ap-11;
    const pi_low = -0x1.d9cceba3f91f2p-65;
    var rounded = pi_high * value;
    var rounded_bits: u64 = @bitCast(rounded);
    rounded_bits &= 0x000f_ffff_ffff_ffff;
    rounded_bits +%= @as(u64, 85) << 51;
    rounded = (rounded + @as(f64, @bitCast(rounded_bits))) - @as(f64, @bitCast(rounded_bits));
    const rounded_low = @mulAdd(f64, pi_high, value, -rounded);
    var middle = pi_middle * value;
    var middle_low = @mulAdd(f64, pi_middle, value, -middle) + pi_low * value;
    middle = addDoubleDouble(middle, middle_low, rounded_low, 0.0, &middle_low);
    middle = addDoubleDouble(middle, middle_low, high, low_part, &low_part);
    rounded = fastTwoSum(rounded, middle, &middle);
    middle = fastTwoSum(middle, low_part, &low_part);

    const middle_bits: u64 = @bitCast(middle);
    if ((middle_bits & 0x000f_ffff_ffff_ffff) == 0) {
        const low_bits: u64 = @bitCast(low_part);
        var adjusted_bits = middle_bits;
        if (((low_bits ^ middle_bits) >> 63) != 0) adjusted_bits -%= 1 else adjusted_bits +%= 1;
        middle = @bitCast(adjusted_bits);
    }
    return rounded + middle;
}

fn sinpiRefine(index_input: i32, z: f64) f64 {
    const x = z * 0x1p-63;
    const square = x * x;
    const square_error = @mulAdd(f64, x, x, -square);
    const sine_coefficients = [_][2]f64{
        .{ 0x1.921fb54442d18p+1, 0x1.1a62633145c06p-53 },
        .{ -0x1.4abbce625be53p-22, 0x1.05511cbc65743p-76 },
        .{ 0x1.466bc6775aae1p-47, -0x1.9c3c168d990ap-114 },
    };
    const cosine_coefficients = [_][2]f64{
        .{ -0x1.3bd3cc9be45dep-22, -0x1.692b71366cc04p-76 },
        .{ 0x1.03c1f081b5ac4p-46, -0x1.32b33fda9113cp-100 },
    };
    var sine_low = -0x1.32d2cc920dcb4p-73 * square;
    var sine_high = polyDoubleDouble(square, square_error, &sine_coefficients, &sine_low);
    sine_high = mulDoubleDouble(sine_high, sine_low, x * 0x1p-12, &sine_low);
    var cosine_low = square * (-0x1.55d3c7e3cbff9p-72 + 0x1.e1f50604fa0ffp-99 * square);
    var cosine_high = polyDoubleDouble(square, square_error, &cosine_coefficients, &cosine_low);
    cosine_high = mulDoubleDoubleAcc(cosine_high, cosine_low, square, square_error, &cosine_low);

    const table = sincosN2(index_input);
    var sine_correction_low: f64 = undefined;
    const sine_correction_high = mulDoubleDoubleAcc(
        cosine_high,
        cosine_low,
        table.sin_high,
        table.sin_low,
        &sine_correction_low,
    );
    var cosine_correction_low: f64 = undefined;
    const cosine_correction_high = mulDoubleDoubleAcc(
        sine_high,
        sine_low,
        table.cos_high,
        table.cos_low,
        &cosine_correction_low,
    );
    var result_low: f64 = undefined;
    var result_high = fastTwoSum(
        cosine_correction_high,
        sine_correction_high,
        &result_low,
    );
    result_low += sine_correction_low + cosine_correction_low;
    var result_low_error: f64 = undefined;
    result_high = fastTwoSum(table.sin_high, result_high, &result_low_error);
    result_low = table.sin_low + result_low + result_low_error;

    const result_bits: u64 = @bitCast(result_low);
    if ((result_bits | (@as(u64, 0xfff) << 52)) == ~@as(u64, 0) or
        (result_bits << 12) == 0)
    {
        const hard_cases = [_]struct { index: i32, x: f64, result: f64, delta: f64 }{
            .{ .index = 76, .x = -0x1.276b3fef466p-2, .result = 0x1.db8a79a80c3a0p-4, .delta = 0x1p-110 },
            .{ .index = 108, .x = -0x1.33caea0f24cp-2, .result = 0x1.5146c0bc45bcep-3, .delta = 0x1p-109 },
            .{ .index = 490, .x = -0x1.8a1e8a3e82cp-1, .result = 0x1.5d6561936b699p-1, .delta = -0x1p-55 },
            .{ .index = 903, .x = -0x1.bdd02d1ad60p-2, .result = 0x1.f72c906962631p-1, .delta = 0x1p-55 },
        };
        const sign_multiplier: f64 = if (index_input > 2048) -1.0 else 1.0;
        const reduced_index = index_input & 0x7ff;
        for (hard_cases) |hard_case| {
            if ((x == hard_case.x and reduced_index == hard_case.index) or
                (x == -hard_case.x and reduced_index == 2048 - hard_case.index))
            {
                return sign_multiplier * hard_case.result + sign_multiplier * hard_case.delta;
            }
        }
    }
    return result_high + result_low;
}

/// Return correctly rounded sinpi for a binary64 value.
///
/// The caller must keep IEEE round-to-nearest, ties-to-even active.
pub fn sinpid(value: f64) f64 {
    @setFloatMode(.strict);

    const bits: u64 = @bitCast(value);
    const magnitude = bits & magnitude_mask;
    if (magnitude == 0) return value;

    const exponent: i32 = @intCast(magnitude >> 52);
    const significand_bits = (bits & 0x000f_ffff_ffff_ffff) | (@as(u64, 1) << 52);
    const sign: i64 = if ((bits & sign_mask) == 0) 0 else -1;
    var mantissa: i64 = @intCast(significand_bits);
    mantissa = (mantissa ^ sign) - sign;
    const shift: i32 = 1063 - exponent;

    if (shift < 0) {
        if (exponent == 0x7ff) {
            if ((bits << 12) == 0) return @bitCast(canonical_nan_bits);
            return @bitCast(bits | quiet_nan_bit);
        }
        const reduction_shift = -shift - 1;
        if (reduction_shift > 10) return copySign(0.0, bits);
        const mantissa_bits: u64 = @bitCast(mantissa);
        const index = (mantissa_bits << @intCast(reduction_shift)) & 0xfff;
        if ((index & 0x7ff) == 0) return copySign(0.0, bits);
        const table = sincosN(@intCast(index));
        return table.sin_high + table.sin_low;
    }

    if (magnitude <= 0x3fa2_0000_0000_0000) {
        const pi_high = 0x1.921fb54442d18p+1;
        const pi_low = 0x1.1a62633145c07p-53;
        var high: f64 = undefined;
        var low: f64 = undefined;
        const absolute_value: f64 = @bitCast(magnitude);
        if (absolute_value < 0x1p-54) {
            if (absolute_value < 0x1p-970) {
                const scaled = value * 0x1p106;
                high = pi_high * scaled;
                low = @mulAdd(f64, pi_high, scaled, -high) + pi_low * scaled;
                const sum = high + low;
                const scaled_sum = sum * 0x1p-106;
                const restored = scaled_sum * 0x1p106;
                return @mulAdd(f64, (high - restored) + low, 0x1p-106, scaled_sum);
            }
            high = pi_high * value;
            low = @mulAdd(f64, pi_high, value, -high) + pi_low * value;
            return high + low;
        }
        high = pi_high * value;
        low = @mulAdd(f64, pi_high, value, -high) + pi_low * value;
        const square = value * value;
        const cube = square * value;
        const fourth = square * square;
        const epsilon = value * (square * 0x1p-47 + 0x1p-102);
        const coefficients = [_]f64{
            -0x1.4abbce625be51p+2,
            0x1.466bc67754b46p+1,
            -0x1.32d2cc12a51f4p-1,
            0x1.5060540058476p-4,
        };
        low += cube * ((coefficients[0] + square * coefficients[1]) +
            fourth * (coefficients[2] + square * coefficients[3]));
        const lower = (low - epsilon) + high;
        const upper = (low + epsilon) + high;
        if (lower == upper) return lower;
        return sinpiZero(value);
    }

    const integer_shift = exponent - 1011;
    if (integer_shift >= 0 and
        (significand_bits << @intCast(integer_shift + 1)) == 0)
    {
        if ((significand_bits << @intCast(integer_shift)) == 0) {
            return copySign(0.0, bits);
        }
        const parity = (significand_bits << @intCast(integer_shift - 1)) >> 63;
        const signed_one = copySign(1.0, bits);
        if (parity == 0) return signed_one;
        return -signed_one;
    }
    const mantissa_bits: u64 = @bitCast(mantissa);
    const quotient = ((mantissa_bits >> @intCast(shift)) & 0x1fff);
    const table_index: i32 = @intCast((quotient + 1) >> 1);
    const k_bits = mantissa_bits << @intCast(exponent - 1000);
    const k: i64 = @bitCast(k_bits);
    const z: f64 = @floatFromInt(k);
    const square = z * z;
    const sine_correction = 0x1.921fb54442d18p-74 + square *
        (-0x1.4abbce625be51p-223 + square * 0x1.466bc6044ba16p-374);
    const cosine_correction = -0x1.3bd3cc9be45dbp-148 + square *
        0x1.03c1f00186416p-298;
    const table = sincosN(table_index);
    const error_bound = 5.5e-19;
    const result = table.sin_low + table.sin_high * (square * cosine_correction) +
        table.cos_high * (z * sine_correction);
    const lower = (result - error_bound) + table.sin_high;
    const upper = (result + error_bound) + table.sin_high;
    if (lower == upper) return lower;
    return sinpiRefine(table_index, z);
}

fn sincosN(input: i32) SinCos {
    var reduced = input & 0x3ff;
    if ((input & 0x400) != 0) reduced = 0x400 - reduced;
    const sine_index: usize = @intCast(reduced >> 5);
    const cosine_index: usize = @intCast(0x20 - (reduced >> 5));
    const low_index: usize = @intCast(reduced & 0x1f);
    const sine_sign = ((input >> 11) & 1) != 0;
    const cosine_sign = (((input + 1024) >> 11) & 1) != 0;

    const sine_base_high = data.sincosn_sn[sine_index][0];
    const sine_base_low = data.sincosn_sn[sine_index][1];
    const cosine_base_high = data.sincosn_sn[cosine_index][0];
    const cosine_base_low = data.sincosn_sn[cosine_index][1];
    const sine_delta_high = data.sincosn_sm[low_index][0];
    const sine_delta_low = data.sincosn_sm[low_index][1];
    const cosine_delta_high = data.sincosn_cm[low_index][0];
    const cosine_delta_low = data.sincosn_cm[low_index][1];

    const sine_base = sine_base_high + sine_base_low;
    const cosine_base = cosine_base_high + cosine_base_low;
    const cosine_high = cosine_base_high * cosine_delta_high - sine_base_high * sine_delta_high;
    const cosine_low = cosine_delta_high * cosine_base_low - sine_delta_high * sine_base_low +
        cosine_base * cosine_delta_low - sine_base * sine_delta_low;
    const sine_high = sine_base_high * cosine_delta_high + cosine_base_high * sine_delta_high;
    const sine_low = sine_delta_high * cosine_base_low + cosine_delta_high * sine_base_low +
        cosine_base * sine_delta_low + sine_base * cosine_delta_low;

    var result_cosine_low: f64 = undefined;
    const result_cosine_high = fastTwoSum(cosine_high, cosine_low, &result_cosine_low);
    var result_sine_low: f64 = undefined;
    const result_sine_high = fastTwoSum(sine_high, sine_low, &result_sine_low);
    return .{
        .sin_high = if (sine_sign) -result_sine_high else result_sine_high,
        .sin_low = if (sine_sign) -result_sine_low else result_sine_low,
        .cos_high = if (cosine_sign) -result_cosine_high else result_cosine_high,
        .cos_low = if (cosine_sign) -result_cosine_low else result_cosine_low,
    };
}

fn sincosN2(input: i32) SinCos {
    var reduced = input & 0x3ff;
    if ((input & 0x400) != 0) reduced = 0x400 - reduced;
    const sine_index: usize = @intCast(reduced >> 5);
    const cosine_index: usize = @intCast(0x20 - (reduced >> 5));
    const low_index: usize = @intCast(reduced & 0x1f);
    const sine_sign = ((input >> 11) & 1) != 0;
    const cosine_sign = (((input + 1024) >> 11) & 1) != 0;

    const sine_base_high = data.sincosn2_sn[sine_index][0];
    const sine_base_low = data.sincosn2_sn[sine_index][1];
    const cosine_base_high = data.sincosn2_sn[cosine_index][0];
    const cosine_base_low = data.sincosn2_sn[cosine_index][1];
    const sine_delta_high = data.sincosn2_sm[low_index][0];
    const sine_delta_low = data.sincosn2_sm[low_index][1];
    const cosine_delta_high = data.sincosn2_cm[low_index][0];
    const cosine_delta_low = data.sincosn2_cm[low_index][1];

    var cosine_product_low: f64 = undefined;
    const cosine_product_high = mulDoubleDoubleAcc(
        cosine_delta_high,
        cosine_delta_low,
        cosine_base_high,
        cosine_base_low,
        &cosine_product_low,
    );
    var sine_product_low: f64 = undefined;
    const sine_product_high = mulDoubleDoubleAcc(
        sine_delta_high,
        sine_delta_low,
        sine_base_high,
        sine_base_low,
        &sine_product_low,
    );
    var sine_cosine_low: f64 = undefined;
    const sine_cosine_high = mulDoubleDoubleAcc(
        cosine_delta_high,
        cosine_delta_low,
        sine_base_high,
        sine_base_low,
        &sine_cosine_low,
    );
    var cosine_sine_low: f64 = undefined;
    const cosine_sine_high = mulDoubleDoubleAcc(
        sine_delta_high,
        sine_delta_low,
        cosine_base_high,
        cosine_base_low,
        &cosine_sine_low,
    );

    var result_cosine_low: f64 = undefined;
    var result_cosine_high = fastTwoSum(sine_product_high, cosine_product_high, &result_cosine_low);
    result_cosine_low += cosine_product_low + sine_product_low;
    var result_sine_low: f64 = undefined;
    var result_sine_high = fastTwoSum(-cosine_sine_high, sine_cosine_high, &result_sine_low);
    result_sine_low += sine_cosine_low - cosine_sine_low;
    var correction: f64 = undefined;
    result_cosine_high = fastTwoSum(cosine_base_high, -result_cosine_high, &correction);
    result_cosine_low = cosine_base_low - result_cosine_low + correction;
    result_sine_high = fastTwoSum(sine_base_high, -result_sine_high, &correction);
    result_sine_low = sine_base_low - result_sine_low + correction;
    return .{
        .sin_high = if (sine_sign) -result_sine_high else result_sine_high,
        .sin_low = if (sine_sign) -result_sine_low else result_sine_low,
        .cos_high = if (cosine_sign) -result_cosine_high else result_cosine_high,
        .cos_low = if (cosine_sign) -result_cosine_low else result_cosine_low,
    };
}
