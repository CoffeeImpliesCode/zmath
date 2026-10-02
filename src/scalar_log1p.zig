// Correctly rounded biased logarithm functions for binary16 and binary32 values.
//
// Copyright (c) 2025 Paul Zimmermann.
// Copyright (c) 2023-2026 Alexei Sibidanov <sibid@uvic.ca>.
//
// This code is derived from CORE-MATH src/binary16/log1p/log1pf16.c and
// src/binary32/log1p/log1pf.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
//
// This file is part of the CORE-MATH project
// (https://core-math.gitlabpages.inria.fr/).
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

const data = @import("scalar_log1ph_data.zig");

const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const quiet_nan_bit_h: u16 = 0x0200;
const canonical_nan_bits_h: u16 = 0x7e00;
const minus_infinity_bits_h: u16 = 0xfc00;

const sign_mask_f: u32 = 0x8000_0000;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const infinity_bits_f: u32 = 0x7f80_0000;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const canonical_nan_bits_f: u32 = 0x7fc0_0000;
const minus_infinity_bits_f: u32 = 0xff80_0000;

fn log1phSmall(value: f32) f32 {
    @setFloatMode(.strict);
    const coefficient3 = fma_impl.fma(value, -0x1p-2, 0x1.555556p-2);
    const coefficient2 = fma_impl.fma(value, coefficient3, -0.5);
    return fma_impl.fma(value, value * coefficient2, value);
}

pub fn log1ph(value: f16) f16 {
    @setFloatMode(.strict);

    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;

    if (magnitude > infinity_bits_h) {
        return @bitCast(bits | quiet_nan_bit_h);
    }
    if (bits == 0xbc00) return @bitCast(minus_infinity_bits_h);
    if ((bits & sign_mask_h) != 0 and magnitude > 0x3c00) {
        return @bitCast(canonical_nan_bits_h);
    }
    if (magnitude == infinity_bits_h) return value;

    const working: f32 = @floatCast(value);
    const input_bits: u32 = @bitCast(working);
    const absolute_bits = input_bits & magnitude_mask_f;

    if (absolute_bits < 0x3d00_0000) {
        const result = log1phSmall(working);

        if (input_bits == 0xbbbf_a000) return @floatCast(-0x1.805ffep-8);
        if (input_bits == 0xbcfb_6000) return @floatCast(-0x1.fea004p-6);
        if (input_bits == 0x3bc0_6000) return @floatCast(0x1.7fa002p-8);

        if (absolute_bits <= 0x3400_0000) {
            if (absolute_bits == 0) return value;
            return @floatCast(working - 0x1p-45);
        }
        return @floatCast(result);
    }

    if (input_bits == 0x3dec_8000) return @floatCast(0x1.bfa006p-4);

    const adjusted: f32 = working + 1.0;
    const adjusted_bits: u32 = @bitCast(adjusted);
    const table_index: usize = @intCast((adjusted_bits >> 20) - 928);
    const reduction = data.half_reduction[table_index];
    const reduced = fma_impl.fma(reduction[0], adjusted, -1.0);
    return @floatCast(reduction[1] + log1phSmall(reduced));
}

fn log1pfSeries(value: f64, value_squared: f64) f64 {
    @setFloatMode(.strict);
    const coefficients = data.float_series;
    const value_fourth = value_squared * value_squared;
    const term0 = coefficients[1] + value * coefficients[2];
    const term1 = coefficients[3] + value * coefficients[4];
    const term2 = coefficients[6] + value * coefficients[7];
    const term3 = coefficients[5] + value * term2;
    return value_squared * ((term0 + value_squared * term1) + value_fourth * term3);
}

pub fn log1pf(value: f32) f32 {
    @setFloatMode(.strict);

    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask_f;

    if (magnitude > infinity_bits_f) {
        return @bitCast(bits | quiet_nan_bit_f);
    }
    if (bits == 0xbf80_0000) return @bitCast(minus_infinity_bits_f);
    if ((bits & sign_mask_f) != 0 and magnitude > 0x3f80_0000) {
        return @bitCast(canonical_nan_bits_f);
    }
    if (magnitude == infinity_bits_f) return value;

    const widened: f64 = @floatCast(value);
    if (magnitude < 0x3c88_0000) {
        if (magnitude < 0x3300_0000) {
            if (magnitude == 0) return value;
            return value;
        }

        const value_squared = widened * widened;
        const correction = log1pfSeries(widened, value_squared);
        var result = widened + correction;
        const result_bits: u64 = @bitCast(result);
        if ((result_bits & 0x0fff_ffff) == 0) {
            result += 0x1p+14 * (correction + (widened - result));
        }
        return @floatCast(result);
    }

    const one_plus_bits: u64 = @bitCast(widened + 1.0);
    const exponent: i32 = @as(i32, @intCast(one_plus_bits >> 52)) - 0x3ff;
    const mantissa = one_plus_bits & 0x000f_ffff_ffff_ffff;
    const table_index: usize = @intCast((one_plus_bits >> 47) & 31);
    const normalized: f64 = @bitCast(mantissa | (@as(u64, 0x3ff) << 52));
    const reduced = normalized * data.float_reciprocals[table_index] - 1.0;

    const exponent_float: f64 = @floatFromInt(exponent);
    const log_two = 0x1.62e42fefa39efp-1;
    const reduced_squared = reduced * reduced;
    const reduced_series = data.float_reduced_series;
    const polynomial0 = reduced_series[0] + reduced * reduced_series[1];
    const polynomial1 = reduced_series[2] + reduced * reduced_series[3];
    const result = (log_two * exponent_float + data.float_log_reciprocals[table_index]) +
        reduced * (polynomial0 + reduced_squared * polynomial1);

    const epsilon = 2.1555e-11;
    var upper_bound: f32 = @floatCast(result);
    const lower_bound: f32 = @floatCast(result - epsilon);
    if (upper_bound != lower_bound) {
        const correction = log1pfSeries(reduced, reduced_squared);
        const log_reciprocal = data.float_log_reciprocals[table_index] - 0x1.7654p-37;
        const log_two_low = 0x1.7f7d1cf79abcap-20;
        const log_two_high = 0x1.62e4p-1;
        const high_product = log_two_high * exponent_float;
        var low_product = log_two_low * exponent_float;
        low_product += reduced;
        const high = high_product + log_reciprocal;
        const low = ((high_product - high) + log_reciprocal) + (low_product + correction);
        const high_float: f32 = @floatCast(high + low);
        const low_residual = (high - @as(f64, high_float)) + low;
        var low_float: f32 = @floatCast(low_residual);
        const twice_low = low_float * 2.0;
        if ((high_float + twice_low) - high_float == twice_low) {
            const residual: f32 = @floatCast(low_residual - @as(f64, low_float));
            const residual_bits: u32 = @bitCast(residual);
            const signed_half: f32 = @bitCast((residual_bits & sign_mask_f) | 0x3f00_0000);
            const low_bits: u32 = @bitCast(low_float);
            const absolute_low: f32 = @bitCast(low_bits & magnitude_mask_f);
            low_float += signed_half * absolute_low;
        }
        upper_bound = high_float + low_float;
    }
    return upper_bound;
}
