const fma_impl = @import("scalar_fma.zig");

const rsqrth_data = @import("scalar_rsqrth_data.zig");

const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const quiet_nan_bit_h: u16 = 0x0200;
const canonical_nan_bits_h: u16 = 0x7e00;

const sign_mask_f: u32 = 0x8000_0000;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const infinity_bits_f: u32 = 0x7f80_0000;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const canonical_nan_bits_f: u32 = 0x7fc0_0000;

const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const fraction_mask_d: u64 = 0x000f_ffff_ffff_ffff;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits_d: u64 = 0x7ff8_0000_0000_0000;

// Correctly rounded reciprocal square root for binary16 values.
//
// Copyright (c) 2025 Maxence Ponsardin and Paul Zimmermann.
//
// This code is derived from CORE-MATH src/binary16/rsqrt/rsqrtf16.c at
// revision 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT
// License:
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
/// Computes the correctly rounded reciprocal square root of a binary16 value.
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn rsqrth(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;

    if (magnitude > infinity_bits_h) return @bitCast(bits | quiet_nan_bit_h);
    if ((bits & sign_mask_h) != 0) {
        if (magnitude == 0) return @bitCast(@as(u16, 0xfc00));
        return @bitCast(canonical_nan_bits_h);
    }
    if (magnitude == infinity_bits_h) return @bitCast(@as(u16, 0x0000));
    if (magnitude == 0) return @bitCast(@as(u16, 0x7c00));

    const exponent_index = bits >> 10;
    const fraction_index = bits & 0x03ff;
    if (exponent_index == 0) {
        const table_value: f32 = @bitCast(rsqrth_data.t3[fraction_index]);
        return @floatCast(table_value);
    }
    const table1_value: f32 = @bitCast(rsqrth_data.t1[exponent_index]);
    const table2_value: f32 = @bitCast(rsqrth_data.t2[fraction_index]);
    return @floatCast(table1_value * table2_value);
}

// Correctly rounded reciprocal square root for binary32 values.
//
// Copyright (c) 2022-2023 Alexei Sibidanov.
//
// This code is derived from CORE-MATH src/binary32/rsqrt/rsqrtf.c at
// revision 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT
// License:
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
noinline fn rsqrtfSpecial(bits: u32) f32 {
    const magnitude = bits & magnitude_mask_f;
    if (magnitude > infinity_bits_f) return @bitCast(bits | quiet_nan_bit_f);
    if ((bits & sign_mask_f) != 0) {
        if (magnitude == 0) return @bitCast(@as(u32, 0xff80_0000));
        return @bitCast(canonical_nan_bits_f);
    }
    if (magnitude == infinity_bits_f) return 0.0;
    return @bitCast(infinity_bits_f);
}

/// Computes the correctly rounded reciprocal square root of a binary32 value.
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn rsqrtf(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    if (bits -% 1 >= infinity_bits_f - 1) return rsqrtfSpecial(bits);

    const shifted_bits = bits << 8;
    if (bits == 0x002f_7e2a or shifted_bits == 0xbd_f8a800 or
        shifted_bits == 0x55_b7bd00)
    {
        if (bits != 0x0055_b7bd) {
            var exponent: u32 = bits >> 23;
            var table_index: usize = 1;
            if (bits == 0x002f_7e2a) exponent = 0xffff_ffff;
            if (shifted_bits == 0x55_b7bd00) table_index = 0;
            const table = [_]u32{ 0x000c_1740, 0x0052_22e0 };
            exponent = ((512 -% exponent) / 2) -% 578;
            const result_bits = table[table_index] | (exponent << 23);
            const correction_bits = (exponent -% 25) << 23;
            const result: f32 = @bitCast(result_bits);
            const correction: f32 = @bitCast(correction_bits);
            return result - correction;
        }
    }

    const widened: f64 = @floatCast(value);
    return @floatCast((1.0 / widened) * @sqrt(widened));
}

// Correctly rounded reciprocal square root for binary64 values.
//
// Copyright (c) 2022-2025 Alexei Sibidanov.
//
// This code is derived from CORE-MATH src/binary64/rsqrt/rsqrt.c at
// revision 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT
// License:
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
noinline fn rsqrtdRefine(initial: f64, input: f64) f64 {
    @setFloatMode(.strict);
    var result_bits: u64 = @bitCast(initial);
    const input_bits: u64 = @bitCast(input);

    if ((input_bits << 11) == (1 << 63)) return initial;

    const exponent_bit: u64 = (input_bits >> 52) & 1;
    const result_mantissa: u64 = ((result_bits << 11) | (1 << 63)) >> 11;
    const input_mantissa: u64 =
        ((input_bits & (~@as(u64, 0) >> 12)) | (1 << 52)) << @intCast(5 - exponent_bit);

    const product: u128 = @as(u128, result_mantissa) * input_mantissa;
    const product_high: u64 = @truncate(product >> 64);
    const product_low: u64 = @truncate(product);

    var remainder_product: u128 = @as(u128, product_low) * result_mantissa;
    const low_word: u64 = @truncate(remainder_product);
    const high_word: u64 = @truncate((remainder_product >> 64) +
        @as(u128, product_high) * result_mantissa);
    remainder_product = (@as(u128, high_word) << 64) | low_word;
    const sign_bit: u128 = remainder_product >> 127;
    const direction: u64 = if (sign_bit == 0) 1 else ~@as(u64, 0);
    const signed_remainder: u128 =
        ((product << 1) ^ (0 -% sign_bit)) +% sign_bit;
    var input_mantissa_negated: u64 = -%input_mantissa;
    const input_mantissa_double = input_mantissa << 1;
    var previous_remainder: u128 = 0;
    // The compensated pair has relative error below 2^-100 before the
    // final rounding. The candidate is at most two encodings from the root;
    // this recurrence moves one encoding toward the root per step.
    var crossed_root = false;
    for (0..3) |_| {
        result_bits -%= direction;
        previous_remainder = remainder_product;
        input_mantissa_negated +%= input_mantissa_double;
        const difference = signed_remainder -% input_mantissa_negated;
        remainder_product -%= difference;
        if (((previous_remainder ^ remainder_product) >> 127) != 0) {
            crossed_root = true;
            break;
        }
    }
    if (!crossed_root) unreachable;

    if ((remainder_product >> 127) == 0) result_bits +%= direction;
    if ((remainder_product >> 127) == 0) remainder_product = previous_remainder;

    const refined_mantissa: u64 = ((result_bits << 11) | (1 << 63)) >> 11;
    const refined_product: u128 = @as(u128, refined_mantissa) * input_mantissa;
    remainder_product +%= @as(u128, input_mantissa) >> 2;
    remainder_product +%= refined_product;
    result_bits +%= @truncate(remainder_product >> 127);
    return @bitCast(result_bits);
}

inline fn rsqrtdFinite(value: f64) f64 {
    @setFloatMode(.strict);
    const reciprocal = 1.0 / value;
    const estimate = @sqrt(value) * reciprocal;
    const negative_half = -0.5 * value;
    const sigma = fma_impl.fma64(negative_half, reciprocal, 0.5);
    const tau = fma_impl.fma64(estimate, estimate, -reciprocal);
    const nu = fma_impl.fma64(negative_half, tau, sigma);
    const refined = fma_impl.fma64(estimate, nu, estimate);
    const correction = fma_impl.fma64(estimate, nu, estimate - refined);

    // With u=2^-53, the seed error is below (3+4u)u and nu approximates
    // (1-value*estimate^2)/2 with absolute error below 8u^2. The Newton
    // remainder and all rounding terms total less than 32u^2=2^-101
    // relative. The fused residual contributes less than 2u^2.
    // Rounding cells are symmetric except below powers of two. Input
    // spacing excludes a root below a power when refined equals that
    // power, so one outward endpoint certifies both directions.
    const outside = refined + fma_impl.fma64(refined, 0x1p-100, @abs(correction));
    if (outside != refined) {
        @branchHint(.unlikely);
        return rsqrtdRefine(refined, value);
    }
    return refined;
}

noinline fn rsqrtdScaled(bits: u64) f64 {
    @setFloatMode(.strict);
    if (bits -% 1 >= infinity_bits_d - 1) {
        const magnitude = bits & magnitude_mask_d;
        if (magnitude > infinity_bits_d) return @bitCast(bits | quiet_nan_bit_d);
        if ((bits & sign_mask_d) != 0) {
            if (magnitude == 0) return @bitCast(@as(u64, 0xfff0_0000_0000_0000));
            return @bitCast(canonical_nan_bits_d);
        }
        if (magnitude == infinity_bits_d) return 0.0;
        return @bitCast(infinity_bits_d);
    }

    var normalized_bits: u64 = undefined;
    var scale_exponent: i16 = undefined;
    const exponent_bits = bits >> 52;
    if (exponent_bits == 0) {
        const leading_zeroes: i32 = @intCast(@clz(bits));
        const highest_bit = 63 - leading_zeroes;
        const parity = highest_bit & 1;
        const normalization_shift: u6 = @intCast(52 - highest_bit);
        normalized_bits = ((bits << normalization_shift) & fraction_mask_d) |
            (@as(u64, @intCast(0x3ff + parity)) << 52);
        scale_exponent = @intCast(@divExact(1074 + parity - highest_bit, 2));
    } else {
        const exponent: i32 = @intCast(exponent_bits);
        const parity = (exponent + 1) & 1;
        normalized_bits = (bits & fraction_mask_d) |
            (@as(u64, @intCast(0x3ff + parity)) << 52);
        scale_exponent = @intCast(@divExact(0x3ff + parity - exponent, 2));
    }

    // Scaling into [1,4) keeps every residual normal. Every finite input
    // has a normal reciprocal root, so restoring the exponent is exact.
    const normalized: f64 = @bitCast(normalized_bits);
    const result = rsqrtdFinite(normalized);
    const exponent_adjustment = @as(u64, @bitCast(@as(i64, scale_exponent))) << 52;
    return @bitCast(@as(u64, @bitCast(result)) +% exponent_adjustment);
}

/// Computes the correctly rounded reciprocal square root of a binary64 value.
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn rsqrtd(value: f64) f64 {
    const bits: u64 = @bitCast(value);
    // The reciprocal and value/2 stay normal. For value <= 2^970, a
    // subnormal square residual is an exact multiple of 2^-1074.
    if (bits -% 0x0020_0000_0000_0000 <= 0x7c90_0000_0000_0000 - 0x0020_0000_0000_0000)
        return rsqrtdFinite(value);
    return rsqrtdScaled(bits);
}
