// Correctly rounded binary32 tangent for nearest-even rounding.
//
// Copyright (c) 2022-2025 Alexei Sibidanov.
//
// This code is derived from CORE-MATH src/binary32/tan/tanf.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT License.
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

const sign_mask: u32 = 0x8000_0000;
const magnitude_mask: u32 = 0x7fff_ffff;
const infinity_bits: u32 = 0x7f80_0000;
const quiet_nan_bit: u32 = 0x0040_0000;

fn roundEvenFinite(value: f64, integer: *i64) f64 {
    @setFloatMode(.strict);
    // reduceSmall gives |value| < 2^28. The bias has unit spacing at this range.
    const rounded = value + 0x1.8p52;
    const rounded_bits: u64 = @bitCast(rounded);
    const integer_low: i32 = @bitCast(@as(u32, @truncate(rounded_bits)));
    integer.* = @as(i64, integer_low);
    return rounded - 0x1.8p52;
}

// For |z| < 2^28, return r such that 2/pi*x = q + r.
fn reduceSmall(value: f32, quotient: *i64) f64 {
    @setFloatMode(.strict);
    const x: f64 = @floatCast(value);
    const low = -0x1.b1bbead603d8bp-32 * x;
    const high = 0x1.45f306ep-1 * x;
    const integral = roundEvenFinite(high, quotient);
    return (high - integral) + low;
}

// For |x| >= 2^28, return r such that 2/pi*x = q + r.
fn reduceLarge(bits: u32, quotient: *i64) f64 {
    const inverse_pi = [_]u64{
        0xfe51_63ab_debb_c562,
        0xdb62_9599_3c43_9041,
        0xfc27_57d1_f534_ddc0,
        0xa2f9_836e_4e44_1529,
    };
    const exponent: i64 = @intCast((bits >> 23) & 0xff);
    const mantissa: u64 = @as(u64, bits & 0x007f_ffff) | 0x0080_0000;
    const product_0: u128 = @as(u128, mantissa) * inverse_pi[0];
    const product_1: u128 = @as(u128, mantissa) * inverse_pi[1] + (product_0 >> 64);
    const product_2: u128 = @as(u128, mantissa) * inverse_pi[2] + (product_1 >> 64);
    const product_3: u128 = @as(u128, mantissa) * inverse_pi[3] + (product_2 >> 64);
    const product_3_high: u64 = @truncate(product_3 >> 64);
    const product_3_low: u64 = @truncate(product_3);
    const product_2_low: u64 = @truncate(product_2);
    const product_1_low: u64 = @truncate(product_1);
    const k = exponent - 127;
    const shift = k - 23;

    var integer_part: u64 = undefined;
    var fractional_part: u64 = undefined;
    if (shift < 64) {
        const amount: u6 = @intCast(shift);
        integer_part = (product_3_high << amount) | (product_3_low >> @as(u6, @intCast(64 - shift)));
        fractional_part = (product_3_low << amount) | (product_2_low >> @as(u6, @intCast(64 - shift)));
    } else if (shift == 64) {
        integer_part = product_3_low;
        fractional_part = product_2_low;
    } else {
        const amount: u6 = @intCast(shift - 64);
        integer_part = (product_3_low << amount) | (product_2_low >> @as(u6, @intCast(128 - shift)));
        fractional_part = (product_2_low << amount) | (product_1_low >> @as(u6, @intCast(128 - shift)));
    }

    const negative = (bits & sign_mask) != 0;
    const sign: i64 = if (negative) -1 else 0;
    const signed_fraction: i64 = @bitCast(fractional_part);
    const signed_adjustment: i64 = signed_fraction >> 63;
    var integer_value: i64 = @bitCast(integer_part);
    integer_value -= signed_adjustment;
    const reduced_fraction: i64 = signed_fraction ^ sign;
    integer_value = (integer_value ^ sign) - sign;
    quotient.* = integer_value;
    return @as(f64, @floatFromInt(reduced_fraction)) * 0x1p-64;
}

/// Return correctly rounded tangent with the caller's nearest-even mode.
///
/// The caller must set the floating-point environment to nearest-even before
/// calling this function.
pub fn tanf(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const exponent = (bits >> 23) & 0xff;
    if (exponent == 0xff) {
        if ((bits & 0x007f_ffff) != 0) return @bitCast(bits | quiet_nan_bit);
        return @bitCast(@as(u32, 0x7fc0_0000));
    }

    var quotient: i64 = undefined;
    var reduced: f64 = undefined;
    if (exponent < 127 + 28) {
        if (exponent < 115) {
            if (exponent < 102) return value;
            const square = value * value;
            return fma_impl.fma(value, 0x1.555556p-2 * square, value);
        }
        reduced = reduceSmall(value, &quotient);
    } else {
        reduced = reduceLarge(bits, &quotient);
    }

    const square = reduced * reduced;
    const square_4 = square * square;
    const numerator_coefficients = [_]f64{
        0x1.921fb54442d18p+0,
        -0x1.fd226e573289fp-2,
        0x1.b7a60c8dac9f6p-6,
        -0x1.725beb40f33e5p-13,
    };
    const denominator_coefficients = [_]f64{
        0x1p+0,
        -0x1.2395347fb829dp+0,
        0x1.2313660f29c36p-3,
        -0x1.9a707ab98d1c1p-9,
    };
    const numerator_1 = numerator_coefficients[0] + square * numerator_coefficients[1];
    const numerator_2 = numerator_coefficients[2] + square * numerator_coefficients[3];
    const numerator = (numerator_1 + square_4 * numerator_2) * reduced;
    const denominator_1 = denominator_coefficients[0] + square * denominator_coefficients[1];
    const denominator_2 = denominator_coefficients[2] + square * denominator_coefficients[3];
    const denominator = denominator_1 + square_4 * denominator_2;
    const sine_sign: f64 = if ((quotient & 1) != 0) 1.0 else 0.0;
    const cosine_sign: f64 = 1.0 - sine_sign;
    const result = (numerator * cosine_sign - denominator * sine_sign) /
        (numerator * sine_sign + denominator * cosine_sign);
    const result_bits: u64 = @bitCast(result);
    const tail = (result_bits +% 7) & 0x1fff_ffff;
    if (tail <= 14) {
        const Correction = struct {
            argument_bits: u32,
            high: f32,
            low: f32,
        };
        const corrections = [_]Correction{
            .{ .argument_bits = 0x3f8a_1f62, .high = 0x1.ddf9f6p+0, .low = -0x1.891d24p-52 },
            .{ .argument_bits = 0x4d56_d355, .high = 0x1.e80304p-3, .low = 0x1.419f46p-58 },
            .{ .argument_bits = 0x57d7_b0ed, .high = 0x1.60d1c8p-2, .low = -0x1.2d6c3ap-55 },
            .{ .argument_bits = 0x5980_445e, .high = 0x1.ca1edp+0, .low = 0x1.f6053p-53 },
            .{ .argument_bits = 0x63fc_86fe, .high = 0x1.597f9cp-1, .low = 0x1.925978p-53 },
            .{ .argument_bits = 0x6a66_2711, .high = -0x1.f33584p+1, .low = 0x1.d7254ap-51 },
            .{ .argument_bits = 0x6ad3_6709, .high = -0x1.c5612ep-1, .low = -0x1.26c33ep-53 },
            .{ .argument_bits = 0x72b5_05bb, .high = -0x1.e42a1ep+0, .low = -0x1.1dc906p-52 },
        };
        const magnitude = bits & magnitude_mask;
        const negative = (bits & sign_mask) != 0;
        for (corrections) |correction| {
            if (correction.argument_bits == magnitude) {
                return if (negative)
                    -correction.high - correction.low
                else
                    correction.high + correction.low;
            }
        }
    }
    return @floatCast(result);
}
