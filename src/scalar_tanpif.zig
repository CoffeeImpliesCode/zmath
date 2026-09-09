// Correctly rounded half-revolution tangent for binary32 values.
//
// Copyright (c) 2022-2025 Alexei Sibidanov.
//
// This code is derived from CORE-MATH src/binary32/tanpi/tanpif.c at revision
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

const magnitude_mask: u32 = 0x7fff_ffff;
const infinity_bits: u32 = 0x7f80_0000;
const quiet_nan_bit: u32 = 0x0040_0000;
const canonical_nan_bits: u32 = 0x7fc0_0000;

inline fn roundEvenSmall(value: f64) f64 {
    // The callers bound |value| below 2^29. Adding 2^52 rounds to an
    // integer with ties-to-even, and subtracting the exact bias restores it.
    const rounded = value + 0x1.8p52;
    return rounded - 0x1.8p52;
}

fn copySign(value: f32, sign_source: f32) f32 {
    const value_bits: u32 = @bitCast(value);
    const sign_bits: u32 = @bitCast(sign_source);
    return @bitCast((value_bits & magnitude_mask) | (sign_bits & 0x8000_0000));
}

/// Return the correctly rounded tangent of pi times `value`.
///
/// The caller must keep IEEE round-to-nearest, ties-to-even active.
pub fn tanpif(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask;
    const exponent = magnitude & infinity_bits;
    if (magnitude > infinity_bits) return @bitCast(bits | quiet_nan_bit);
    if (magnitude == infinity_bits) return @bitCast(canonical_nan_bits);

    // The cutoff bounds |4*value| below 2^26, so both reductions fit the
    // bias method's exact |value| < 2^29 range.
    if (exponent > (150 << 23)) return @bitCast(bits & 0x8000_0000);
    const value_times_four: f32 = 4.0 * value;
    const rounded_times_four: f32 = @floatCast(roundEvenSmall(@floatCast(value_times_four)));
    const delta_times_four: f32 = value_times_four - rounded_times_four;
    const rounded_value: f32 = @floatCast(roundEvenSmall(@floatCast(value)));
    const reduced: f32 = value - rounded_value;

    if (delta_times_four == 0.0) {
        var quadrant: i64 = @intFromFloat(@min(@max(value_times_four, -0x1p+26), 0x1p+26));
        if ((quadrant & 1) != 0) return copySign(1.0, reduced);
        quadrant &= 6;
        if (quadrant == 0) return copySign(0.0, value);
        if (quadrant == 4) return -copySign(0.0, value);
        if (quadrant == 2) return @bitCast(@as(u32, 0x7f80_0000));
        return @bitCast(@as(u32, 0xff80_0000));
    }

    const reduced_double: f64 = @floatCast(reduced);
    const square = reduced_double * reduced_double;
    const square_four = square * square;
    const numerator_coefficients = [_]f64{
        0x1.921fb54442d19p-1,
        -0x1.1f458b3e1f8d6p-2,
        0x1.68a34bd0b8f6ap-6,
        -0x1.e4866f7a25f99p-13,
    };
    const denominator_coefficients = [_]f64{
        0x1p+0,
        -0x1.4b4b98d2df3a7p-1,
        0x1.8e9926d2bb901p-4,
        -0x1.a6f77fd847eep-9,
    };
    const numerator = (numerator_coefficients[0] + square * numerator_coefficients[1] +
        square_four * (numerator_coefficients[2] + square * numerator_coefficients[3])) *
        (reduced_double - reduced_double * square);
    const denominator = (denominator_coefficients[0] + square * denominator_coefficients[1] +
        square_four * (denominator_coefficients[2] + square * denominator_coefficients[3])) *
        (0.25 - square);
    const result: f64 = numerator / denominator;

    const reduced_bits: u32 = @bitCast(reduced);
    const reduced_magnitude = reduced_bits & magnitude_mask;
    if (reduced_magnitude == 0x3e93_3802) {
        return copySign(0x1.44cfbap+0, reduced) + copySign(0x1p-25, reduced);
    }
    if (reduced_magnitude == 0x38f2_6685) {
        return copySign(0x1.7cc304p-12, reduced) + copySign(0x1p-37, reduced);
    }
    return @floatCast(result);
}
