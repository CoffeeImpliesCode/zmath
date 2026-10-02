// Correctly rounded binary64 tangent for nearest-even rounding.
//
// Copyright (c) 2022-2025 Paul Zimmermann and Tom Hubrecht.
//
// This code is derived from CORE-MATH src/binary64/tan/tan.c at revision
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
//
// The finite Dint arithmetic is shared with the pinned CORE-MATH sine port.
const fma_impl = @import("scalar_fma.zig");

const sind_impl = @import("scalar_sind.zig");

const magnitude_mask: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits: u64 = infinity_bits | quiet_nan_bit;

pub fn tand(value: f64) f64 {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    const magnitude = bits & magnitude_mask;
    if (magnitude > infinity_bits) return @bitCast(bits | quiet_nan_bit);
    if (magnitude == infinity_bits) return @bitCast(canonical_nan_bits);
    if (magnitude <= 0x3e4d_12ed_0af1_a27e) return value;
    if (magnitude == 0x3e9d_ffff_ffff_ff1f) {
        return @bitCast(@as(u64, 0x3e9e_0000_0000_0151) | (bits & 0x8000_0000_0000_0000));
    }
    if (magnitude == 0x3ead_ffff_ffff_fc7c) {
        return @bitCast(@as(u64, 0x3eae_0000_0000_0546) | (bits & 0x8000_0000_0000_0000));
    }
    if (magnitude <= 0x3fc0_0000_0000_0000) return tanSmall(value);
    var high: f64 = undefined;
    var low: f64 = undefined;
    const error_bound = sind_impl.tanFast(&high, &low, value);
    const left = high + (low - error_bound);
    const right = high + (low + error_bound);
    if (left != right) {
        @branchHint(.cold);
        return sind_impl.tanAccurate(value);
    }
    return left;
}

inline fn tanSmall(value: f64) f64 {
    @setFloatMode(.strict);
    // On |x| <= 1/8, the positive Taylor tail after degree seventeen is
    // bounded by c19*|x|^19/(1-4*x^2/pi^2). Together with evaluation
    // error, this is below 2u*|x^3|, where u=2^-53.
    const square = value * value;
    const cubic = value * square;
    var polynomial: f64 = 0x1.3558248036744p-11;
    polynomial = fma_impl.fma64(polynomial, square, 0x1.7da36452b75e3p-10);
    polynomial = fma_impl.fma64(polynomial, square, 0x1.d6d3d0e157de0p-9);
    polynomial = fma_impl.fma64(polynomial, square, 0x1.226e355e6c23dp-7);
    polynomial = fma_impl.fma64(polynomial, square, 0x1.664f4882c10fap-6);
    polynomial = fma_impl.fma64(polynomial, square, 0x1.ba1ba1ba1ba1cp-5);
    polynomial = fma_impl.fma64(polynomial, square, 0x1.1111111111111p-3);
    polynomial = fma_impl.fma64(polynomial, square, 0x1.5555555555555p-2);
    const correction = cubic * polynomial;
    const high = value + correction;
    const low = correction - (high - value);
    const error_bound = @abs(cubic) * 0x1p-52;
    const left = high + (low - error_bound);
    const right = high + (low + error_bound);
    if (left == right) return left;
    return sind_impl.tanAccurate(value);
}
