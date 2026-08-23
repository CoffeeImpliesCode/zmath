// Correctly rounded half-revolution tangent for binary16 values.
//
// Copyright (c) 2025 Paul Zimmermann.
//
// This code is derived from CORE-MATH src/binary16/tanpi/tanpif16.c at revision
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

const data = @import("scalar_tanpih_data.zig");

const magnitude_mask: u16 = 0x7fff;
const infinity_bits: u16 = 0x7c00;
const quiet_nan_bit: u16 = 0x0200;
const canonical_nan_bits: u16 = 0x7e00;

/// Return the correctly rounded tangent of pi times `value`.
///
/// The caller must keep IEEE round-to-nearest, ties-to-even active.
pub fn tanpih(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask;

    if (magnitude >= 0x3800) {
        if (magnitude >= infinity_bits) {
            if (magnitude > infinity_bits) return @bitCast(bits | quiet_nan_bit);
            return @bitCast(canonical_nan_bits);
        }

        const masks = [_]u16{
            0x0fff,
            0x03ff,
            0x01ff,
            0x00ff,
            0x007f,
            0x003f,
            0x001f,
            0x000f,
            0x0007,
            0x0003,
            0x0001,
            0x0000,
            0x0000,
            0x0000,
            0x0000,
            0x0000,
            0x0000,
        };
        const exponent: u16 = magnitude >> 10;
        if ((magnitude & masks[exponent - 14]) == 0) {
            const even_magnitude = exponent >= 26 or
                (magnitude & (@as(u16, 1) << @intCast(25 - exponent))) == 0;
            const positive = magnitude == bits;
            return @bitCast(if (even_magnitude == positive) @as(u16, 0x0000) else 0x8000);
        }
        if ((magnitude & masks[exponent - 13]) == 0) {
            const odd_half_magnitude = exponent == 14 or
                (magnitude & (@as(u16, 1) << @intCast(25 - exponent))) == 0;
            const positive = magnitude == bits;
            return @bitCast(if (odd_half_magnitude == positive) @as(u16, 0x7c00) else 0xfc00);
        }
    }

    if (bits == 0x8000) return value;

    const table_1_index: usize = @intCast(bits >> 5);
    const table_2_index: usize = @intCast(((bits >> 10) << 5) | (bits & 0x001f));
    const sine_1: f32 = @bitCast(data.s1[table_1_index]);
    const cosine_1: f32 = @bitCast(data.c1[table_1_index]);
    const sine_2: f32 = @bitCast(data.s2[table_2_index]);
    const cosine_2: f32 = @bitCast(data.c2[table_2_index]);
    const sine = @mulAdd(f32, sine_1, cosine_2, cosine_1 * sine_2);
    const cosine = @mulAdd(f32, cosine_1, cosine_2, -(sine_1 * sine_2));
    return @floatCast(sine / cosine);
}
