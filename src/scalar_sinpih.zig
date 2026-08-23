// Correctly rounded half-revolution sine for binary16 nearest-even rounding.
//
// Copyright (c) 2025 Paul Zimmermann.
//
// This code is derived from CORE-MATH src/binary16/sinpi/sinpif16.c at revision
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

const data = @import("scalar_sinpih_data.zig");

const sign_mask: u16 = 0x8000;
const magnitude_mask: u16 = 0x7fff;
const infinity_bits: u16 = 0x7c00;

/// Return correctly rounded sinpi for a binary16 value.
///
/// The caller must keep IEEE round-to-nearest, ties-to-even active.
pub fn sinpih(value: f16) f16 {
    @setFloatMode(.strict);

    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask;

    if (magnitude >= infinity_bits) {
        if (magnitude > infinity_bits) return @bitCast(bits | 0x0200);
        return @bitCast(@as(u16, 0x7e00));
    }

    if (magnitude >= 0x3c00) {
        const exponent: usize = magnitude >> 10;
        const integer_mask = [_]u16{
            0x3ff, 0x1ff, 0xff, 0x7f, 0x3f, 0x1f, 0xf, 0x7, 0x3, 0x1, 0x0,
            0x0,   0x0,   0x0,  0x0,  0x0,
        };
        if ((magnitude & integer_mask[exponent - 15]) == 0) {
            return @bitCast(if (magnitude == bits) @as(u16, 0) else @as(u16, sign_mask));
        }
    }

    if (magnitude == 0) return value;

    // This input is the one exceptional case retained by the upstream table.
    if (magnitude == 0x3bf7) {
        const result: f32 = if ((bits & sign_mask) == 0)
            0x1.c45ffep-7
        else
            -0x1.c45ffep-7;
        return @floatCast(result);
    }

    const table1_index: usize = bits >> 5;
    const table2_index: usize = ((bits >> 10) << 5) | (bits & 0x1f);
    const result = @mulAdd(
        f32,
        data.s1[table1_index],
        data.c2[table2_index],
        data.c1[table1_index] * data.s2[table2_index],
    );
    return @floatCast(result);
}
