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
const data = @import("scalar_log2_data.zig");

// Correctly rounded binary16 log2, ported from CORE-MATH
// src/binary16/log2/log2f16.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
// Copyright (c) 2025 Maxence Ponsardin and Paul Zimmermann.
//
// This code is derived from CORE-MATH and used under the MIT License.

const sign_mask: u16 = 0x8000;
const magnitude_mask: u16 = 0x7fff;
const infinity_bits: u16 = 0x7c00;
const quiet_nan_bit: u16 = 0x0200;
const minus_infinity_bits: u16 = 0xfc00;
const domain_nan_bits: u16 = 0x7e00;

pub fn log2h(value: f16) f16 {
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask;
    const table_1_index: usize = @intCast(bits >> 10);
    const table_2_index: usize = @intCast(bits & 0x03ff);

    if ((bits & sign_mask) != 0) {
        if (bits == sign_mask) return @bitCast(minus_infinity_bits);
        if (magnitude > infinity_bits) return @bitCast(bits | quiet_nan_bit);
        return @bitCast(domain_nan_bits);
    }
    if (magnitude > infinity_bits) return @bitCast(bits | quiet_nan_bit);
    if (magnitude == infinity_bits) return value;
    if (table_1_index == 0) return @floatCast(@as(f32, @bitCast(data.t3[table_2_index])));
    const first: f32 = @bitCast(data.t1[table_1_index]);
    const second: f32 = @bitCast(data.t2[table_2_index]);
    return @floatCast(first + second);
}
