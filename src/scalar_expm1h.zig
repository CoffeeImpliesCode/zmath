const data = @import("scalar_expm1h_data.zig");

const sign_mask: u16 = 0x8000;
const magnitude_mask: u16 = 0x7fff;
const infinity_bits: u16 = 0x7c00;
const quiet_nan_bit: u16 = 0x0200;

/// Correctly-rounded exponential function biased by one for binary16 values.
///
/// Copyright (c) 2025 Paul Zimmermann.
///
/// This file is part of the CORE-MATH project
/// (https://core-math.gitlabpages.inria.fr/).
///
/// Ported from src/binary16/expm1/expm1f16.c at revision
/// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in all
/// copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
/// SOFTWARE.
///
/// The caller must keep IEEE round-to-nearest, ties-to-even active.
pub fn expm1h(value: f16) f16 {
    @setFloatMode(.strict);

    const input_bits: u16 = @bitCast(value);
    const absolute_bits = input_bits & magnitude_mask;
    if (absolute_bits > infinity_bits) return @bitCast(input_bits | quiet_nan_bit);
    if (absolute_bits == infinity_bits) {
        if ((input_bits & sign_mask) == 0) return value;
    }
    if (absolute_bits == 0) return value;

    if (absolute_bits < 0x26d0) {
        const value_f32: f32 = value;
        if (absolute_bits <= 0x11a8) {
            if (input_bits == 0x11a8) return @floatCast(0x1.6a2002p-11);
            if (absolute_bits == 0x1) {
                if ((input_bits & sign_mask) == 0) return @floatCast(0x1.000002p-24);
                return @floatCast(-0x1.fffffep-25);
            }
            if (absolute_bits == 0x2) {
                if ((input_bits & sign_mask) == 0) return @floatCast(0x1.000002p-23);
                return @floatCast(-0x1.fffffep-24);
            }
        }
        const coefficient: f32 = @mulAdd(f32, value_f32, 0x1.555556p-3, 0.5);
        return @floatCast(@mulAdd(f32, coefficient * value_f32, value_f32, value_f32));
    }

    const table_1_index: usize = @intCast(input_bits >> 5);
    const table_2_index: usize = @intCast(((input_bits >> 10) << 5) | (input_bits & 0x1f));
    const first: f32 = @bitCast(data.t1[table_1_index]);
    const second: f32 = @bitCast(data.t2[table_2_index]);
    return @floatCast(@mulAdd(f32, first, second, -1.0));
}
