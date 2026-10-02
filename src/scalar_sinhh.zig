const fma_impl = @import("scalar_fma.zig");

const data = @import("scalar_sinhh_data.zig");

const sign_mask: u16 = 0x8000;
const magnitude_mask: u16 = 0x7fff;
const infinity_bits: u16 = 0x7c00;
const quiet_nan_bit: u16 = 0x0200;

/// Correctly-rounded hyperbolic sine function for binary16 values.
///
/// Copyright (c) 2025 Paul Zimmermann.
///
/// This file is part of the CORE-MATH project
/// (https://core-math.gitlabpages.inria.fr/).
///
/// Ported from src/binary16/sinh/sinhf16.c at revision
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
pub fn sinhh(value: f16) f16 {
    @setFloatMode(.strict);
    const input_bits: u16 = @bitCast(value);
    const magnitude = input_bits & magnitude_mask;
    if (magnitude > infinity_bits) return @bitCast(input_bits | quiet_nan_bit);
    if (magnitude == infinity_bits) return value;

    const table_1_index: usize = @intCast(input_bits >> 5);
    const table_2_index: usize = @intCast(((input_bits >> 10) << 5) | (input_bits & 0x001f));
    const sinh_1: f32 = @bitCast(data.s1[table_1_index]);
    const cosh_1: f32 = @bitCast(data.c1[table_1_index]);
    const sinh_2: f32 = @bitCast(data.s2[table_2_index]);
    const cosh_2: f32 = @bitCast(data.c2[table_2_index]);
    return @floatCast(fma_impl.fma(sinh_1, cosh_2, cosh_1 * sinh_2));
}
