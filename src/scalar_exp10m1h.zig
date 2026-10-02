const fma_impl = @import("scalar_fma.zig");

const data = @import("scalar_exp10m1h_data.zig");

const sign_mask: u16 = 0x8000;
const magnitude_mask: u16 = 0x7fff;
const infinity_bits: u16 = 0x7c00;
const quiet_nan_bit: u16 = 0x0200;

/// Correctly-rounded base-10 exponential biased by one for binary16 values.
///
/// Copyright (c) 2025 Paul Zimmermann.
///
/// This file is part of the CORE-MATH project
/// (https://core-math.gitlabpages.inria.fr/).
///
/// Ported from src/binary16/exp10m1/exp10m1f16.c at revision
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
pub fn exp10m1h(value: f16) f16 {
    @setFloatMode(.strict);

    const input_bits: u16 = @bitCast(value);
    const absolute_bits = input_bits & magnitude_mask;
    if (absolute_bits > infinity_bits) return @bitCast(input_bits | quiet_nan_bit);
    if (absolute_bits == infinity_bits) {
        if ((input_bits & sign_mask) == 0) return value;
        return -1.0;
    }
    if (absolute_bits == 0) return value;

    if (absolute_bits < 0x2400) {
        switch (input_bits) {
            0x8757 => return @floatCast(-0x1.0e5ffep-12),
            0x9788 => return @floatCast(-0x1.14dffep-8),
            0xa1fe => return @floatCast(-0x1.b39feap-6),
            0x2147 => return @floatCast(0x1.89800ap-6),
            0x2394 => return @floatCast(0x1.1c001ap-5),
            0x23fc => return @floatCast(0x1.2b801p-5),
            else => {},
        }

        const value_f32: f32 = value;
        const scaled = 0x1.26bb1cp+1 * value_f32;
        const quadratic = fma_impl.fma(scaled, 0x1.555556p-2 * scaled, scaled);
        return @floatCast(fma_impl.fma(quadratic, 0.5 * scaled, scaled));
    }

    switch (input_bits) {
        0xa4a8 => return @floatCast(-0x1.50000ap-5),
        0xa4f6 => return @floatCast(-0x1.658004p-5),
        0x286c => return @floatCast(0x1.531ff4p-4),
        0x24e6 => return @floatCast(0x1.70ffd6p-5),
        0x2ce8 => return @floatCast(0x1.8b5ffap-3),
        0x2d71 => return @floatCast(0x1.badffcp-3),
        0x2615 => return @floatCast(0x1.cc9ffap-5),
        0x2633 => return @floatCast(0x1.d5bffcp-5),
        0x2ffa => return @floatCast(0x1.545ffep-2),
        else => {},
    }

    const table_1_index: usize = @intCast(input_bits >> 5);
    const table_2_index: usize = @intCast(((input_bits >> 10) << 5) | (input_bits & 0x1f));
    const first: f32 = @bitCast(data.t1[table_1_index]);
    const second: f32 = @bitCast(data.t2[table_2_index]);
    return @floatCast(fma_impl.fma(first, second, -1.0));
}
