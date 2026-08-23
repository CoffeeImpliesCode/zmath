// Correctly rounded half-revolution sine for binary32 nearest-even rounding.
//
// Copyright (c) 2022-2025 Alexei Sibidanov.
//
// This code is derived from CORE-MATH src/binary32/sinpi/sinpif.c at revision
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

const data = @import("scalar_sinpif_data.zig");

const sign_mask: u32 = 0x8000_0000;
const magnitude_mask: u32 = 0x7fff_ffff;
const infinity_bits: u32 = 0x7f80_0000;
const quiet_nan_bit: u32 = 0x0040_0000;
const canonical_nan_bits: u32 = infinity_bits | quiet_nan_bit;

/// Return correctly rounded sinpi for a binary32 value.
///
/// The caller must keep IEEE round-to-nearest, ties-to-even active.
pub fn sinpif(value: f32) f32 {
    @setFloatMode(.strict);

    const bits: u32 = @bitCast(value);
    const exponent: i32 = @intCast((bits >> 23) & 0xff);
    if (exponent == 0xff) {
        if ((bits << 9) == 0) return @bitCast(canonical_nan_bits);
        return @bitCast(bits | quiet_nan_bit);
    }

    var mantissa: i32 = @intCast((bits & 0x007f_ffff) | 0x0080_0000);
    const sign: i32 = if ((bits & sign_mask) == 0) 0 else -1;
    mantissa = (mantissa ^ sign) - sign;
    const shift: i32 = 143 - exponent;

    if (shift < 0) {
        if (shift < -6) return @bitCast(bits & sign_mask);
        const index_bits: u32 = @bitCast(mantissa);
        const index: usize = @intCast((index_bits << @intCast(-shift - 1)) & 127);
        if (index == 0 or index == 64) return @bitCast(bits & sign_mask);
        return @floatCast(data.s[index]);
    }

    if (shift > 30) {
        const z: f64 = @floatCast(value);
        const z_squared = z * z;
        const result = z *
            (0x1.921fb54442d18p+1 + z_squared * -0x1.4abbce625be53p+2);
        return @floatCast(result);
    }

    const magnitude_bits: u32 = @bitCast(mantissa);
    const integer_shift = 25 - shift;
    if (integer_shift >= 0 and
        (magnitude_bits << @intCast(integer_shift)) == 0)
    {
        return @bitCast(bits & sign_mask);
    }

    const k_bits = magnitude_bits << @intCast(31 - shift);
    const k: i32 = @bitCast(k_bits);
    const z: f64 = @floatFromInt(k);
    const z_squared = z * z;
    const sine_correction = data.sn[0] + z_squared *
        (data.sn[1] + z_squared * data.sn[2]);
    const cosine_correction = data.cn[0] + z_squared *
        (data.cn[1] + z_squared * data.cn[2]);
    const quotient_bits: u32 = @bitCast(mantissa >> @intCast(shift));
    const quotient = (quotient_bits +% 1) >> 1;
    const sine_index: usize = @intCast(quotient & 127);
    const cosine_index: usize = @intCast((quotient + 32) & 127);
    const table_sine = data.s[sine_index];
    const table_cosine = data.s[cosine_index];
    const result = table_sine + (table_sine * z_squared) * cosine_correction +
        (table_cosine * z) * sine_correction;
    return @floatCast(result);
}
