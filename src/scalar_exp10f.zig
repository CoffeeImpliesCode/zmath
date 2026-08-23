// Correctly rounded 10^x function for binary32 values.
//
// Copyright (c) 2023-2025 Alexei Sibidanov.
//
// This file is part of the CORE-MATH project
// (https://core-math.gitlabpages.inria.fr/).
//
// Ported from src/binary32/exp10/exp10f.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
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

const data = @import("scalar_exp10f_data.zig");

const sign_mask: u32 = 0x8000_0000;
const magnitude_mask: u32 = 0x7fff_ffff;
const infinity_bits: u32 = 0x7f80_0000;
const quiet_nan_bit: u32 = 0x0040_0000;

const exact_powers = [_]f32{
    10.0,
    100.0,
    1000.0,
    10000.0,
    100000.0,
    1000000.0,
    10000000.0,
    100000000.0,
    1000000000.0,
    10000000000.0,
};

fn roundEvenFinite(value: f64) f64 {
    @setFloatMode(.strict);

    const truncated: i64 = @intFromFloat(value);
    const truncated_value: f64 = @floatFromInt(truncated);
    const fraction = value - truncated_value;
    const odd = (truncated & 1) != 0;
    if (fraction > 0.5 or (fraction == 0.5 and odd)) {
        return @floatFromInt(truncated + 1);
    }
    if (fraction < -0.5 or (fraction == -0.5 and odd)) {
        return @floatFromInt(truncated - 1);
    }
    return truncated_value;
}

fn exp10fSpecial(value: f32, bits: u32) f32 {
    @setFloatMode(.strict);

    const unsigned_magnitude = bits << 1;
    if (unsigned_magnitude > 0xff00_0000) return @bitCast(bits | quiet_nan_bit);
    if (unsigned_magnitude == 0xff00_0000) {
        if ((bits & sign_mask) != 0) return 0.0;
        return @bitCast(infinity_bits);
    }
    if ((bits & sign_mask) != 0) {
        const value_wide: f64 = value;
        var result = 0x1p-149 + (value_wide + 0x1.66d3e7bd9a403p+5) * 0x1.a934f0979a37p-149;
        if (result < 0x1p-151) result = 0x1p-151;
        return @floatCast(result);
    }
    return @bitCast(infinity_bits);
}

fn exp10fFinite(value: f32, bits: u32) f32 {
    @setFloatMode(.strict);

    const unsigned_magnitude = bits << 1;
    if (unsigned_magnitude > 0x8434_4134 or unsigned_magnitude < 0x72ad_f1c6) {
        if (unsigned_magnitude < 0x72ad_f1c6) {
            const value_wide: f64 = value;
            return @floatCast(1.0 + value_wide *
                (0x1.26bb1bbb55516p+1 + value_wide *
                    (0x1.53524c73cea69p+1 + value_wide * 0x1.0470591de2ca4p+1)));
        }
        if (unsigned_magnitude >= 0xff00_0000) return exp10fSpecial(value, bits);
        if (bits > 0xc233_69f4) {
            const value_wide: f64 = value;
            var result = 0x1p-149 +
                (value_wide + 0x1.66d3e7bd9a403p+5) * 0x1.a934f0979a37p-149;
            if (result < 0x1p-151) result = 0x1p-151;
            return @floatCast(result);
        }
        if ((bits & sign_mask) == 0) return @bitCast(infinity_bits);
    }

    if ((bits << 12) == 0) {
        const exponent = bits >> 20;
        if (exponent >= 1016) {
            const shift_value = exponent - 1016;
            if (shift_value <= 26) {
                const shift: u5 = @intCast(shift_value);
                const table_bit: u32 = @as(u32, 1) << shift;
                const mask: u32 = 0x0755_1101;
                if ((table_bit & mask) != 0) {
                    const index: usize = @intCast(@popCount(mask & (table_bit - 1)));
                    return exact_powers[index];
                }
            }
        }
    }

    const value_wide: f64 = value;
    const iln102: f64 = 0x1.a934f0979a371p+6;
    const iln102h: f64 = 0x1.a934f09p+1;
    const iln102l: f64 = 0x1.e68dc57f2496p-29;
    const scaled = iln102 * value_wide;
    const rounded = roundEvenFinite(scaled);
    const reduced = scaled - rounded;
    const rounded_integer: i64 = @intFromFloat(rounded);
    const table_index: usize = @intCast(rounded_integer & 0x1f);
    const table_scale_bits = data.table[table_index];
    const exponent_bits = @as(u64, @bitCast(rounded_integer >> 5)) << 52;
    const scale: f64 = @bitCast(table_scale_bits +% exponent_bits);
    const reduced_squared = reduced * reduced;
    const fast_coefficients = [_]f64{
        1.0,
        0x1.62e42fef4c4e7p-6,
        0x1.ebfd1b232f475p-13,
        0x1.c6b19384ecd93p-20,
    };
    var result = scale * ((fast_coefficients[0] + reduced * fast_coefficients[1]) +
        reduced_squared * (fast_coefficients[2] + reduced * fast_coefficients[3]));
    const error_bound: f64 = 1.45e-10;
    var upper: f32 = @floatCast(result);
    const lower: f32 = @floatCast(result - result * error_bound);

    if (upper != lower) {
        const accurate_reduced = (iln102h * value_wide - rounded * 0.03125) + iln102l * value_wide;
        const accurate_squared = accurate_reduced * accurate_reduced;
        const accurate_product = scale * accurate_reduced;
        result = scale + accurate_product *
            ((0x1.62e42fefa39efp-1 + accurate_reduced * 0x1.ebfbdff82c58fp-3) +
                accurate_squared *
                    ((0x1.c6b08d702e0edp-5 + accurate_reduced * 0x1.3b2ab6fb92e5ep-7) +
                        accurate_squared *
                            (0x1.5d886e6d54203p-10 +
                                accurate_reduced * 0x1.430976b8ce6efp-13)));
        upper = @floatCast(result);
    }
    return upper;
}

/// Computes 10^value in binary32 with caller-supplied nearest-even rounding.
pub fn exp10f(value: f32) f32 {
    @setFloatMode(.strict);

    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask;
    if (magnitude > infinity_bits) return @bitCast(bits | quiet_nan_bit);
    if (magnitude == infinity_bits) {
        if ((bits & sign_mask) != 0) return 0.0;
        return @bitCast(infinity_bits);
    }
    return exp10fFinite(value, bits);
}
