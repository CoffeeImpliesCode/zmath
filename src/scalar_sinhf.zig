const sign_mask: u32 = 0x8000_0000;
const magnitude_mask: u32 = 0x7fff_ffff;
const infinity_bits: u32 = 0x7f80_0000;
const quiet_nan_bit: u32 = 0x0040_0000;

fn roundEvenFinite(value: f64) i64 {
    const truncated: i64 = @intFromFloat(value);
    const truncated_float: f64 = @floatFromInt(truncated);
    const fraction = value - truncated_float;
    const odd = (truncated & 1) != 0;
    if (fraction > 0.5 or (fraction == 0.5 and odd)) return truncated + 1;
    if (fraction < -0.5 or (fraction == -0.5 and odd)) return truncated - 1;
    return truncated;
}

/// Correctly-rounded hyperbolic sine function for binary32 values.
///
/// Copyright (c) 2022-2025 Alexei Sibidanov.
///
/// This file is part of the CORE-MATH project
/// (https://core-math.gitlabpages.inria.fr/).
///
/// Ported from src/binary32/sinh/sinhf.c at revision
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
pub fn sinhf(value: f32) f32 {
    @setFloatMode(.strict);

    const c = [4]f64{
        1.0,
        0x1.62e42fef4c4e7p-6,
        0x1.ebfd1b232f475p-13,
        0x1.c6b19384ecd93p-20,
    };
    const accurate_c = [7]f64{
        1.0,
        0x1.62e42fefa39efp-6,
        0x1.ebfbdff82c58fp-13,
        0x1.c6b08d702e0edp-20,
        0x1.3b2ab6fb92e5ep-27,
        0x1.5d886e6d54203p-35,
        0x1.430976b8ce6efp-43,
    };
    const table = [32]u64{
        0x3fe0000000000000, 0x3fe059b0d3158574, 0x3fe0b5586cf9890f,
        0x3fe11301d0125b51, 0x3fe172b83c7d517b, 0x3fe1d4873168b9aa,
        0x3fe2387a6e756238, 0x3fe29e9df51fdee1, 0x3fe306fe0a31b715,
        0x3fe371a7373aa9cb, 0x3fe3dea64c123422, 0x3fe44e086061892d,
        0x3fe4bfdad5362a27, 0x3fe5342b569d4f82, 0x3fe5ab07dd485429,
        0x3fe6247eb03a5585, 0x3fe6a09e667f3bcd, 0x3fe71f75e8ec5f74,
        0x3fe7a11473eb0187, 0x3fe82589994cce13, 0x3fe8ace5422aa0db,
        0x3fe93737b0cdc5e5, 0x3fe9c49182a3f090, 0x3fea5503b23e255d,
        0x3feae89f995ad3ad, 0x3feb7f76f2fb5e47, 0x3fec199bdd85529c,
        0x3fecb720dcef9069, 0x3fed5818dcfba487, 0x3fedfc97337b9b5f,
        0x3feea4afa2a490da, 0x3fef50765b6e4540,
    };

    const input_bits: u32 = @bitCast(value);
    const shifted_bits = input_bits << 1;
    if (shifted_bits > 0x8565_a9f8) {
        if (shifted_bits >= 0xff00_0000) {
            if ((shifted_bits << 8) != 0) return @bitCast(input_bits | quiet_nan_bit);
            return value;
        }
        return @bitCast((input_bits & sign_mask) | infinity_bits);
    }
    if (shifted_bits < 0x7c00_0000) {
        if (shifted_bits <= 0x7425_0bfe) {
            if (shifted_bits < 0x6600_0000) {
                if ((input_bits & magnitude_mask) == 0) return value;
                const absolute_value: f32 = @bitCast(input_bits & magnitude_mask);
                return @mulAdd(f32, value, absolute_value, value);
            }
            if (shifted_bits == 0x7425_0bfe) {
                const sign: f32 = if ((input_bits & sign_mask) == 0) 1.0 else -1.0;
                return sign * 0x1.250bfep-11 + sign * 0x1p-36;
            }
            return ((value * 0x1.555556p-3) * (value * value)) + value;
        }
        const z: f64 = value;
        const z_squared = z * z;
        const z_fourth = z_squared * z_squared;
        const cp = [4]f64{
            0x1.5555555555555p-3,
            0x1.11111111146e1p-7,
            0x1.a01a00930dda6p-13,
            0x1.71f92198aa6e9p-19,
        };
        return @floatCast(z + (z_squared * z) *
            ((cp[0] + z_squared * cp[1]) +
                z_fourth * (cp[2] + z_squared * cp[3])));
    }

    const z: f64 = value;
    const inverse_ln2: f64 = 0x1.71547652b82fep+5;
    const scaled = inverse_ln2 * z;
    const rounded_integer = roundEvenFinite(scaled);
    const rounded: f64 = @floatFromInt(rounded_integer);
    const h = scaled - rounded;
    const h_squared = h * h;
    const rounding_constant: f64 = 0x1.8p+52;
    const ja_bits: u64 = @bitCast(rounded + rounding_constant);
    const positive_index: i64 = @bitCast(ja_bits);
    const negative_index = -positive_index;
    const positive_index_bits: u64 = @bitCast(positive_index);
    const negative_index_bits: u64 = @bitCast(negative_index);
    const positive_bits = table[@intCast(positive_index_bits & 0x1f)] +%
        ((@as(u64, @bitCast(positive_index >> 5))) << 52);
    const negative_bits = table[@intCast(negative_index_bits & 0x1f)] +%
        ((@as(u64, @bitCast(negative_index >> 5))) << 52);
    const positive_scale: f64 = @bitCast(positive_bits);
    const negative_scale: f64 = @bitCast(negative_bits);
    const even = c[0] + h_squared * c[2];
    const odd = c[1] + h_squared * c[3];
    const positive = positive_scale * (even + h * odd);
    const negative = negative_scale * (even - h * odd);
    const result = positive - negative;
    var upper: f32 = @floatCast(result);
    const lower: f32 = @floatCast(result - 1.52e-10 * result);
    if (upper != lower) {
        const inverse_ln2_high: f64 = 0x1.7154765p+5;
        const inverse_ln2_low: f64 = 0x1.5c17f0bbbe88p-26;
        const corrected_h = (inverse_ln2_high * z - rounded) + inverse_ln2_low * z;
        const corrected_squared = corrected_h * corrected_h;
        const accurate_even = accurate_c[0] + corrected_squared * accurate_c[2] +
            (corrected_squared * corrected_squared) *
                (accurate_c[4] + corrected_squared * accurate_c[6]);
        const accurate_odd = accurate_c[1] + corrected_squared *
            (accurate_c[3] + corrected_squared * accurate_c[5]);
        const accurate_result = positive_scale *
            (accurate_even + corrected_h * accurate_odd) -
            negative_scale * (accurate_even - corrected_h * accurate_odd);
        upper = @floatCast(accurate_result);
    }
    return upper;
}
