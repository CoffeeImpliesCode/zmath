const sign_mask: u32 = 0x8000_0000;
const magnitude_mask: u32 = 0x7fff_ffff;
const infinity_bits: u32 = 0x7f80_0000;
const quiet_nan_bit: u32 = 0x0040_0000;

const table: [32]f64 = .{
    0x1p+0,
    0x1.059b0d3158574p+0,
    0x1.0b5586cf9890fp+0,
    0x1.11301d0125b51p+0,
    0x1.172b83c7d517bp+0,
    0x1.1d4873168b9aap+0,
    0x1.2387a6e756238p+0,
    0x1.29e9df51fdee1p+0,
    0x1.306fe0a31b715p+0,
    0x1.371a7373aa9cbp+0,
    0x1.3dea64c123422p+0,
    0x1.44e086061892dp+0,
    0x1.4bfdad5362a27p+0,
    0x1.5342b569d4f82p+0,
    0x1.5ab07dd485429p+0,
    0x1.6247eb03a5585p+0,
    0x1.6a09e667f3bcdp+0,
    0x1.71f75e8ec5f74p+0,
    0x1.7a11473eb0187p+0,
    0x1.82589994cce13p+0,
    0x1.8ace5422aa0dbp+0,
    0x1.93737b0cdc5e5p+0,
    0x1.9c49182a3f09p+0,
    0x1.a5503b23e255dp+0,
    0x1.ae89f995ad3adp+0,
    0x1.b7f76f2fb5e47p+0,
    0x1.c199bdd85529cp+0,
    0x1.cb720dcef9069p+0,
    0x1.d5818dcfba487p+0,
    0x1.dfc97337b9b5fp+0,
    0x1.ea4afa2a490dap+0,
    0x1.f50765b6e454p+0,
};

fn roundEvenFinite(value: f64) i64 {
    const truncated: i64 = @intFromFloat(value);
    const truncated_float: f64 = @floatFromInt(truncated);
    const fraction = value - truncated_float;
    const odd = (truncated & 1) != 0;
    if (fraction > 0.5 or (fraction == 0.5 and odd)) return truncated + 1;
    if (fraction < -0.5 or (fraction == -0.5 and odd)) return truncated - 1;
    return truncated;
}

/// Correctly-rounded exponential function biased by one for binary32 values.
///
/// Copyright (c) 2022-2025 Alexei Sibidanov.
///
/// This file is part of the CORE-MATH project
/// (https://core-math.gitlabpages.inria.fr/).
///
/// Ported from src/binary32/expm1/expm1f.c at revision
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
pub fn expm1f(value: f32) f32 {
    @setFloatMode(.strict);

    const input_bits: u32 = @bitCast(value);
    const absolute_bits = input_bits & magnitude_mask;
    if (absolute_bits > infinity_bits) return @bitCast(input_bits | quiet_nan_bit);
    if (absolute_bits == infinity_bits) {
        if ((input_bits & sign_mask) == 0) return @bitCast(infinity_bits);
        return @bitCast(@as(u32, 0xbf80_0000));
    }

    const value_wide: f64 = value;
    if (absolute_bits < 0x3e20_0000) {
        if (absolute_bits < 0x33b5_04f4) {
            if (absolute_bits == 0) return value;
            const absolute_value: f32 = @bitCast(absolute_bits);
            return @mulAdd(f32, absolute_value, 0x1p-25, value);
        }
        const z_squared = value_wide * value_wide;
        const z_fourth = z_squared * z_squared;
        const b = [8]f64{
            0x1.fffffffffffc2p-2,
            0x1.55555555555fep-3,
            0x1.555555559767fp-5,
            0x1.1111111098dc1p-7,
            0x1.6c16bca988aa9p-10,
            0x1.a01a07658483fp-13,
            0x1.a05b04d2c3503p-16,
            0x1.71de3a960b5e3p-19,
        };
        return @floatCast(value_wide + z_squared *
            ((b[0] + value_wide * b[1]) + z_squared *
                (b[2] + value_wide * b[3]) + z_fourth *
                ((b[4] + value_wide * b[5]) + z_squared *
                    (b[6] + value_wide * b[7]))));
    }

    if (absolute_bits >= 0x42b1_7218) {
        if ((input_bits & sign_mask) != 0) return @bitCast(@as(u32, 0xbf80_0000));
        return @bitCast(infinity_bits);
    }

    const c = [4]f64{
        1.0,
        0x1.62e42fef4c4e7p-6,
        0x1.ebfd1b232f475p-13,
        0x1.c6b19384ecd93p-20,
    };
    const accurate_c = [6]f64{
        0x1.62e42fefa39efp-6,
        0x1.ebfbdff82c58fp-13,
        0x1.c6b08d702e0edp-20,
        0x1.3b2ab6fb92e5ep-27,
        0x1.5d886e6d54203p-35,
        0x1.430976b8ce6efp-43,
    };
    const inverse_ln2: f64 = 0x1.71547652b82fep+5;
    const rounding_constant: f64 = 0x1.8p+52;
    const scaled = inverse_ln2 * value_wide;
    const rounded_integer = roundEvenFinite(scaled);
    const rounded: f64 = @floatFromInt(rounded_integer);
    const remainder = scaled - rounded;
    const remainder_squared = remainder * remainder;
    const scale_integer: u64 = @bitCast(rounded + rounding_constant);
    const table_bits: u64 = @bitCast(table[@intCast(scale_integer & 0x1f)]);
    const scale_bits = table_bits +% ((scale_integer >> 5) << 52);
    const scale: f64 = @bitCast(scale_bits);
    const c0 = c[0] + remainder * c[1];
    const c2 = c[2] + remainder * c[3];
    var result = (c0 + remainder_squared * c2) * scale - 1.0;
    var upper: f32 = @floatCast(result);
    const lower: f32 = @floatCast(result - scale * 0x1.3b3p-33);
    if (upper != lower) {
        if (input_bits > 0xc18a_a123) return @bitCast(@as(u32, 0xbf80_0000));
        const inverse_ln2_high: f64 = 0x1.7154765p+5;
        const inverse_ln2_low: f64 = 0x1.5c17f0bbbe88p-26;
        const corrected_remainder = (inverse_ln2_high * value_wide - rounded) +
            inverse_ln2_low * value_wide;
        const corrected_squared = corrected_remainder * corrected_remainder;
        result = (scale - 1.0) + scale * corrected_remainder *
            ((accurate_c[0] + corrected_remainder * accurate_c[1]) +
                corrected_squared * ((accurate_c[2] + corrected_remainder * accurate_c[3]) +
                    corrected_squared * (accurate_c[4] + corrected_remainder * accurate_c[5])));
        upper = @floatCast(result);
    }
    return upper;
}
