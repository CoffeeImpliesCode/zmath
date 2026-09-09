const magnitude_mask: u32 = 0x7fff_ffff;
const infinity_bits: u32 = 0x7f80_0000;
const quiet_nan_bit: u32 = 0x0040_0000;

noinline fn exp2m1fSmall(value: f32, input_bits: u32, absolute_bits: u32) f32 {
    @setFloatMode(.strict);
    const value_f64: f64 = value;
    // SAFETY: Every magnitude branch initializes the polynomial before it is used.
    var polynomial: f64 = undefined;
    const value_squared = value_f64 * value_f64;
    if (absolute_bits < 0x3d67_a4cc) {
        if (absolute_bits < 0x3caa_2fee) {
            if (absolute_bits < 0x3bac_1405) {
                if (absolute_bits < 0x3a35_8876) {
                    if (absolute_bits < 0x37d3_2ef6) {
                        if (absolute_bits < 0x331f_dd82) {
                            if (absolute_bits < 0x2538_aa3b) {
                                polynomial = 0x1.62e42fefa39efp-1;
                            } else {
                                polynomial = 0x1.62e42fefa39fp-1 +
                                    value_f64 * 0x1.ebfbdff82c58fp-3;
                            }
                        } else {
                            if (input_bits == 0xb3d8_5005)
                                return @floatCast(-0x1.2bdf76p-24 - 0x1.8p-77);
                            if (input_bits == 0x3338_428d)
                                return @floatCast(0x1.fee08ap-26 + 0x1p-80);
                            const c = [3]f64{
                                0x1.62e42fefa39efp-1,
                                0x1.ebfbdff8548fdp-3,
                                0x1.c6b08d704a06dp-5,
                            };
                            polynomial = c[0] + value_f64 * (c[1] + value_f64 * c[2]);
                        }
                    } else {
                        if (input_bits == 0x388b_ca4f)
                            return @floatCast(0x1.839702p-15 - 0x1.8p-68);
                        const c = [4]f64{
                            0x1.62e42fefa39efp-1,
                            0x1.ebfbdff82c58fp-3,
                            0x1.c6b08dc82b347p-5,
                            0x1.3b2ab6fbad172p-7,
                        };
                        polynomial = (c[0] + value_f64 * c[1]) +
                            value_squared * (c[2] + value_f64 * c[3]);
                    }
                } else {
                    const c = [5]f64{
                        0x1.62e42fefa39efp-1,
                        0x1.ebfbdff82c068p-3,
                        0x1.c6b08d704a6dcp-5,
                        0x1.3b2ac262c3eedp-7,
                        0x1.5d87fe7af779ap-10,
                    };
                    polynomial = (c[0] + value_f64 * c[1]) + value_squared *
                        (c[2] + value_f64 * (c[3] + value_f64 * c[4]));
                }
            } else {
                const c = [6]f64{
                    0x1.62e42fefa39fp-1,
                    0x1.ebfbdff82c58dp-3,
                    0x1.c6b08d7011d13p-5,
                    0x1.3b2ab6fbd267dp-7,
                    0x1.5d88a81cea49ep-10,
                    0x1.430912ea9b963p-13,
                };
                polynomial = (c[0] + value_f64 * c[1]) + value_squared *
                    ((c[2] + value_f64 * c[3]) +
                        value_squared * (c[4] + value_f64 * c[5]));
            }
        } else {
            const c = [7]f64{
                0x1.62e42fefa39efp-1,
                0x1.ebfbdff82c639p-3,
                0x1.c6b08d7049f1cp-5,
                0x1.3b2ab6f5243bdp-7,
                0x1.5d87fe80a9e6cp-10,
                0x1.430d0b9257fa8p-13,
                0x1.ffcbfc4cf0952p-17,
            };
            polynomial = (c[0] + value_f64 * c[1]) + value_squared *
                ((c[2] + value_f64 * c[3]) + value_squared *
                    (c[4] + value_f64 * (c[5] + value_f64 * c[6])));
        }
    } else {
        const c = [8]f64{
            0x1.62e42fefa39efp-1,
            0x1.ebfbdff82c591p-3,
            0x1.c6b08d704cf6bp-5,
            0x1.3b2ab6fba00cep-7,
            0x1.5d87fdfdaadb4p-10,
            0x1.4309137333066p-13,
            0x1.ffe5e90daf7ddp-17,
            0x1.62c0220eed731p-20,
        };
        polynomial = ((c[0] + value_f64 * c[1]) +
            value_squared * (c[2] + value_f64 * c[3])) +
            (value_squared * value_squared) *
                ((c[4] + value_f64 * c[5]) +
                    value_squared * (c[6] + value_f64 * c[7]));
    }
    return @floatCast(polynomial * value_f64);
}

/// Correctly-rounded base-2 exponential biased by one for binary32 values.
///
/// Copyright (c) 2022-2025 Alexei Sibidanov.
///
/// This file is part of the CORE-MATH project
/// (https://core-math.gitlabpages.inria.fr/).
///
/// Ported from src/binary32/exp2m1/exp2m1f.c at revision
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
pub fn exp2m1f(value: f32) f32 {
    @setFloatMode(.strict);

    const input_bits: u32 = @bitCast(value);
    const absolute_bits = input_bits & magnitude_mask;
    if (absolute_bits > infinity_bits) return @bitCast(input_bits | quiet_nan_bit);
    if (absolute_bits == infinity_bits) {
        if (input_bits == infinity_bits) return value;
        return -1.0;
    }
    if (absolute_bits == 0) return value;

    const value_f64: f64 = value;
    if (input_bits >= 0xc1c8_0000) {
        return -1.0 + 0x1p-26;
    }
    if (absolute_bits >= 0x4300_0000) {
        return 0x1.fffffep+127 + 0x1.fffffep+127;
    }
    if (absolute_bits < 0x3df9_5f1f) return exp2m1fSmall(value, input_bits, absolute_bits);

    const c = [6]f64{
        0x1.62e42fefa398bp-5,
        0x1.ebfbdff84555ap-11,
        0x1.c6b08d4ad86d3p-17,
        0x1.3b2ad1b1716a2p-23,
        0x1.5d7472718ce9dp-30,
        0x1.4a1d7f457ac56p-37,
    };
    const table = [16]f64{
        0x1p+0,
        0x1.0b5586cf9890fp+0,
        0x1.172b83c7d517bp+0,
        0x1.2387a6e756238p+0,
        0x1.306fe0a31b715p+0,
        0x1.3dea64c123422p+0,
        0x1.4bfdad5362a27p+0,
        0x1.5ab07dd485429p+0,
        0x1.6a09e667f3bcdp+0,
        0x1.7a11473eb0187p+0,
        0x1.8ace5422aa0dap+0,
        0x1.9c49182a3f09p+0,
        0x1.ae89f995ad3adp+0,
        0x1.c199bdd85529cp+0,
        0x1.d5818dcfba487p+0,
        0x1.ea4afa2a490dap+0,
    };
    const scaled = 16.0 * value_f64;
    const integral = @floor(scaled);
    const integral_bits: u64 = @bitCast(integral + 0x1.8p52);
    const remainder = scaled - integral;
    const remainder_squared = remainder * remainder;
    const table_index: u4 = @truncate(integral_bits);
    const exponent_bits = ((integral_bits >> 4) +% 0x3ff) << 52;
    const scale: f64 = table[table_index] * @as(f64, @bitCast(exponent_bits));
    const c0 = (c[0] + remainder * c[1]) + remainder_squared *
        ((c[2] + remainder * c[3]) + remainder_squared * (c[4] + remainder * c[5]));
    return @floatCast((scale - 1.0) + (scale * remainder) * c0);
}
