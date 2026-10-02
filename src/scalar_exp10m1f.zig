const fma_impl = @import("scalar_fma.zig");

const sign_mask: u32 = 0x8000_0000;
const magnitude_mask: u32 = 0x7fff_ffff;
const infinity_bits: u32 = 0x7f80_0000;
const quiet_nan_bit: u32 = 0x0040_0000;

const coefficients: [6]f64 = .{
    0x1.62e42fefa398bp-5,
    0x1.ebfbdff84555ap-11,
    0x1.c6b08d4ad86d3p-17,
    0x1.3b2ad1b1716a2p-23,
    0x1.5d7472718ce9dp-30,
    0x1.4a1d7f457ac56p-37,
};

const table: [16]f64 = .{
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

const exceptional_rounding: [2][2]f32 = .{
    .{ 0x1.fffffep+127, 0x1.fffffep+127 },
    .{ -1.0, 0x1p-26 },
};

/// Correctly-rounded base-10 exponential function biased by one for binary32 values.
///
/// Copyright (c) 2022-2025 Alexei Sibidanov, Paul Zimmermann.
///
/// This file is part of the CORE-MATH project
/// (https://core-math.gitlabpages.inria.fr/).
///
/// Ported from src/binary32/exp10m1/exp10m1f.c at revision
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
pub fn exp10m1f(value: f32) f32 {
    @setFloatMode(.strict);

    const input_bits: u32 = @bitCast(value);
    const absolute_bits = input_bits & magnitude_mask;
    if (absolute_bits > infinity_bits) return @bitCast(input_bits | quiet_nan_bit);
    if (absolute_bits == infinity_bits) {
        if ((input_bits & sign_mask) == 0) return value;
        return exceptional_rounding[1][0];
    }

    const value_wide: f64 = value;
    if (input_bits > 0xc0f0_d2f1) {
        return exceptional_rounding[1][0] + exceptional_rounding[1][1];
    }
    if (absolute_bits > 0x421a_209a) {
        return exceptional_rounding[0][0] + exceptional_rounding[0][1];
    }

    if (absolute_bits < 0x3d89_c604) {
        const squared = value_wide * value_wide;
        var result_wide: f64 = undefined;
        if (absolute_bits < 0x3d16_22fb) {
            if (absolute_bits < 0x3c8b_76a3) {
                if (absolute_bits < 0x3bcc_ed04) {
                    if (absolute_bits < 0x3acf_33eb) {
                        if (absolute_bits < 0x395a_966b) {
                            if (absolute_bits < 0x36fe_4a4b) {
                                if (absolute_bits < 0x3240_7f39) {
                                    if (absolute_bits < 0x245e_5bd9) {
                                        result_wide = 0x1.26bb1bbb55516p+1;
                                    } else {
                                        if (input_bits == 0x2c99_4b7b) {
                                            return @as(f32, 0x1.60f974p-37) -
                                                @as(f32, 0x1p-90);
                                        }
                                        result_wide = fma_impl.fma64(
                                            value_wide,
                                            0x1.53524c73cea69p+1,
                                            0x1.26bb1bbb55516p+1,
                                        );
                                    }
                                } else {
                                    if (input_bits == 0xb6fa_215b) {
                                        return -@as(f32, 0x1.1ff87ep-16) +
                                            @as(f32, 0x1p-68);
                                    }
                                    const inner = fma_impl.fma64(
                                        value_wide,
                                        0x1.0470591de2c75p+1,
                                        0x1.53524c73cea62fp+1,
                                    );
                                    result_wide = fma_impl.fma64(
                                        value_wide,
                                        inner,
                                        0x1.26bb1bbb55516p+1,
                                    );
                                }
                            } else {
                                const cp = [4]f64{
                                    0x1.26bb1bbb55515p+1,
                                    0x1.53524c73cea69p+1,
                                    0x1.0470595038cc2p+1,
                                    0x1.2bd7609fe1561p+0,
                                };
                                const low = fma_impl.fma64(value_wide, cp[1], cp[0]);
                                const high = fma_impl.fma64(value_wide, cp[3], cp[2]);
                                result_wide = fma_impl.fma64(squared, high, low);
                            }
                        } else {
                            const cp = [5]f64{
                                0x1.26bb1bbb55516p+1,
                                0x1.53524c73ce6dbp+1,
                                0x1.0470591de3024p+1,
                                0x1.2bd76b79060e6p+0,
                                0x1.1429ffd3a963dp-1,
                            };
                            const low = fma_impl.fma64(value_wide, cp[1], cp[0]);
                            const tail = fma_impl.fma64(value_wide, cp[4], cp[3]);
                            const high = fma_impl.fma64(value_wide, tail, cp[2]);
                            result_wide = fma_impl.fma64(squared, high, low);
                        }
                    } else {
                        const cp = [6]f64{
                            0x1.26bb1bbb55516p+1,
                            0x1.53524c73cea67p+1,
                            0x1.0470591dc2953p+1,
                            0x1.2bd760a004d64p+0,
                            0x1.142a85da6f072p-1,
                            0x1.a7ed70725b00ep-3,
                        };
                        const low = fma_impl.fma64(value_wide, cp[1], cp[0]);
                        const middle = fma_impl.fma64(value_wide, cp[3], cp[2]);
                        const high = fma_impl.fma64(value_wide, cp[5], cp[4]);
                        const upper = fma_impl.fma64(squared, high, middle);
                        result_wide = fma_impl.fma64(squared, upper, low);
                    }
                } else {
                    const cp = [7]f64{
                        0x1.26bb1bbb55516p+1,
                        0x1.53524c73ceadep+1,
                        0x1.0470591de2bb4p+1,
                        0x1.2bd76099a9d33p+0,
                        0x1.1429ffd829b0bp-1,
                        0x1.a7f2a6a0f7dc8p-3,
                        0x1.16e4dfbce0f56p-4,
                    };
                    const low = fma_impl.fma64(value_wide, cp[1], cp[0]);
                    const middle = fma_impl.fma64(value_wide, cp[3], cp[2]);
                    const tail = fma_impl.fma64(value_wide, cp[6], cp[5]);
                    const high = fma_impl.fma64(value_wide, tail, cp[4]);
                    const upper = fma_impl.fma64(squared, high, middle);
                    result_wide = fma_impl.fma64(squared, upper, low);
                }
            } else {
                const cp = [8]f64{
                    0x1.26bb1bbb55515p+1,
                    0x1.53524c73cea6ap+1,
                    0x1.0470591de476p+1,
                    0x1.2bd7609fd4ee2p+0,
                    0x1.1429ff70a9b48p-1,
                    0x1.a7ed71259ba5bp-3,
                    0x1.16f3004fb3ac1p-4,
                    0x1.4116b0388aa9fp-6,
                };
                const fourth = squared * squared;
                const low_0 = fma_impl.fma64(value_wide, cp[1], cp[0]);
                const low_2 = fma_impl.fma64(value_wide, cp[3], cp[2]);
                const low = fma_impl.fma64(squared, low_2, low_0);
                const high_0 = fma_impl.fma64(value_wide, cp[5], cp[4]);
                const high_2 = fma_impl.fma64(value_wide, cp[7], cp[6]);
                const high = fma_impl.fma64(squared, high_2, high_0);
                result_wide = fma_impl.fma64(fourth, high, low);
            }
        } else {
            const cp = [9]f64{
                0x1.26bb1bbb55515p+1,
                0x1.53524c73cea42p+1,
                0x1.0470591de2d1dp+1,
                0x1.2bd760a010a53p+0,
                0x1.1429ffd16170cp-1,
                0x1.a7ed6b2a0d97fp-3,
                0x1.16e4e37fa51e4p-4,
                0x1.4147fe4c1676fp-6,
                0x1.4897c4b3e329ap-8,
            };
            const fourth = squared * squared;
            const low_0 = fma_impl.fma64(value_wide, cp[1], cp[0]);
            const low_2 = fma_impl.fma64(value_wide, cp[3], cp[2]);
            const low = fma_impl.fma64(squared, low_2, low_0);
            const high_0 = fma_impl.fma64(value_wide, cp[5], cp[4]);
            const high_tail = fma_impl.fma64(value_wide, cp[8], cp[7]);
            const high_2 = fma_impl.fma64(value_wide, high_tail, cp[6]);
            const high = fma_impl.fma64(squared, high_2, high_0);
            result_wide = fma_impl.fma64(fourth, high, low);
        }
        result_wide *= value_wide;
        return @floatCast(result_wide);
    }

    if ((input_bits << 11) == 0) {
        const power_index = (input_bits >> 21) -% 0x1fc;
        if (power_index <= 0xb) {
            if (power_index == 0) return @as(f32, 10.0) - 1.0;
            if (power_index == 4) return @as(f32, 100.0) - 1.0;
            if (power_index == 6) return @as(f32, 1000.0) - 1.0;
            if (power_index == 8) return @as(f32, 10_000.0) - 1.0;
            if (power_index == 9) return @as(f32, 100_000.0) - 1.0;
            if (power_index == 10) return @as(f32, 1_000_000.0) - 1.0;
            if (power_index == 11) return @as(f32, 10_000_000.0) - 1.0;
        }
    }

    const log2_10_high: f64 = 0x1.a934f09p+1 * 16.0;
    const log2_10_low: f64 = 0x1.e68dc57f2496p-29 * 16.0;
    const scaled = log2_10_high * value_wide;
    const integral = @floor(scaled);
    const remainder = fma_impl.fma64(log2_10_low, value_wide, scaled - integral);
    // The finite guards bound integral to [-400, 2047]. The bias sum is exact,
    // and its low 32 encoding bits hold the signed integer in two's complement.
    const integral_bits: u64 = @bitCast(integral + 0x1.8p52);
    const integral_integer: i64 = @as(i32, @bitCast(@as(u32, @truncate(integral_bits))));
    const table_index: u4 = @truncate(@as(u64, @bitCast(integral_integer)));
    const table_bits: u64 = @bitCast(table[table_index]);
    // Each table value is in [1, 2), and the exponent is in [-25, 127].
    // Adjusting its exponent is the same exact normal scaling as multiplication.
    const exponent_bits = @as(u64, @bitCast(integral_integer >> 4)) << 52;
    const scale: f64 = @bitCast(table_bits +% exponent_bits);

    const remainder_squared = remainder * remainder;
    const c0 = fma_impl.fma64(remainder, coefficients[1], coefficients[0]);
    const c2 = fma_impl.fma64(remainder, coefficients[3], coefficients[2]);
    const c4 = fma_impl.fma64(remainder, coefficients[5], coefficients[4]);
    const polynomial = fma_impl.fma64(
        remainder_squared,
        fma_impl.fma64(remainder_squared, c4, c2),
        c0,
    );
    const correction = scale * remainder;
    const result_wide = fma_impl.fma64(correction, polynomial, scale - 1.0);
    return @floatCast(result_wide);
}
