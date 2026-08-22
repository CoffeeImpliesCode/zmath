// Correctly rounded cube roots for binary16, binary32, and binary64.
//
// The finite kernels are ported from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa:
//   src/binary16/cbrt/cbrtf16.c (Copyright 2025 Maxence Ponsardin).
//   src/binary32/cbrt/cbrtf.c (Copyright 2023 Alexei Sibidanov).
//   src/binary64/cbrt/cbrt.c (Copyright 2021-2022 Alexei Sibidanov).
//
// The upstream code is distributed under the MIT License:
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

const sign_mask_f: u32 = 0x8000_0000;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const infinity_bits_f: u32 = 0x7f80_0000;
const quiet_nan_bit_f: u32 = 0x0040_0000;

const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;

const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const quiet_nan_bit_h: u16 = 0x0200;

fn absD(value: f64) f64 {
    return if (value < 0.0) -value else value;
}

fn copySignD(value: f64, sign: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    const sign_bits: u64 = @bitCast(sign);
    return @bitCast((value_bits & ~sign_mask_d) | (sign_bits & sign_mask_d));
}

pub fn cbrth(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;
    if (magnitude > infinity_bits_h) return @bitCast(bits | quiet_nan_bit_h);
    if (magnitude == infinity_bits_h or magnitude == 0) return value;

    // CORE-MATH's binary16 kernel resolves its exact cases directly and uses
    // the correctly-rounded binary32 kernel for the remaining finite inputs.
    const value_f: f32 = @floatCast(value);
    const value_f_bits: u32 = @bitCast(value_f);
    const table_index: usize = @intCast((value_f_bits >> 19) % 16);
    const exact_mantissas = [_]u16{
        0x00, 0xff, 0xff, 0xff,
        0x33, 0x1c, 0x32, 0xff,
        0xff, 0xff, 0xff, 0x00,
        0xff, 0xff, 0xff, 0x10,
    };
    const exact_exponents = [_]i32{
        0, 0, 0, 0,
        1, 2, 0, 0,
        0, 0, 0, 1,
        0, 0, 0, 0,
    };
    const exact_results = [_]u16{
        0x3c00, 0,      0,      0,
        0x3d80, 0x3f00, 0x3c80, 0,
        0,      0,      0,      0x3e00,
        0,      0,      0,      0x3d00,
    };

    const mantissa_match = exact_mantissas[table_index] == @as(u16, @intCast((value_f_bits >> 13) % 64));
    if (mantissa_match) {
        const exponent: i32 = @intCast((value_f_bits & magnitude_mask_f) >> 23);
        const residue: i32 = @mod(exponent + 2, 3);
        if (exact_exponents[table_index] == residue) {
            const scaled_exponent: i32 = @divTrunc(exponent - 127 - exact_exponents[table_index], 3);
            const result_exponent: u16 = @truncate(@as(u32, @bitCast(scaled_exponent)));
            const result_bits = (bits & sign_mask_h) +% (result_exponent << 10) +% exact_results[table_index];
            return @bitCast(result_bits);
        }
    }

    if ((bits & 0x03ff) == 0x0151) {
        const exponent: i32 = @intCast((bits & magnitude_mask_h) >> 10);
        if (@mod(exponent, 3) == 1 and exponent < 31) {
            var corrected_bits: u16 = @intCast((@divTrunc(exponent - 16, 3) + 15) << 10);
            corrected_bits +%= 0x018b;
            corrected_bits +%= (bits >> 15) << 15;
            const corrected: f16 = @bitCast(corrected_bits);
            const sign_adjustment: f32 = @as(f32, @floatFromInt(corrected_bits >> 15)) - 0.5;
            const correction: f32 = 0x1p-16 * sign_adjustment;
            return @floatCast(@as(f32, corrected) + correction);
        }
    }

    return @floatCast(cbrtf(value_f));
}

pub fn cbrtf(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask_f;
    if (magnitude > infinity_bits_f) return @bitCast(bits | quiet_nan_bit_f);
    if (magnitude == infinity_bits_f or magnitude == 0) return value;

    const escale = [_]f64{
        1.0,
        0x1.428a2f98d728bp+0,
        0x1.965fea53d6e3dp+0,
    };
    var shifted_bits = bits << 1;
    const sign: u32 = bits >> 31;
    var exponent: u32 = shifted_bits >> 24;
    if (shifted_bits < (1 << 24)) {
        const normalization_shift: u32 = @clz(shifted_bits) - 7;
        shifted_bits <<= @intCast(normalization_shift);
        exponent -%= normalization_shift - 1;
    }

    const mantissa = shifted_bits & 0x00ff_ffff;
    const normalized_bits = (@as(u64, mantissa) << 28) | (@as(u64, 0x3ff) << 52);
    const normalized: f64 = @bitCast(normalized_bits);
    exponent +%= 899;
    const exponent_third: u32 = exponent / 3;
    const exponent_remainder: u32 = exponent % 3;

    var scale_bits: u64 = @bitCast(escale[exponent_remainder]);
    scale_bits +%= (@as(u64, exponent_third) -% 342) << 52;
    scale_bits |= @as(u64, sign) << 63;
    const scale: f64 = @bitCast(scale_bits);

    const coefficients = [_]f64{
        0x1.2319d352ea5d5p-1,
        0x1.67ad8ee258d1ap-1,
        -0x1.9342edf9cbad9p-2,
        0x1.b6388fc510a75p-3,
        -0x1.6002455599e2fp-4,
        0x1.7b096936192c4p-6,
        -0x1.e5577187e8bf8p-9,
        0x1.169ef81d6c34ep-12,
    };
    const reciprocal_correction = -0x1.9931c6c2d19d1p-6 / normalized;
    const normalized_squared = normalized * normalized;
    const normalized_fourth = normalized_squared * normalized_squared;
    const polynomial = ((coefficients[0] + normalized * coefficients[1]) +
        normalized_squared * (coefficients[2] + normalized * coefficients[3])) +
        normalized_fourth * ((coefficients[4] + normalized * coefficients[5]) +
            normalized_squared * (coefficients[6] + normalized * coefficients[7])) + reciprocal_correction;
    var result = polynomial * scale;
    const upper: f32 = @floatCast(result);
    const lower: f32 = @floatCast(result - scale * 1.4182e-9);
    if (upper == lower) return upper;

    const correction_constant = -0x1.ab16ec65d138fp+3;
    const error_value = polynomial * polynomial * polynomial - normalized;
    const corrected_polynomial = polynomial -
        (polynomial * reciprocal_correction * correction_constant) * error_value;
    result = corrected_polynomial * scale;

    var result_bits: u64 = @bitCast(result);
    var rounded_result: f64 = result;
    const m0: i64 = @bitCast(result_bits << 19);
    const m1: i64 = m0 >> 63;
    const near_rounding_boundary = @as(u64, @bitCast(m0 ^ m1)) < (@as(u64, 1) << 31);
    if (near_rounding_boundary) {
        result_bits = (result_bits + (@as(u64, 1) << 31)) & 0xffff_ffff_0000_0000;
        rounded_result = @bitCast(result_bits);
    }
    return @floatCast(rounded_result);
}

pub fn cbrtd(value: f64) f64 {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    const magnitude = bits & magnitude_mask_d;
    if (magnitude > infinity_bits_d) return @bitCast(bits | quiet_nan_bit_d);
    if (magnitude == infinity_bits_d or magnitude == 0) return value;

    const escale = [_]f64{
        1.0,
        0x1.428a2f98d728bp+0,
        0x1.965fea53d6e3dp+0,
    };
    const coefficients = [_]f64{
        0x1.1b0babccfef9cp-1,
        0x1.2c9a3e94d1da5p-1,
        -0x1.4dc30b1a1ddbap-3,
        0x1.7a8d3e4ec9b07p-6,
    };
    const correction_u0 = 0x1.5555555555555p-2;
    const correction_u1 = 0x1.c71c71c71c71cp-3;
    const reciprocal_scales = [_]f64{ 1.0, -1.0, 0.5, -0.5, 0.25, -0.25 };

    const mantissa_mask: u64 = 0x000f_ffff_ffff_ffff;
    var mantissa = bits & mantissa_mask;
    var exponent: u64 = (bits >> 52) & 0x7ff;
    if (exponent == 0) {
        const absolute_bits = bits & magnitude_mask_d;
        const normalization_shift: u64 = @clz(absolute_bits) - 11;
        mantissa <<= @intCast(normalization_shift);
        mantissa &= mantissa_mask;
        exponent -%= normalization_shift - 1;
    }
    exponent +%= 3072;

    const exponent_third: u64 = exponent / 3;
    const exponent_remainder: u64 = exponent % 3;
    const normalized_bits = mantissa | (@as(u64, 0x3ff) << 52);
    var scaled_normalized_bits = normalized_bits;
    scaled_normalized_bits +%= exponent_remainder << 52;
    scaled_normalized_bits |= (bits >> 63) << 63;
    const scaled_normalized: f64 = @bitCast(scaled_normalized_bits);

    var scale_bits: u64 = @bitCast(escale[exponent_remainder]);
    scale_bits |= (bits >> 63) << 63;
    const scale: f64 = @bitCast(scale_bits);

    const normalized: f64 = @bitCast(normalized_bits);
    const sign: u64 = bits >> 63;
    const reciprocal = 1.0 / normalized;
    const reciprocal_scale = reciprocal_scales[(exponent_remainder << 1) | sign];
    const reciprocal_scaled = reciprocal * reciprocal_scale;
    const normalized_squared = normalized * normalized;
    const coefficient_0 = coefficients[0] + normalized * coefficients[1];
    const coefficient_2 = coefficients[2] + normalized * coefficients[3];
    var approximation = coefficient_0 + normalized_squared * coefficient_2;
    const approximation_squared = approximation * approximation;
    var error_value = approximation_squared * (approximation * reciprocal) - 1.0;
    approximation -= (error_value * approximation) * (correction_u0 - correction_u1 * error_value);
    approximation *= scale;

    var approximation_squared_high = approximation * approximation;
    const approximation_squared_low = @mulAdd(f64, approximation, approximation, -approximation_squared_high);
    const approximation_cubed_high = approximation_squared_high * approximation;
    const approximation_cubed_low = @mulAdd(f64, approximation, approximation_squared_high, -approximation_cubed_high) +
        approximation * approximation_squared_low;
    error_value = ((approximation_cubed_high - scaled_normalized) + approximation_cubed_low) *
        reciprocal_scaled;
    var correction = error_value * (approximation * correction_u0);
    var rounded = approximation - correction;
    correction = (approximation - rounded) - correction;
    var adjustment = absD(correction);
    var adjustment_zero = absD(adjustment - 0x1p-53);
    var adjustment_one = absD(adjustment - (0x1p-52 + 0x1p-53));

    if (adjustment_zero < 0x1p-75 or adjustment_one < 0x1p-75) {
        approximation_squared_high = rounded * rounded;
        const rounded_squared_low = @mulAdd(f64, rounded, rounded, -approximation_squared_high);
        const rounded_cubed_high = approximation_squared_high * rounded;
        const rounded_cubed_low = @mulAdd(f64, rounded, approximation_squared_high, -rounded_cubed_high) +
            rounded * rounded_squared_low;
        error_value = ((rounded_cubed_high - scaled_normalized) + rounded_cubed_low) *
            reciprocal_scaled;
        const next_correction = error_value * (rounded * correction_u0);
        const next_rounded = rounded - next_correction;
        correction = (rounded - next_rounded) - next_correction;
        rounded = next_rounded;
        adjustment = absD(correction);
        adjustment_zero = absD(adjustment - 0x1p-53);
        adjustment_one = absD(adjustment - (0x1p-52 + 0x1p-53));

        if (adjustment_zero < 0x1p-98 or adjustment_one < 0x1p-98) {
            const absolute_scaled_normalized = absD(scaled_normalized);
            if (absolute_scaled_normalized == 0x1.9b78223aa307cp+1) {
                rounded = copySignD(0x1.79d15d0e8d59cp+0, scaled_normalized);
            }
            if (absolute_scaled_normalized == 0x1.a202bfc89ddffp+2) {
                rounded = copySignD(0x1.de87aa837820fp+0, scaled_normalized);
            }
        }
    }

    var result_bits: u64 = @bitCast(rounded);
    result_bits +%= (@as(u64, exponent_third) -% 1365) << 52;
    const boundary_bits = result_bits << 30;
    const boundary_signed: i64 = @bitCast(boundary_bits);
    const boundary_sign: i64 = boundary_signed >> 63;
    const near_boundary = @as(u64, @bitCast(boundary_signed ^ boundary_sign)) <= (@as(u64, 1) << 30);
    if (near_boundary) {
        const rounded_candidate_bits =
            (@as(u64, @bitCast(rounded)) + (@as(u64, 1) << 15)) & 0xffff_ffff_ffff_0000;
        const rounded_candidate: f64 = @bitCast(rounded_candidate_bits);
        const rounded_result_bits = (result_bits + (@as(u64, 1) << 15)) & 0xffff_ffff_ffff_0000;
        if (absD((rounded_candidate - rounded) - correction) < 0x1p-60 or absD(scaled_normalized) == 1.0) {
            result_bits = rounded_result_bits;
        }
    }
    return @bitCast(result_bits);
}
