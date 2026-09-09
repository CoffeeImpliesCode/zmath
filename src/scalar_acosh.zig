// Correctly rounded inverse hyperbolic cosine for binary16, binary32, and binary64.
//
// The finite kernels are ported from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa:
//   src/binary16/acosh/acoshf16.c (Copyright 2025 Paul Zimmermann).
//   src/binary32/acosh/acoshf.c (Copyright 2023 Alexei Sibidanov).
//   src/binary64/acosh/acosh.c (Copyright 2023 Alexei Sibidanov).
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

const acoshf_data = @import("scalar_acoshf_data.zig");

const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const quiet_nan_bit_h: u16 = 0x0200;
const canonical_nan_bits_h: u16 = 0x7e00;

const acoshh_polynomials: [16][6]f32 = .{
    .{
        0x1.91940ap+0, -0x1.a1ee52p-3,  0x1.0dc64ap-4, -0x1.453108p-6,
        0x1.040a46p-8, -0x1.803198p-12,
    },
    .{
        0x1.8eb5fap+0,  -0x1.6c52a2p-3,  0x1.4dddcp-5, -0x1.e9e95cp-8,
        0x1.b697b2p-11, -0x1.5b2df4p-15,
    },
    .{
        0x1.854f5ap+0, -0x1.11e128p-3,  0x1.362606p-6, -0x1.fb23bp-10,
        0x1.e2d03p-14, -0x1.8d21eap-19,
    },
    .{
        0x1.717322p+0,  -0x1.608d74p-4,  0x1.c3ceb2p-8, -0x1.8859d4p-12,
        0x1.828e08p-17, -0x1.44eab2p-23,
    },
    .{
        0x1.52b1b4p+0,  -0x1.8e60acp-5,  0x1.12e502p-9, -0x1.ee8866p-15,
        0x1.f10d92p-21, -0x1.a71f3ep-28,
    },
    .{
        0x1.1b4176p+1,   0x1.c58fdp-4,   -0x1.3d491cp-9, 0x1.244e58p-15,
        -0x1.2b59e6p-22, 0x1.029682p-30,
    },
    .{
        0x1.744536p+1,   0x1.c4d402p-5,  -0x1.3c787ap-11, 0x1.235b58p-18,
        -0x1.2a375cp-26, 0x1.017ceap-35,
    },
    .{
        0x1.cd16eap+1,   0x1.c490cp-6,   -0x1.3c27dep-13, 0x1.22f70cp-21,
        -0x1.29b90ep-30, 0x1.00fd26p-40,
    },
    .{
        0x1.12ec74p+2,   0x1.c475e4p-7,  -0x1.3c05dap-15, 0x1.22cb2cp-24,
        -0x1.29804ap-34, 0x1.00c29ap-45,
    },
    .{
        0x1.3f4c7ep+2,  0x1.c45ddap-8,  -0x1.3be5p-17, 0x1.229e72p-27,
        -0x1.2943fp-38, 0x1.008224p-50,
    },
    .{
        0x1.6baa5p+2,    0x1.c45576p-9, -0x1.3bd988p-19, 0x1.228f6ep-30,
        -0x1.29308ep-42, 0x1.006ddp-55,
    },
    .{
        0x1.980882p+2,   0x1.c449a8p-10, -0x1.3bc944p-21, 0x1.2278ep-33,
        -0x1.2911bap-46, 0x1.004d3p-60,
    },
    .{
        0x1.c46664p+2,   0x1.c44006p-11, -0x1.3bbbdp-23, 0x1.22665ep-36,
        -0x1.28f8a4p-50, 0x1.00326p-65,
    },
    .{
        0x1.f0c4d6p+2,   0x1.c43246p-12, -0x1.3ba86ap-25, 0x1.224b62p-39,
        -0x1.28d39ep-54, 0x1.000a46p-70,
    },
    .{
        0x1.0e912cp+3,  0x1.c42b4ep-13, -0x1.3b9eb4p-27, 0x1.223e1p-42,
        -0x1.28c1ap-58, 0x1.ffee3cp-76,
    },
    .{
        0x1.24bfccp+3,   0x1.c4263p-14,  -0x1.3b97b8p-29, 0x1.22347ap-45,
        -0x1.28b4c4p-62, 0x1.ffd2ecp-81,
    },
};

/// Correctly rounded binary16 inverse hyperbolic cosine under caller round-to-nearest-even.
pub fn acoshh(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;
    if (magnitude > infinity_bits_h) return @bitCast(bits | quiet_nan_bit_h);
    if (magnitude == infinity_bits_h) {
        if (bits == infinity_bits_h) return value;
        return @bitCast(canonical_nan_bits_h);
    }
    if ((bits & 0x8000) != 0 or magnitude < 0x3c00) return @bitCast(canonical_nan_bits_h);
    if (bits == 0x3c00) return 0.0;

    const input: f32 = @floatCast(value);
    const input_squared = input * input;
    const exponent_index: usize = @intCast((magnitude >> 10) - 15);
    const coefficients = acoshh_polynomials[exponent_index];
    const coefficient_4 = @mulAdd(f32, coefficients[5], input, coefficients[4]);
    const coefficient_2 = @mulAdd(f32, coefficients[3], input, coefficients[2]);
    const coefficient_2_adjusted = @mulAdd(f32, coefficient_4, input_squared, coefficient_2);
    const coefficient_0 = @mulAdd(f32, coefficients[1], input, coefficients[0]);
    var result = @mulAdd(f32, coefficient_2_adjusted, input_squared, coefficient_0);
    if (exponent_index < 5) result = @sqrt(input - 1.0) * result;

    switch (bits) {
        0x4211 => result = @bitCast(@as(u32, 0x3fe3_1ffe)),
        0x48e0 => result = @bitCast(@as(u32, 0x403d_f003)),
        0x4afb => result = @bitCast(@as(u32, 0x4055_0001)),
        0x4dff => result = @bitCast(@as(u32, 0x4077_aff6)),
        0x4e02 => result = @bitCast(@as(u32, 0x4077_cffb)),
        0x4ede => result = @bitCast(@as(u32, 0x4080_300b)),
        0x4e9e => result = @bitCast(@as(u32, 0x407e_0015)),
        0x51ba => result = @bitCast(@as(u32, 0x4090_900e)),
        0x522e => result = @bitCast(@as(u32, 0x4093_0009)),
        0x5448 => result = @bitCast(@as(u32, 0x409d_700c)),
        0x54eb => result = @bitCast(@as(u32, 0x40a1_dffa)),
        0x59e9 => result = @bitCast(@as(u32, 0x40bd_f003)),
        0x5c42 => result = @bitCast(@as(u32, 0x40c9_a002)),
        0x5fcd => result = @bitCast(@as(u32, 0x40dd_0007)),
        0x61d8 => result = @bitCast(@as(u32, 0x40e9_f006)),
        0x6dbc => result = @bitCast(@as(u32, 0x4115_f004)),
        0x739d => result = @bitCast(@as(u32, 0x4125_9004)),
        else => {},
    }
    return @floatCast(result);
}

const magnitude_mask_f: u32 = 0x7fff_ffff;
const infinity_bits_f: u32 = 0x7f80_0000;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const canonical_nan_bits_f: u32 = 0x7fc0_0000;

/// Correctly rounded binary32 inverse hyperbolic cosine under caller round-to-nearest-even.
pub fn acoshf(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask_f;
    if (magnitude > infinity_bits_f) return @bitCast(bits | quiet_nan_bit_f);
    if (magnitude == infinity_bits_f) {
        if (bits == infinity_bits_f) return value;
        return @bitCast(canonical_nan_bits_f);
    }
    if ((bits & 0x8000_0000) != 0 or magnitude < 0x3f80_0000) return @bitCast(canonical_nan_bits_f);
    if (bits == 0x3f80_0000) return 0.0;

    if (bits <= 0x3f99_db23) {
        const zf = value - 1.0;
        const z: f64 = @floatCast(zf);
        const root = @sqrt(2.0 * z);
        const coefficients = [_]f64{
            -0x1.555555555491ep-4,
            0x1.333333319c2p-6,
            -0x1.6db6da26e5e0ep-8,
            0x1.f1c698b7100cep-10,
            -0x1.6e7ca9e152dc2p-11,
            0x1.1b4cf187158bdp-12,
            -0x1.b5137f11a88e5p-14,
            0x1.056c6c57152b2p-15,
        };
        const z_squared = z * z;
        const z_fourth = z_squared * z_squared;
        const polynomial = ((coefficients[0] + z * coefficients[1]) +
            z_squared * (coefficients[2] + z * coefficients[3])) +
            z_fourth * ((coefficients[4] + z * coefficients[5]) +
                z_squared * (coefficients[6] + z * coefficients[7]));
        return @floatCast(root + (root * z) * polynomial);
    }

    const input: f64 = @floatCast(value);
    const input_squared = input * input;
    const root_argument = input_squared - 1.0;
    const root = @sqrt(root_argument);
    const sum = input + root;
    const sum_bits: u64 = @bitCast(sum);
    const mantissa = sum_bits & 0x000f_ffff_ffff_ffff;
    const j: usize = @intCast((mantissa + (1 << 44)) >> 45);
    const exponent: i64 = @as(i64, @intCast(sum_bits >> 52)) - 0x3ff;
    const normalized_bits = mantissa | 0x3ff0_0000_0000_0000;
    const normalized: f64 = @bitCast(normalized_bits);
    const z = normalized * acoshf_data.ix[j] - 1.0;
    const coefficients = [_]f64{
        0x1.0000000066947p+0,
        -0x1.00007f053d8cbp-1,
        0x1.555280111d914p-2,
    };
    const z_squared = z * z;
    var result = ((acoshf_data.lix[128] * @as(f64, @floatFromInt(exponent)) + acoshf_data.lix[j]) +
        z * coefficients[0]) + z_squared * (coefficients[1] + z * coefficients[2]);
    const result_bits: u64 = @bitCast(result);
    if (((result_bits + 259000) & 0x0fff_ffff) < 260000) {
        const accurate_coefficients = [_]f64{
            1.0,
            -0.5,
            0x1.55555555030bcp-2,
            -0x1.ffffffff2b4e5p-3,
            0x1.999b5076a42f2p-3,
            -0x1.55570c45a647dp-3,
        };
        const coefficient_0 = accurate_coefficients[0] + z * accurate_coefficients[1];
        const coefficient_2 = accurate_coefficients[2] + z * accurate_coefficients[3];
        const coefficient_4 = accurate_coefficients[4] + z * accurate_coefficients[5];
        const accurate_polynomial = coefficient_0 +
            z_squared * (coefficient_2 + z_squared * coefficient_4);
        const low_log_two = 0x1.7f7d1cf79abcap-20;
        const high_log_two = 0x1.62e4p-1;
        const high_part = high_log_two * @as(f64, @floatFromInt(exponent));
        const low_part = low_log_two * @as(f64, @floatFromInt(exponent));
        result = (z * accurate_polynomial + (low_part + acoshf_data.lix[j])) + high_part;
        if ((@as(u64, @bitCast(result)) & 0x0fff_ffff) == 0) {
            const correction = (z * accurate_polynomial + (low_part + acoshf_data.lix[j])) +
                (high_part - result);
            result += 64.0 * correction;
        }
    }
    return @floatCast(result);
}
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits_d: u64 = 0x7ff8_0000_0000_0000;
const acoshd_data = @import("scalar_acoshd_data.zig");

fn acoshdFastTwoSum(x: f64, y: f64, error_value: *f64) f64 {
    const sum = x + y;
    const z = sum - x;
    error_value.* = y - z;
    return sum;
}

fn acoshdTwoSum(a: f64, b: f64, error_value: *f64) f64 {
    const sum = a + b;
    const a_prime = sum - b;
    const b_prime = sum - a_prime;
    const delta_a = a - a_prime;
    const delta_b = b - b_prime;
    error_value.* = delta_a + delta_b;
    return sum;
}

fn acoshdAddDoubleDouble(x_high: f64, x_low: f64, c_high: f64, c_low: f64, low: *f64) f64 {
    const sum = acoshdTwoSum(x_high, c_high, low);
    low.* += x_low + c_low;
    return sum;
}

fn acoshdMulDoubleDoubleAccurate(
    x_high: f64,
    x_low: f64,
    c_high: f64,
    c_low: f64,
    low: *f64,
) f64 {
    const a_high_low = c_high * x_low;
    const a_low_high = c_low * x_high;
    const a_high_high = c_high * x_high;
    const a_high_error = @mulAdd(f64, c_high, x_high, -a_high_high);
    const middle = a_high_error + (a_low_high + a_high_low);
    return acoshdFastTwoSum(a_high_high, middle, low);
}

fn acoshdMulDoubleDouble(x_high: f64, x_low: f64, c: f64, low: *f64) f64 {
    const a_high_low = c * x_low;
    const a_high_high = c * x_high;
    const a_high_error = @mulAdd(f64, c, x_high, -a_high_high);
    const middle = a_high_error + a_high_low;
    const high = a_high_high + middle;
    low.* = (a_high_high - high) + middle;
    return high;
}

fn acoshdPolynomial(
    comptime coefficient_count: usize,
    x_high: f64,
    x_low: f64,
    coefficients: []const [2]f64,
    low: *f64,
) f64 {
    var index: usize = coefficient_count - 1;
    var coefficient_low: f64 = low.*;
    var coefficient_high = acoshdFastTwoSum(
        coefficients[index][0],
        coefficient_low,
        &coefficient_low,
    );
    coefficient_low += coefficients[index][1];
    while (index > 0) {
        index -= 1;
        coefficient_high = acoshdMulDoubleDoubleAccurate(
            x_high,
            x_low,
            coefficient_high,
            coefficient_low,
            &coefficient_low,
        );
        var term_low: f64 = undefined;
        const term_high = acoshdFastTwoSum(
            coefficients[index][0],
            coefficient_high,
            &term_low,
        );
        coefficient_high = term_high;
        coefficient_low += term_low + coefficients[index][1];
    }
    low.* = coefficient_low;
    return coefficient_high;
}

fn acoshdOne(x: f64, square_root_high: f64, square_root_low: f64) f64 {
    var low = x * (acoshd_data.acosh_one_cl[0] + x *
        (acoshd_data.acosh_one_cl[1] + x *
            (acoshd_data.acosh_one_cl[2] + x *
                (acoshd_data.acosh_one_cl[3] + x *
                    (acoshd_data.acosh_one_cl[4] + x * acoshd_data.acosh_one_cl[5])))));
    var high = acoshdPolynomial(10, x, 0.0, acoshd_data.acosh_one_ch[0..], &low);
    high = acoshdMulDoubleDouble(high, low, x, &low);
    var y1: f64 = undefined;
    var y0 = acoshdFastTwoSum(1.0, high, &y1);
    y1 += low;
    y0 = acoshdMulDoubleDoubleAccurate(y0, y1, square_root_high, square_root_low, &y1);
    return y0 + y1;
}

fn acoshdDatabase(x: f64, approximation: f64) f64 {
    var first: i64 = 0;
    var last: i64 = acoshd_data.database.len - 1;
    while (first <= last) {
        const middle: i64 = @divTrunc(first + last, 2);
        const candidate = acoshd_data.database[@intCast(middle)];
        if (candidate[0] < x) {
            first = middle + 1;
        } else if (candidate[0] > x) {
            last = middle - 1;
        } else {
            return candidate[1] + candidate[2];
        }
    }
    return approximation;
}

fn acoshdRefine(x: f64, approximation: f64) f64 {
    const input_bits: u64 = @bitCast(x);
    var high: f64 = undefined;
    var low: f64 = undefined;
    var exponent_adjustment: i64 = 0;
    if (input_bits < 0x4190_0000_0000_0000) {
        const x_squared_high = x * x;
        const x_squared_low = @mulAdd(f64, x, x, -x_squared_high);
        var root_argument_low: f64 = undefined;
        var root_argument_high = x_squared_high - 1.0;
        root_argument_high = acoshdFastTwoSum(
            root_argument_high,
            x_squared_low,
            &root_argument_low,
        );
        const root_high = @sqrt(root_argument_high);
        const root_low = (root_argument_low -
            @mulAdd(f64, root_high, root_high, -root_argument_high)) / (2.0 * root_high);
        high = acoshdFastTwoSum(x, root_high, &low);
        low += root_low;
        high = acoshdFastTwoSum(high, low, &low);
    } else if (input_bits < 0x4330_0000_0000_0000) {
        high = 2.0 * x;
        low = -0.5 / x;
    } else {
        high = x;
        low = 0.0;
        exponent_adjustment = 1;
    }

    var high_bits: u64 = @bitCast(high);
    const exponent: i64 = @as(i64, @intCast(high_bits >> 52)) - 0x3ff + exponent_adjustment;
    high_bits &= 0x000f_ffff_ffff_ffff;
    high_bits |= 0x3ff0_0000_0000_0000;
    const normalized: f64 = @bitCast(high_bits);
    const exponent_float: f64 = @floatFromInt(exponent);
    const approximation_bits: u64 = @bitCast(approximation - exponent_float + 0x1.00008p+0);
    const index: u64 = (approximation_bits - 0x3ff0_0000_0000_0000) >> 36;
    const index_1: usize = @intCast((index >> 12) & 0x1f);
    const index_2: usize = @intCast((index >> 8) & 0xf);
    const index_3: usize = @intCast((index >> 4) & 0xf);
    const index_4: usize = @intCast(index & 0xf);

    const log_two_low_low = -0x1.9ff0342542fc3p-91;
    const log_two_low = 0x1.ef35793c768p-46;
    const log_two_high = 0x1.62e42fefa38p-2;
    const log_two_low_low_scaled = log_two_low_low * exponent_float;
    const log_two_low_scaled = log_two_low * exponent_float;
    const log_two_high_scaled = log_two_high * exponent_float;
    var logarithm: [3]f64 = undefined;
    logarithm[0] = acoshd_data.refine_ll[0][index_1][0] +
        acoshd_data.refine_ll[1][index_2][0] +
        (acoshd_data.refine_ll[2][index_3][0] + acoshd_data.refine_ll[3][index_4][0]);
    logarithm[1] = acoshd_data.refine_ll[0][index_1][1] +
        acoshd_data.refine_ll[1][index_2][1] +
        (acoshd_data.refine_ll[2][index_3][1] + acoshd_data.refine_ll[3][index_4][1]);
    logarithm[2] = acoshd_data.refine_ll[0][index_1][2] +
        acoshd_data.refine_ll[1][index_2][2] +
        (acoshd_data.refine_ll[2][index_3][2] + acoshd_data.refine_ll[3][index_4][2]);
    logarithm[0] += log_two_high_scaled;

    const t12 = acoshd_data.refine_t1[index_1] * acoshd_data.refine_t2[index_2];
    const t34 = acoshd_data.refine_t3[index_3] * acoshd_data.refine_t4[index_4];
    const table_high = t12 * t34;
    const table_low = @mulAdd(f64, t12, t34, -table_high);
    const product_high = table_high * normalized;
    const product_low = @mulAdd(f64, table_high, normalized, -product_high);
    const product_low_high = table_low * normalized;
    const product_low_low = @mulAdd(f64, table_low, normalized, -product_low_high);
    var reduced_low: f64 = undefined;
    var reduced_high = acoshdFastTwoSum(product_high - 1.0, product_low, &reduced_low);
    if (low != 0.0) {
        var low_bits: u64 = @bitCast(low);
        const exponent_shift = @as(u64, @intCast(exponent)) << 52;
        low_bits -%= exponent_shift;
        reduced_low += table_high * @as(f64, @bitCast(low_bits));
    }
    reduced_high = acoshdAddDoubleDouble(
        reduced_high,
        reduced_low,
        product_low_high,
        product_low_low,
        &reduced_low,
    );

    var series_low = reduced_high * (acoshd_data.refine_cl[0] +
        reduced_high * (acoshd_data.refine_cl[1] + reduced_high * acoshd_data.refine_cl[2]));
    var series_high = acoshdPolynomial(3, reduced_high, reduced_low, acoshd_data.refine_ch[0..], &series_low);
    series_high = acoshdMulDoubleDoubleAccurate(
        reduced_high,
        reduced_low,
        series_high,
        series_low,
        &series_low,
    );
    series_high = acoshdAddDoubleDouble(
        series_high,
        series_low,
        log_two_low_scaled,
        log_two_low_low_scaled,
        &series_low,
    );
    series_high = acoshdAddDoubleDouble(
        series_high,
        series_low,
        logarithm[1],
        logarithm[2],
        &series_low,
    );

    var second_low: f64 = undefined;
    var result_high = acoshdFastTwoSum(logarithm[0], series_high, &second_low);
    var result_middle = acoshdFastTwoSum(second_low, series_low, &second_low);
    result_high *= 2.0;
    result_middle *= 2.0;
    second_low *= 2.0;
    var result_middle_bits: u64 = @bitCast(result_middle);
    if ((result_middle_bits & 0x000f_ffff_ffff_ffff) == 0) {
        const second_low_bits: u64 = @bitCast(second_low);
        if (((second_low_bits ^ result_middle_bits) >> 63) != 0) {
            result_middle_bits -%= 1;
        } else {
            result_middle_bits +%= 1;
        }
        result_middle = @bitCast(result_middle_bits);
    }
    const result_high_bits: u64 = @bitCast(result_high);
    const error_region = (result_middle_bits + 7) & 0x000f_ffff_ffff_ffff;
    const exponent_difference = ((result_high_bits >> 52) & 0x7ff) -
        ((result_middle_bits >> 52) & 0x7ff);
    const result = result_high + result_middle;
    if (exponent_difference > 102 or error_region < 15) {
        return acoshdDatabase(x, result);
    }
    return result;
}

/// Correctly rounded binary64 inverse hyperbolic cosine under caller round-to-nearest-even.
pub fn acoshd(value: f64) f64 {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    if (bits >= infinity_bits_d) {
        const magnitude = bits & magnitude_mask_d;
        if (magnitude > infinity_bits_d) return @bitCast(bits | quiet_nan_bit_d);
        if (bits == infinity_bits_d) return value;
        return @bitCast(canonical_nan_bits_d);
    }
    if (bits <= 0x3ff0_0000_0000_0000) {
        if (bits == 0x3ff0_0000_0000_0000) return 0.0;
        return @bitCast(canonical_nan_bits_d);
    }

    var approximation: f64 = undefined;
    var error_bound: f64 = undefined;
    var offset: i64 = 0x3fe;
    var normalized_bits = bits;
    if (bits < 0x3ff1_e83e_425a_ee63) {
        const z = value - 1.0;
        const inverse = -0.25 / z;
        const doubled = 2.0 * z;
        const square_root_high = @sqrt(doubled);
        const square_root_low = @mulAdd(f64, square_root_high, square_root_high, -doubled) *
            (square_root_high * inverse);
        const coefficients = [_]f64{
            -0x1.5555555555555p-4,
            0x1.3333333332f95p-6,
            -0x1.6db6db6d5534cp-8,
            0x1.f1c71c1e04356p-10,
            -0x1.6e8b8e3e40d58p-11,
            0x1.1c4ba825ac4fep-12,
            -0x1.c9045534e6d9ep-14,
            0x1.71fedae26a76bp-15,
            -0x1.f1f4f8cc65342p-17,
        };
        const z_squared = z * z;
        const z_fourth = z_squared * z_squared;
        const coefficient_1 = @mulAdd(f64, z, coefficients[2], coefficients[1]);
        const coefficient_3 = @mulAdd(f64, z, coefficients[4], coefficients[3]);
        const coefficient_5 = @mulAdd(f64, z, coefficients[6], coefficients[5]);
        const coefficient_7 = @mulAdd(f64, z, coefficients[8], coefficients[7]);
        const first = @mulAdd(f64, z_squared, coefficient_3, coefficient_1);
        const second = @mulAdd(f64, z_squared, coefficient_7, coefficient_5);
        const polynomial = @mulAdd(
            f64,
            z,
            @mulAdd(f64, z_fourth, second, first),
            coefficients[0],
        );
        const delta = @mulAdd(f64, square_root_high * z, polynomial, square_root_low);
        const epsilon = delta * 0x1.00p-50 - 0x1p-104 * square_root_high;
        const lower = square_root_high + (delta - epsilon);
        const upper = square_root_high + (delta + epsilon);
        if (lower == upper) return lower;
        return acoshdOne(z, square_root_high, square_root_low);
    } else if (bits < 0x405b_f000_0000_0000) {
        offset = 0x3ff;
        const squared_high = value * value;
        const squared_low = @mulAdd(f64, value, value, -squared_high);
        const root_argument_low = squared_low;
        const root_argument_high = squared_high - 1.0;
        const square_root_high = @sqrt(root_argument_high);
        const inverse_root_argument = 0.5 / root_argument_high;
        const square_root_low = (root_argument_low -
            @mulAdd(f64, square_root_high, square_root_high, -root_argument_high)) *
            (square_root_high * inverse_root_argument);
        var sum_low: f64 = undefined;
        const sum_high = acoshdFastTwoSum(value, square_root_high, &sum_low);
        sum_low += square_root_low;
        normalized_bits = @bitCast(sum_high);
        approximation = sum_low / sum_high;
        error_bound = 0x1.81p-63;
    } else if (bits < 0x4087_1000_0000_0000) {
        const coefficients = [_]f64{
            0x1.5c4b6148816e2p-66,
            -0x1.000000000005cp-2,
            -0x1.7fffffebf3e6cp-4,
            -0x1.aab6691f2bae7p-5,
        };
        const inverse_squared = 1.0 / (value * value);
        approximation = @mulAdd(
            f64,
            inverse_squared,
            @mulAdd(
                f64,
                inverse_squared,
                @mulAdd(f64, inverse_squared, coefficients[3], coefficients[2]),
                coefficients[1],
            ),
            coefficients[0],
        );
        error_bound = 0x1.c3p-63;
    } else if (bits < 0x40e0_1000_0000_0000) {
        const coefficients = [_]f64{
            -0x1.7f77c8429c6c6p-67,
            -0x1.ffffffffff214p-3,
            -0x1.8000268641bfep-4,
        };
        const inverse_squared = 1.0 / (value * value);
        approximation = @mulAdd(
            f64,
            inverse_squared,
            @mulAdd(f64, inverse_squared, coefficients[2], coefficients[1]),
            coefficients[0],
        );
        error_bound = 0x1.9ap-63;
    } else if (bits < 0x41ea_0000_0000_0000) {
        const coefficients = [_]f64{
            0x1.7a0ed2effdd1p-67,
            -0x1.000000017d048p-2,
        };
        const inverse_squared = 1.0 / (value * value);
        approximation = @mulAdd(f64, inverse_squared, coefficients[1], coefficients[0]);
        error_bound = 0x1.99p-63;
    } else {
        approximation = 0.0;
        error_bound = 0x1.b2p-63;
    }

    const exponent = @as(i64, @intCast(normalized_bits >> 52)) - offset;
    normalized_bits &= 0x000f_ffff_ffff_ffff;
    const normalized: f64 = @bitCast(normalized_bits | 0x3ff0_0000_0000_0000);
    const exponent_float: f64 = @floatFromInt(exponent);
    const significand_index: usize = @intCast(normalized_bits >> 47);
    const low_significand: i64 = @intCast(normalized_bits & 0x0000_7fff_ffff_ffff);
    const interpolation_index_value: i64 = @as(i64, @intCast(normalized_bits)) +
        (@as(i64, acoshd_data.b_c0[significand_index]) << 33) +
        (@as(i64, acoshd_data.b_c1[significand_index]) * (low_significand >> 16));
    const interpolation_index: u64 = @intCast(@divTrunc(interpolation_index_value, 1 << 42));
    const interpolation_index_1: usize = @intCast(interpolation_index >> 5);
    const interpolation_index_2: usize = @intCast(interpolation_index & 0x1f);
    const reciprocal = acoshd_data.r1[interpolation_index_1] *
        acoshd_data.r2[interpolation_index_2];
    const delta = @mulAdd(f64, reciprocal, normalized, -1.0);
    const delta_squared = delta * delta;
    const first = @mulAdd(f64, delta, acoshd_data.log_poly[1], acoshd_data.log_poly[0]);
    const second = @mulAdd(f64, delta, acoshd_data.log_poly[3], acoshd_data.log_poly[2]);
    const tail = @mulAdd(f64, delta_squared, acoshd_data.log_poly[4], second);
    const log_correction = delta_squared * @mulAdd(f64, delta_squared, tail, first);
    const log_two_high = 0x1.62e42fefa38p-1;
    const log_two_low = 0x1.ef35793c7673p-45;
    const logarithm_high = @mulAdd(
        f64,
        log_two_high,
        exponent_float,
        acoshd_data.l1[interpolation_index_1][1] +
            acoshd_data.l2[interpolation_index_2][1],
    );
    const logarithm_low = log_two_low * exponent_float +
        (acoshd_data.l1[interpolation_index_1][0] + acoshd_data.l2[interpolation_index_2][0]);
    const logarithm_correction = log_correction + logarithm_low;
    const logarithm_delta = delta + (approximation + logarithm_correction);
    const lower = logarithm_high + (logarithm_delta - error_bound);
    const upper = logarithm_high + (logarithm_delta + error_bound);
    if (lower == upper) return lower;
    return acoshdRefine(value, 0x1.71547652b82fep+0 * lower);
}
