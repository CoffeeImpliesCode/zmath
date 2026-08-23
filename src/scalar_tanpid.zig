const std = @import("std");

// Correctly rounded half-revolution tangent for binary64 values.
//
// Copyright (c) 2023-2025 Alexei Sibidanov.
//
// This code is derived from CORE-MATH src/binary64/tanpi/tanpi.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT License.
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

const magnitude_mask: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits: u64 = 0x7ff8_0000_0000_0000;

fn copySign(value: f64, sign_source: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    const sign_bits: u64 = @bitCast(sign_source);
    return @bitCast((value_bits & magnitude_mask) | (sign_bits & 0x8000_0000_0000_0000));
}

fn fastTwoSum(left: f64, right: f64, residual: *f64) f64 {
    const sum = left + right;
    const correction = sum - left;
    residual.* = right - correction;
    return sum;
}

fn fastTwoSub(left: f64, right: f64, residual: *f64) f64 {
    const difference = left - right;
    const correction = left - difference;
    residual.* = correction - right;
    return difference;
}

fn mulDdAcc(high: f64, low: f64, coefficient_high: f64, coefficient_low: f64, low_out: *f64) f64 {
    const high_low = coefficient_high * low;
    const low_high = coefficient_low * high;
    const high_high = coefficient_high * high;
    var high_error = @mulAdd(f64, coefficient_high, high, -high_high);
    high_error += low_high + high_low;
    return fastTwoSum(high_high, high_error, low_out);
}

fn mulDdd(high: f64, low: f64, coefficient: f64, low_out: *f64) f64 {
    const high_low = coefficient * low;
    const high_high = coefficient * high;
    var high_error = @mulAdd(f64, coefficient, high, -high_high);
    high_error += high_low;
    const sum = high_high + high_error;
    low_out.* = (high_high - sum) + high_error;
    return sum;
}

fn fastSumDdd(left_high: f64, left_low: f64, right_high: f64, low_out: *f64) f64 {
    var temporary_low: f64 = undefined;
    const sum_high = fastTwoSum(left_high, right_high, &temporary_low);
    low_out.* = left_low + temporary_low;
    return sum_high;
}

fn polyDd(value_high: f64, value_low: f64, coefficients: []const [2]f64, low_out: *f64) f64 {
    var index = coefficients.len - 1;
    var high = coefficients[index][0] + low_out.*;
    var low = ((coefficients[index][0] - high) + low_out.*) + coefficients[index][1];
    while (index > 0) {
        index -= 1;
        high = mulDdAcc(value_high, value_low, high, low, &low);
        const sum = high + coefficients[index][0];
        const correction = (coefficients[index][0] - sum) + high;
        high = sum;
        low += correction + coefficients[index][1];
    }
    low_out.* = low;
    return high;
}

const database = [_][3]f64{
    .{ 0x1.b44453e2404e7p-55, 0x1.56a4aa740a5a7p-53, 0x1p-107 },
    .{ 0x1.45435c15060dcp-54, 0x1.feec101084f24p-53, -0x1p-107 },
    .{ 0x1.ac84c88f979a2p-54, 0x1.508ecb38f52f9p-52, 0x1p-106 },
    .{ 0x1.b44453e2404e7p-54, 0x1.56a4aa740a5a7p-52, 0x1p-106 },
    .{ 0x1.5f095755008fp-52, 0x1.13b412ea2182dp-50, 0x1p-104 },
    .{ 0x1.736d8818a1ef3p-48, 0x1.23b8007cd2448p-46, -0x1p-100 },
    .{ 0x1.c7c00d31634a6p-48, 0x1.65f203d6bb18ep-46, 0x1p-100 },
    .{ 0x1.df6d9f7748faep-45, 0x1.788ac313fc061p-43, 0x1p-97 },
    .{ 0x1.292b56d105932p-44, 0x1.d2cadcfe9f85fp-43, -0x1p-97 },
    .{ 0x1.45a65501d94d3p-44, 0x1.ff878736a1bcdp-43, 0x1p-97 },
    .{ 0x1.e4a84b3095d13p-43, 0x1.7ca626cfcf41p-41, -0x1p-95 },
    .{ 0x1.9eb0b7b0519dap-42, 0x1.45b26ebd51c87p-40, 0x1p-94 },
    .{ 0x1.1d0798b3403ebp-41, 0x1.bfb93b8757f7cp-40, -0x1p-94 },
    .{ 0x1.05d3aa746e163p-38, 0x1.9b46cfb0f10cap-37, 0x1p-91 },
    .{ 0x1.d3f710390db87p-36, 0x1.6f89f6df8258dp-34, 0x1p-88 },
    .{ 0x1.9e2371e233d1bp-34, 0x1.45437a2ebe657p-32, 0x1p-86 },
    .{ 0x1.386de1f01225ap-32, 0x1.eac33f52580c1p-31, 0x1p-85 },
    .{ 0x1.520cc578bf7ccp-30, 0x1.0980f68058c7fp-28, -0x1p-82 },
    .{ 0x1.a5bdf2b61615bp-23, 0x1.4b3c41e10a1d4p-21, 0x1p-75 },
    .{ 0x1.ed835e8ecc245p-23, 0x1.839ab56b14ffdp-21, 0x1p-75 },
    .{ 0x1.5b79e364c2d17p-22, 0x1.10e838704c63dp-20, 0x1p-74 },
    .{ 0x1.a96d0004e5b8bp-21, 0x1.4e20edbb1c21bp-19, 0x1p-73 },
    .{ 0x1.7b544ccbafe6ap-18, 0x1.29ecae32f9ca3p-16, 0x1p-70 },
    .{ 0x1.a2f422a07a3d4p-16, 0x1.490ba1569549fp-14, 0x1p-68 },
    .{ 0x1.ee10ecdb9b092p-15, 0x1.8409e333ea2e8p-13, -0x1p-67 },
    .{ 0x1.ee3f52b4af0dcp-15, 0x1.842e540b07ac1p-13, 0x1p-67 },
    .{ 0x1.1c2dc0a9340efp-13, 0x1.be630d1090fe6p-12, 0x1p-66 },
    .{ 0x1.75fdc341dc7a1p-13, 0x1.25bb6997b4387p-11, -0x1p-65 },
    .{ 0x1.a572fa6b3d117p-12, 0x1.4b016bd111ed5p-10, -0x1p-64 },
    .{ 0x1.a537a6de092a1p-11, 0x1.4ad2f60d5240bp-9, 0x1p-63 },
    .{ 0x1.db74ba716abb6p-11, 0x1.756c5a7eec959p-9, 0x1p-63 },
    .{ 0x1.5d89dc03193cdp-10, 0x1.12874c9307181p-8, -0x1p-62 },
    .{ 0x1.54ae0fa9f24b9p-7, 0x1.0baa21da88501p-5, -0x1p-59 },
    .{ 0x1.8582f89424fcep-7, 0x1.321061924ca03p-5, 0x1p-59 },
    .{ 0x1.14fc0e0c2dacep-6, 0x1.b37ef5831a537p-5, 0x1p-59 },
    .{ 0x1.309ce21202ccap-5, 0x1.e0ac5527e7f0dp-4, 0x1p-58 },
    .{ 0x1.bb6715a0b1f81p-5, 0x1.5fa4b5e7ce93fp-3, 0x1p-57 },
    .{ 0x1.1e38ee7c71b98p-4, 0x1.c8f62be40c3a1p-3, -0x1p-57 },
    .{ 0x1.983c8c74ade5dp-4, 0x1.4b88c3e5b04d0p-2, 0x1p-56 },
    .{ 0x1.a334561b8192cp-4, 0x1.5513de069862bp-2, 0x1p-56 },
    .{ 0x1.a33741abea942p-4, 0x1.55166a4d7d836p-2, 0x1p-56 },
    .{ 0x1.c176d0b4666d3p-4, 0x1.6fb0d2df2ae41p-2, 0x1p-56 },
    .{ 0x1.ca50e673ede6bp-3, 0x1.b1ed1cc270a2dp-1, 0x1p-55 },
    .{ 0x1.df3c63292f1f3p-3, 0x1.cef531bfb44bfp-1, -0x1p-55 },
    .{ 0x1.f551afda4c896p-3, 0x1.ef7dd0e5c6350p-1, 0x1p-55 },
    .{ 0x1.8b9be647e02a9p-2, 0x1.5717df9abccd3p+1, 0x1p-53 },
    .{ 0x1.b648e5b0804dap-2, 0x1.16265f83e8121p+2, -0x1p-52 },
    .{ 0x1.dd9d1ecffc9a0p-2, 0x1.2e342ca2440eep+3, 0x1p-51 },
    .{ 0x1.ee846d3680301p-2, 0x1.2a04fa64d35ffp+4, 0x1p-50 },
};

fn tanpiDatabase(value: f64, result: f64) f64 {
    const input_bits: u64 = @bitCast(value);
    const absolute_bits = input_bits & magnitude_mask;
    var period: u64 = 0;
    var absolute_value = copySign(value, 1.0);
    const exponent = @as(i32, @intCast(absolute_bits >> 52)) - 1022;
    if (exponent >= 0) {
        const shift: u6 = @intCast(12 + exponent);
        const truncated_bits = absolute_bits & ~(@as(u64, std.math.maxInt(u64)) >> shift);
        absolute_value -= @as(f64, @bitCast(truncated_bits));
        period = ((absolute_bits | (@as(u64, 1) << 52)) >> @intCast(52 - exponent)) & 1;
        if (period != 0) absolute_value = 0.5 - absolute_value;
    }
    const sign: f64 = if ((period ^ (input_bits >> 63)) != 0) -1.0 else 1.0;
    var lower: i32 = 0;
    var upper: i32 = @intCast(database.len - 1);
    while (lower <= upper) {
        const middle: i32 = @divTrunc(lower + upper, 2);
        const middle_index: usize = @intCast(middle);
        if (database[middle_index][0] < absolute_value) {
            lower = middle + 1;
        } else if (database[middle_index][0] == absolute_value) {
            return sign * database[middle_index][1] + sign * database[middle_index][2];
        } else {
            upper = middle - 1;
        }
    }
    return result;
}

const correction_table = [_][2]f64{
    .{ 0x0p+0, 0x0p+0 },
    .{ 0x1.ada13ceebab9dp-64, 0x1.927278a3b1162p-5 },
    .{ 0x1.87d054f689d7ap-58, 0x1.936bb8c5b2da2p-4 },
    .{ 0x1.2cfb5a746f62cp-58, 0x1.2fcac73a6064p-3 },
    .{ 0x1.ef5d367441946p-61, 0x1.975f5e0553158p-3 },
    .{ 0x1.5a6d6c3c8b6a7p-57, 0x1.007fa758626aep-2 },
    .{ -0x1.d704d1bfdb6e8p-57, 0x1.36a08355c63dcp-2 },
    .{ 0x1.78e945dc3913cp-57, 0x1.6e649f7d78649p-2 },
    .{ 0x1.08b2fb1366ea9p-56, 0x1.a827999fcef32p-2 },
    .{ 0x1.41522f15e53c5p-56, 0x1.e450e0d273e7ap-2 },
    .{ -0x1.e564bcd1e635dp-55, 0x1.11ab7190834ecp-1 },
    .{ -0x1.10b4421e6a4f8p-55, 0x1.32e1889047ffdp-1 },
    .{ 0x1.7a8c52172b675p-55, 0x1.561b82ab7f99p-1 },
    .{ -0x1.aa7538e44e996p-55, 0x1.7bb99ed2990cfp-1 },
    .{ -0x1.a6db80fe796fep-56, 0x1.a43002ae4285p-1 },
    .{ 0x1.8dcad85e60fbep-56, 0x1.d00cbc7384d2ep-1 },
    .{ 0x0p+0, 0x1p+0 },
    .{ 0x1.946cc0b66979fp-54, 0x1.1a73d55278c4bp+0 },
    .{ -0x1.87e5ad9735569p-56, 0x1.37efd8d87607ep+0 },
    .{ 0x1.6a085e3bc3af6p-55, 0x1.592d11142fa55p+0 },
    .{ 0x1.419fa6954928fp-54, 0x1.7f218e25a7461p+0 },
    .{ -0x1.b6fb77846d763p-56, 0x1.ab1c35d8a74eap+0 },
    .{ 0x1.0fb3e75c7098ep-55, 0x1.def13b73c1406p+0 },
    .{ 0x1.7ce6cb463c972p-54, 0x1.0ea21d716fbf7p+1 },
    .{ 0x1.21165f626cdd5p-53, 0x1.3504f333f9de6p+1 },
    .{ 0x1.aca947bfb1dccp-56, 0x1.65bc6cc825147p+1 },
    .{ -0x1.b7a14d0d691bep-54, 0x1.a5f59e90600ddp+1 },
    .{ -0x1.89fcd637fbf3bp-55, 0x1.ff01305ecd8dcp+1 },
    .{ 0x1.10706fed37f0ep-55, 0x1.41bfee2424771p+2 },
    .{ -0x1.ae397239c5a0dp-53, 0x1.af73f4ca3310fp+2 },
    .{ 0x1.351daea79411dp-51, 0x1.44e6c595afdccp+3 },
    .{ -0x1.b6e6b4de0cd24p-56, 0x1.45affed201b55p+4 },
};

const correction_coefficients = [_]f64{
    0x1.4abbce625be51p-207,
    0x1.466bc6776a9b1p-345,
    0x1.45fff6eb26045p-483,
    0x1.4627663861052p-621,
};

const ch_large = [_][2]f64{
    .{ 0x1.921fb54442d18p-6, 0x1.1a62633145c07p-60 },
    .{ 0x1.4abbce625be53p-18, -0x1.05511c6842515p-72 },
    .{ 0x1.466bc6775aae2p-30, -0x1.6dc0d93fb2ecep-85 },
    .{ 0x1.45fff9b48e95ep-42, 0x1.b226b250d2cc3p-98 },
    .{ 0x1.45f472e3af011p-54, 0x1.90612a0755449p-110 },
};

const cl_large = [_]f64{
    0x1.45f32f25dab7fp-66,
    0x1.45f3030c82af4p-78,
    0x1.464490600a978p-90,
};

const ch_small = [_][2]f64{
    .{ 0x1.4abbce625be53p+3, -0x1.05511c68476a8p-51 },
    .{ 0x1.466bc6775aae2p+5, -0x1.6dc0cbddc0e69p-50 },
    .{ 0x1.45fff9b48e95ep+7, 0x1.a5047910ae0efp-49 },
};

const cl_small = [_]f64{
    0x1.45f472e3aed7dp+9,
    0x1.45f33be0e9598p+11,
};

/// Return the correctly rounded tangent of pi times `value`.
///
/// The caller must keep IEEE round-to-nearest, ties-to-even active.
pub fn tanpid(value: f64) f64 {
    @setFloatMode(.strict);
    const input_bits: u64 = @bitCast(value);
    const absolute_bits = input_bits & magnitude_mask;
    var high: f64 = undefined;
    var low: f64 = undefined;
    var result: f64 = undefined;

    if (absolute_bits >= (0x3f3 << 52)) {
        if (absolute_bits >= (0x42d << 52)) {
            if (absolute_bits >= infinity_bits) {
                if (absolute_bits > infinity_bits) return @bitCast(input_bits | quiet_nan_bit);
                return @bitCast(canonical_nan_bits);
            }
            const exponent: i32 = @intCast(absolute_bits >> 52);
            const shift_signed = exponent - 1069;
            if (shift_signed > 6) return @bitCast(input_bits & 0x8000_0000_0000_0000);
            const shift: u6 = @intCast(shift_signed);
            const sign_bits: u64 = if ((input_bits >> 63) != 0) std.math.maxInt(u64) else 0;
            const quotient_bits = (((absolute_bits ^ sign_bits) -% sign_bits) << shift) & 127;
            if ((quotient_bits & 31) == 0) {
                const quadrant: u64 = quotient_bits >> 5;
                if ((quadrant & 1) != 0) {
                    if ((quadrant & 2) != 0) return @bitCast(@as(u64, 0xfff0_0000_0000_0000));
                    return @bitCast(@as(u64, 0x7ff0_0000_0000_0000));
                }
                if (((quadrant ^ sign_bits) & 2) != 0) {
                    return @bitCast(@as(u64, 0x8000_0000_0000_0000));
                }
                return 0.0;
            }
            const table_bits = quotient_bits & 31;
            if ((quotient_bits & 32) != 0) {
                const reflected_index: usize = @intCast(32 - table_bits);
                return -(correction_table[reflected_index][1] +
                    correction_table[reflected_index][0]);
            }
            const table_index: usize = @intCast(table_bits);
            return correction_table[table_index][1] + correction_table[table_index][0];
        }

        const exponent: i32 = @intCast(absolute_bits >> 52);
        const shift: u6 = @intCast(1068 - exponent);
        const shift_1: u6 = @intCast(exponent - 1011);
        const mantissa = (absolute_bits & (@as(u64, std.math.maxInt(u64)) >> 12)) |
            (@as(u64, 1) << 52);
        const mantissa_shifted = mantissa << shift_1;
        const half_bit: i64 = @as(i64, @bitCast(mantissa_shifted)) >> 63;
        var quotient = ((mantissa ^ @as(u64, @bitCast(half_bit))) >> shift) & 63;
        quotient = (quotient + 1) >> 1;
        const sign_bits: i64 = @as(i64, @bitCast(input_bits)) >> 63;
        const half_sign: i64 = half_bit ^ sign_bits;
        const signed_mantissa: i64 = @bitCast(
            (mantissa ^ @as(u64, @bitCast(sign_bits))) -% @as(u64, @bitCast(sign_bits)),
        );
        const low_bits = @as(u64, @bitCast(signed_mantissa)) << @intCast(exponent - 1005);
        var reduced: f64 = @floatFromInt(@as(i64, @bitCast(low_bits)));

        if ((low_bits << 1) == 0) {
            if (low_bits == 0) {
                if ((quotient & 31) == 0) {
                    const quadrant: i64 = signed_mantissa >> @intCast(shift + 6);
                    if ((quadrant & 1) != 0) {
                        if ((quadrant & 2) != 0) return @bitCast(@as(u64, 0xfff0_0000_0000_0000));
                        return @bitCast(@as(u64, 0x7ff0_0000_0000_0000));
                    }
                    if (((@as(u64, @bitCast(quadrant)) ^ @as(u64, @bitCast(sign_bits))) & 2) != 0) {
                        return @bitCast(@as(u64, 0x8000_0000_0000_0000));
                    }
                    return 0.0;
                }
                const quarter = (mantissa << shift_1) >> 58;
                if (quarter == 0x10) return copySign(1.0, value);
                if (quarter == 0x30) return -copySign(1.0, value);
            }
            reduced = copySign(1.0, value) * reduced;
        }

        const reduced_square = reduced * reduced;
        const reduced_fourth = reduced_square * reduced_square;
        const reduced_cube = reduced * reduced_square;
        const polynomial_tail_low = correction_coefficients[0] +
            reduced_square * correction_coefficients[1];
        const polynomial_tail_high = correction_coefficients[2] +
            reduced_square * correction_coefficients[3];
        const polynomial_tail = polynomial_tail_low + reduced_fourth * polynomial_tail_high;
        const polynomial = reduced_cube * polynomial_tail;
        var bound_error = reduced_cube * 0x1p-256 + copySign(0x1p-104, reduced);
        high = mulDdd(0x1.921fb54442d18p-69, 0x1.1a62633145c07p-123, reduced, &low);
        high = fastSumDdd(high, low, polynomial, &low);
        if (quotient == 32) {
            const inverse = -1.0 / high;
            low = (@mulAdd(f64, inverse, high, 1.0) + low * inverse) * inverse;
            high = inverse;
        } else {
            var coefficient_low = correction_table[quotient][0];
            var coefficient_high = correction_table[quotient][1];
            const sign: f64 = if (half_sign + 1 == 0) -1.0 else 1.0;
            coefficient_high *= sign;
            coefficient_low *= sign;
            var product_low: f64 = undefined;
            var product_high = mulDdAcc(high, low, coefficient_high, coefficient_low, &product_low);
            var difference_low: f64 = undefined;
            product_high = fastTwoSub(1.0, product_high, &difference_low);
            product_low = difference_low - product_low;
            var sum_low: f64 = undefined;
            coefficient_high = fastTwoSum(coefficient_high, high, &sum_low);
            coefficient_low += sum_low + low;
            const inverse = 1.0 / product_high;
            high = coefficient_high * inverse;
            low = @mulAdd(f64, coefficient_high, inverse, -high) +
                (coefficient_low + coefficient_high * (@mulAdd(
                    f64,
                    -product_high,
                    inverse,
                    1.0,
                ) - product_low * inverse)) * inverse;
        }
        bound_error += bound_error * (high * high);
        const lower_bound = high + (low - bound_error);
        const upper_bound = high + (low + bound_error);
        if (lower_bound == upper_bound) return lower_bound;
        reduced *= 0x1p-63;
        const reduced_square_error = @mulAdd(f64, reduced, reduced, -(reduced * reduced));
        low = reduced * reduced *
            (cl_large[0] + reduced * reduced * (cl_large[1] + reduced * reduced * cl_large[2]));
        high = polyDd(reduced * reduced, reduced_square_error, ch_large[0..], &low);
        high = mulDdd(high, low, reduced, &low);
        if (quotient == 32) {
            const inverse = -1.0 / high;
            low = (@mulAdd(f64, inverse, high, 1.0) + low * inverse) * inverse;
            high = inverse;
        } else {
            var coefficient_low = correction_table[quotient][0];
            var coefficient_high = correction_table[quotient][1];
            const sign: f64 = if (half_sign + 1 == 0) -1.0 else 1.0;
            coefficient_high *= sign;
            coefficient_low *= sign;
            var product_low: f64 = undefined;
            var product_high = mulDdAcc(high, low, coefficient_high, coefficient_low, &product_low);
            var difference_low: f64 = undefined;
            product_high = fastTwoSub(1.0, product_high, &difference_low);
            product_low = difference_low - product_low;
            var sum_low: f64 = undefined;
            coefficient_high = fastTwoSum(coefficient_high, high, &sum_low);
            coefficient_low += sum_low + low;
            const inverse = 1.0 / product_high;
            high = coefficient_high * inverse;
            low = @mulAdd(f64, coefficient_high, inverse, -high) +
                (coefficient_low + coefficient_high * (@mulAdd(
                    f64,
                    -product_high,
                    inverse,
                    1.0,
                ) - product_low * inverse)) * inverse;
        }
        high = fastTwoSum(high, low, &low);
        result = high;
    } else {
        if (absolute_bits == 0) return value;
        const pi_high = 0x1.921fb54442d18p+1;
        const pi_low = 0x1.1a62633145c07p-53;
        if (absolute_bits < (0x3ca << 52)) {
            if (absolute_bits < (0x6b << 52)) {
                const exponent: i32 = @intCast(absolute_bits >> 52);
                const scale_bits: u64 = @intCast(@as(i64, 2045 - exponent) << 52);
                const inverse_scale_bits: u64 = @intCast(@as(i64, 1 + exponent) << 52);
                const scaled = value * @as(f64, @bitCast(scale_bits));
                high = mulDdd(pi_high, pi_low, scaled, &low);
                result = high * @as(f64, @bitCast(inverse_scale_bits));
                if (@abs(result) < 0x1p-1022) {
                    const smallest_normal = copySign(0x1p-1022, value);
                    var bounded = (smallest_normal + result) * @as(f64, @bitCast(scale_bits));
                    const high_scaled = result * @as(f64, @bitCast(scale_bits));
                    low += high - high_scaled;
                    bounded += low;
                    return bounded * @as(f64, @bitCast(inverse_scale_bits)) - smallest_normal;
                }
            } else {
                high = mulDdd(pi_high, pi_low, value, &low);
                result = high;
            }
        } else {
            const small_coefficients = [_]f64{
                0x1.4abbce625be53p+3,
                0x1.466bc6775aae2p+5,
                0x1.46000158496c2p+7,
            };
            const value_square = value * value;
            const value_cube = value * value_square;
            const polynomial = value_cube * (small_coefficients[0] + value_square *
                (small_coefficients[1] + value_square * small_coefficients[2]));
            var product_low: f64 = undefined;
            const product_high = mulDdd(pi_high, pi_low, value, &product_low);
            high = fastSumDdd(product_high, product_low, polynomial, &low);
            const bound_error = value * (value_square * 0x1.1p-47 + 0x1p-101);
            const lower_bound = high + (low - bound_error);
            const upper_bound = high + (low + bound_error);
            if (lower_bound == upper_bound) return lower_bound;
            const small_ch = [_][2]f64{
                .{ 0x1.4abbce625be53p+3, -0x1.05511c68476a8p-51 },
                .{ 0x1.466bc6775aae2p+5, -0x1.6dc0cbddc0e69p-50 },
                .{ 0x1.45fff9b48e95ep+7, 0x1.a5047910ae0efp-49 },
            };
            const small_cl = [_]f64{ 0x1.45f472e3aed7dp+9, 0x1.45f33be0e9598p+11 };
            const value_square_error = @mulAdd(f64, value, value, -value_square);
            const value_cube_error = @mulAdd(
                f64,
                value_square,
                value,
                -value_cube,
            ) + value_square_error * value;
            low = value_square * (small_cl[0] + value_square * small_cl[1]);
            high = polyDd(value_square, value_square_error, small_ch[0..], &low);
            high = mulDdAcc(value_cube, value_cube_error, high, low, &low);
            var sum_error: f64 = undefined;
            high = fastTwoSum(product_high, high, &sum_error);
            low += product_low + sum_error;
            high = fastTwoSum(high, low, &low);
            result = high;
        }
    }

    const low_bits: u64 = @bitCast(low);
    const high_bits: u64 = @bitCast(high);
    const error_round = (low_bits +% 6) & (@as(u64, std.math.maxInt(u64)) >> 12);
    const exponent_difference = ((high_bits -% low_bits) >> 52) & 0x7ff;
    if (error_round <= 12 or exponent_difference > 102) {
        return tanpiDatabase(value, result);
    }
    return result;
}
