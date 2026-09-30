// Correctly rounded Euclidean distance functions ported from CORE-MATH:
// https://gitlab.inria.fr/core-math/core-math.git
// Revision 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
//
// Copyright (c) 2022-2025 Alexei Sibidanov.
// Copyright (c) 2025-2026 Maxence Ponsardin and Paul Zimmermann.
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

const fraction_mask_d: u64 = 0x000f_ffff_ffff_ffff;

const magnitude_mask_h: u16 = 0x7fff;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;

const infinity_bits_h: u16 = 0x7c00;
const infinity_bits_f: u32 = 0x7f80_0000;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;

const quiet_nan_bit_h: u16 = 0x0200;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;

const hidden_bit_d: u64 = 0x0010_0000_0000_0000;
const exponent_bias_d: u64 = 0x3ff;
const exponent_mask_d: u64 = 0x7ff0_0000_0000_0000;

const TwoSum = struct {
    high: f64,
    low: f64,
};

fn fastTwoSum(x: f64, y: f64) TwoSum {
    @setFloatMode(.strict);
    const sum = x + y;
    const z = sum - x;
    return .{ .high = sum, .low = y - z };
}

fn quietNanH(bits: u16) f16 {
    return @bitCast(bits | quiet_nan_bit_h);
}

fn quietNanF(bits: u32) f32 {
    return @bitCast(bits | quiet_nan_bit_f);
}

fn quietNanD(bits: u64) f64 {
    return @bitCast(bits | quiet_nan_bit_d);
}

fn isNanH(bits: u16) bool {
    return (bits & magnitude_mask_h) > infinity_bits_h;
}

fn isNanD(bits: u64) bool {
    return (bits & magnitude_mask_d) > infinity_bits_d;
}

fn isSignalingNanH(bits: u16) bool {
    return isNanH(bits) and (bits & quiet_nan_bit_h) == 0;
}

fn isSignalingNanD(bits: u64) bool {
    return isNanD(bits) and (bits & quiet_nan_bit_d) == 0;
}

fn hypotSpecialH(x_bits: u16, y_bits: u16) ?f16 {
    const x_magnitude = x_bits & magnitude_mask_h;
    const y_magnitude = y_bits & magnitude_mask_h;
    const x_nan = x_magnitude > infinity_bits_h;
    const y_nan = y_magnitude > infinity_bits_h;
    const x_inf = x_magnitude == infinity_bits_h;
    const y_inf = y_magnitude == infinity_bits_h;
    if (x_nan and y_nan) return quietNanH(x_bits);
    if (x_nan) {
        if (isSignalingNanH(x_bits)) return quietNanH(x_bits);
        if (y_inf) return @bitCast(infinity_bits_h);
        return quietNanH(x_bits);
    }
    if (y_nan) {
        if (isSignalingNanH(y_bits)) return quietNanH(y_bits);
        if (x_inf) return @bitCast(infinity_bits_h);
        return quietNanH(y_bits);
    }
    if (x_inf or y_inf) return @bitCast(infinity_bits_h);
    return null;
}

noinline fn hypotSpecialF(x_bits: u32, y_bits: u32) f32 {
    const x_magnitude = x_bits & magnitude_mask_f;
    const y_magnitude = y_bits & magnitude_mask_f;

    const x_nan = x_magnitude > infinity_bits_f;
    const y_nan = y_magnitude > infinity_bits_f;
    if (x_nan and y_nan) return quietNanF(x_bits);
    if (x_nan) {
        if ((x_bits & quiet_nan_bit_f) == 0) return quietNanF(x_bits);
        if (y_magnitude == infinity_bits_f) return @bitCast(infinity_bits_f);
        return quietNanF(x_bits);
    }
    if (y_nan) {
        if ((y_bits & quiet_nan_bit_f) == 0) return quietNanF(y_bits);
        if (x_magnitude == infinity_bits_f) return @bitCast(infinity_bits_f);
        return quietNanF(y_bits);
    }
    return @bitCast(infinity_bits_f);
}

fn hypotSpecialD(x_bits: u64, y_bits: u64) ?f64 {
    const x_magnitude = x_bits & magnitude_mask_d;
    const y_magnitude = y_bits & magnitude_mask_d;
    const x_nan = x_magnitude > infinity_bits_d;
    const y_nan = y_magnitude > infinity_bits_d;
    const x_inf = x_magnitude == infinity_bits_d;
    const y_inf = y_magnitude == infinity_bits_d;
    if (x_nan and y_nan) return quietNanD(x_bits);
    if (x_nan) {
        if (isSignalingNanD(x_bits)) return quietNanD(x_bits);
        if (y_inf) return @bitCast(infinity_bits_d);
        return quietNanD(x_bits);
    }
    if (y_nan) {
        if (isSignalingNanD(y_bits)) return quietNanD(y_bits);
        if (x_inf) return @bitCast(infinity_bits_d);
        return quietNanD(y_bits);
    }
    if (x_inf or y_inf) return @bitCast(infinity_bits_d);
    return null;
}

// Correctly rounded binary16 hypot, ported from CORE-MATH
// src/binary16/hypot/hypotf16.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa under the MIT License.
pub fn hypoth(x: f16, y: f16) f16 {
    @setFloatMode(.strict);
    const x_bits: u16 = @bitCast(x);
    const y_bits: u16 = @bitCast(y);
    if (hypotSpecialH(x_bits, y_bits)) |special| return special;

    const x_magnitude = x_bits & magnitude_mask_h;
    const y_magnitude = y_bits & magnitude_mask_h;
    const larger_bits = @max(x_magnitude, y_magnitude);
    const smaller_bits = @min(x_magnitude, y_magnitude);
    const larger: f32 = @floatCast(@as(f16, @bitCast(larger_bits)));
    const smaller: f32 = @floatCast(@as(f16, @bitCast(smaller_bits)));

    // For n <= m/64, hypot(m, n) - m <= m/8192, below every upper f16 midpoint.
    // This includes positive zero, subnormals, and the finite-to-infinity boundary.
    if (smaller * 64.0 <= larger) return @bitCast(larger_bits);

    // Half products are exact in f32; only the sum and square root round.
    // Eight adjacent f32 values on each side enclose the exact root, including
    // at binade boundaries. Equal narrowed endpoints certify the f16 result.
    const root: f32 = @sqrt(larger * larger + smaller * smaller);
    const root_bits: u32 = @bitCast(root);
    const lower: f16 = @floatCast(@as(f32, @bitCast(root_bits - 8)));
    const upper: f16 = @floatCast(@as(f32, @bitCast(root_bits + 8)));
    if (@as(u16, @bitCast(lower)) == @as(u16, @bitCast(upper))) return lower;

    const x_abs: f16 = @bitCast(x_magnitude);
    const y_abs: f16 = @bitCast(y_magnitude);
    const x_wide: f64 = @floatCast(x_abs);
    const y_wide: f64 = @floatCast(y_abs);
    return @floatCast(@sqrt(x_wide * x_wide + y_wide * y_wide));
}

noinline fn hypotfRefine(
    root: f64,
    squared_sum: f64,
    x_squared: f64,
    y_squared: f64,
) f32 {
    @setFloatMode(.strict);
    const inverse_squared_sum = 0.5 / squared_sum;
    const squared_sum_error = (x_squared - squared_sum) + y_squared;
    const root_inverse = root * inverse_squared_sum;
    const squared_error = squared_sum_error - @mulAdd(f64, root, root, -squared_sum);
    const correction = root_inverse * squared_error;
    const refined = root + correction;
    const residual = correction + (root - refined);
    var refined_bits: u64 = @bitCast(refined);
    if ((refined_bits & 0x0fff_ffff) == 0) {
        if (residual > 0.0) refined_bits +%= 1;
        if (residual < 0.0) refined_bits -%= 1;
    }
    return @floatCast(@as(f64, @bitCast(refined_bits)));
}

// Correctly rounded binary32 hypot, ported from CORE-MATH
// src/binary32/hypot/hypotf.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa under the MIT License.
pub fn hypotf(x: f32, y: f32) f32 {
    @setFloatMode(.strict);
    const x_bits: u32 = @bitCast(x);
    const y_bits: u32 = @bitCast(y);
    const x_magnitude = x_bits & magnitude_mask_f;
    const y_magnitude = y_bits & magnitude_mask_f;
    if (@max(x_magnitude, y_magnitude) >= infinity_bits_f) return hypotSpecialF(x_bits, y_bits);

    const xd: f64 = x;
    const yd: f64 = y;
    const x_squared = xd * xd;
    const y_squared = yd * yd;
    const squared_sum = x_squared + y_squared;

    // Subnormal outputs have subnormal inputs. In units of 2^-149 their
    // squared sum is an exact integer; its root is more than 2^-27 from
    // a half-integer, while binary64 rounding is below 2^-30 units.
    // Thus only normal binary32 midpoints can need refinement.

    const root: f64 = @sqrt(squared_sum);
    const root_bits: u64 = @bitCast(root);
    const result: f32 = @floatCast(root);
    if (root_bits > 0x47ef_ffff_e000_0000) return result;
    if (((root_bits +% 0x1000_0001) & 0x1fff_ffff) > 2) return result;

    const larger_squared = @max(x_squared, y_squared);
    const smaller_squared = @min(x_squared, y_squared);
    const result_wide: f64 = @floatCast(result);
    if ((result_wide * result_wide - larger_squared) - smaller_squared == 0.0) return result;
    return hypotfRefine(root, squared_sum, larger_squared, smaller_squared);
}

fn hypotdDenorm(x_bits: u64, y_bits: u64) f64 {
    @setFloatMode(.strict);
    const x_value: f64 = @floatFromInt(x_bits);
    const y_value: f64 = @floatFromInt(y_bits);
    const approximate_root: f64 = @sqrt(x_value * x_value + y_value * y_value);
    // Both inputs are subnormal significands below 2^52, so each square rounds
    // to at most 2^104 - 2^51, their sum is exact and stays below 2^105, and
    // the root is nonnegative and below 2^52.5. The clamp and the guard make
    // that conversion bound explicit.
    const root_bound: f64 = 0x1.8p+52;
    const bounded_root: f64 = @min(@max(approximate_root, 0.0), root_bound);
    const root_value: u64 = if (bounded_root >= 0.0 and bounded_root <= root_bound)
        @intFromFloat(bounded_root)
    else
        unreachable;
    var root: u64 = root_value << 1;
    const x_scaled: u128 = @as(u128, x_bits) << 1;
    const y_scaled: u128 = @as(u128, y_bits) << 1;
    const squared_sum: u128 = x_scaled * x_scaled + y_scaled * y_scaled;

    while (squared_sum > @as(u128, root) * root + 2 * root) {
        root += 1;
    }
    while (squared_sum < @as(u128, root) * root) {
        root -= 1;
    }

    const remainder = squared_sum - @as(u128, root) * root;
    const round_bit: u64 = root & 1;
    const sticky_bit = remainder != 0;
    var result = root >> 1;
    if (sticky_bit) {
        result += round_bit;
    } else if (round_bit != 0) {
        result += result & 1;
    }
    return @bitCast(result);
}

fn hypotdHard(x: f64, y: f64) f64 {
    @setFloatMode(.strict);
    const x_bits: u64 = @bitCast(x);
    const y_bits: u64 = @bitCast(y);
    var x_mantissa = (x_bits & fraction_mask_d) | hidden_bit_d;
    const y_mantissa = (y_bits & fraction_mask_d) | hidden_bit_d;
    const x_exponent: u64 = x_bits >> 52;
    const y_exponent: u64 = y_bits >> 52;
    const root_bits: u64 = @bitCast(@sqrt(x * x + y * y));
    var root_mantissa = (root_bits & fraction_mask_d) | hidden_bit_d;
    var root_exponent: i64 = @as(i64, @intCast(root_bits >> 52)) - 0x3ff;

    var iteration: u8 = 0;
    while (iteration < 3) : (iteration += 1) {
        if (root_mantissa == hidden_bit_d) {
            root_mantissa = 0x001f_ffff_ffff_ffff;
            root_exponent -= 1;
        } else {
            root_mantissa -= 1;
        }
    }

    x_mantissa <<= 2;
    const exponent_difference: i64 = @as(i64, @intCast(x_exponent)) -
        @as(i64, @intCast(y_exponent));
    const alignment = 2 - exponent_difference;
    var squared_sum: u128 = @as(u128, x_mantissa) * x_mantissa;
    if (alignment >= 0) {
        const shifted_y_mantissa = y_mantissa << @as(u6, @intCast(alignment));
        squared_sum += @as(u128, shifted_y_mantissa) * shifted_y_mantissa;
    } else {
        const y_square = @as(u128, y_mantissa) * y_mantissa;
        const shift: u7 = @intCast(-2 * alignment);
        squared_sum += y_square >> shift;
        if ((y_square << @as(u7, @intCast(128 - @as(u8, shift)))) != 0) squared_sum |= 1;
    }

    const shift: u6 = @intCast(2 + root_exponent);
    while (true) {
        const increment: u64 = 1 + @as(u64, @intFromBool(root_mantissa >= (1 << 53)));
        root_mantissa += increment;
        const candidate = @as(u128, root_mantissa) << shift;
        if (squared_sum <= candidate * candidate) break;
    }

    if (squared_sum != blk: {
        const candidate = @as(u128, root_mantissa) << shift;
        break :blk candidate * candidate;
    }) {
        const previous = (@as(u128, root_mantissa) << shift) -
            (@as(u128, 1) << (shift - @intFromBool(root_mantissa <= (1 << 53))));
        const previous_square = previous * previous;
        if (squared_sum < previous_square) {
            root_mantissa -= 1;
        } else if (squared_sum == previous_square) {
            root_mantissa -= root_mantissa & 1;
        }
    }

    if (root_mantissa >= (1 << 53)) {
        root_mantissa >>= 1;
        root_exponent += 1;
    }
    const output_exponent: u64 = @intCast(@as(i64, @intCast(x_exponent)) - 1 + root_exponent);
    const output_bits = (output_exponent << 52) + root_mantissa;
    return @bitCast(output_bits);
}

// Correctly rounded binary64 hypot, ported from CORE-MATH
// src/binary64/hypot/hypot.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa under the MIT License.
pub fn hypotd(x: f64, y: f64) f64 {
    @setFloatMode(.strict);
    const x_bits: u64 = @bitCast(x);
    const y_bits: u64 = @bitCast(y);
    if (hypotSpecialD(x_bits, y_bits)) |special| return special;

    var x_magnitude = x_bits & magnitude_mask_d;
    var y_magnitude = y_bits & magnitude_mask_d;
    if (x_magnitude < y_magnitude) {
        const temporary = x_magnitude;
        x_magnitude = y_magnitude;
        y_magnitude = temporary;
    }
    const original_y_magnitude = y_magnitude;
    if (y_magnitude == 0) return @bitCast(x_magnitude);
    if ((y_magnitude >> 52) == 0) {
        if ((x_magnitude >> 52) == 0) return hypotdDenorm(x_magnitude, y_magnitude);
        const leading_zeroes: u7 = @clz(y_magnitude);
        var normalized_y = y_magnitude << @as(u6, @intCast(leading_zeroes - 11));
        normalized_y &= fraction_mask_d;
        normalized_y -%= @as(u64, leading_zeroes - 12) << 52;
        y_magnitude = normalized_y;
    }

    if (x_magnitude -% y_magnitude > (27 << 52)) {
        const x_value: f64 = @bitCast(x_magnitude);
        const y_value: f64 = @bitCast(original_y_magnitude);
        return @mulAdd(f64, 0x1p-27, y_value, x_value);
    }

    const offset: i64 = @as(i64, @intCast(exponent_bias_d << 52)) -
        @as(i64, @intCast(x_magnitude & exponent_mask_d));
    const offset_bits: u64 = @bitCast(offset);
    const normalized_x_bits = x_magnitude +% offset_bits;
    const normalized_y_bits = y_magnitude +% offset_bits;
    const normalized_x: f64 = @bitCast(normalized_x_bits);
    const normalized_y: f64 = @bitCast(normalized_y_bits);

    const x_squared = normalized_x * normalized_x;
    const x_squared_error = @mulAdd(f64, normalized_x, normalized_x, -x_squared);
    const y_squared = normalized_y * normalized_y;
    const y_squared_error = @mulAdd(f64, normalized_y, normalized_y, -y_squared);
    const squared_sum = x_squared + y_squared;
    const inverse_squared_sum = 0.5 / squared_sum;
    const squared_sum_error = ((x_squared - squared_sum) + y_squared) +
        (x_squared_error + y_squared_error);
    var root = @sqrt(squared_sum);
    const root_inverse = root * inverse_squared_sum;
    const squared_error = squared_sum_error - @mulAdd(f64, root, root, -squared_sum);
    const correction = root_inverse * squared_error;
    const sum = fastTwoSum(root, correction);
    root = sum.high;
    const residual = sum.low;

    const root_bits: u64 = @bitCast(root);
    const residual_bits: u64 = @bitCast(@abs(residual));
    const root_exponent = root_bits & exponent_mask_d;
    const aidr = residual_bits +% (0x3fe << 52) -% root_exponent;
    const midpoint = (aidr -% 0x3c90_0000_0000_0000 +% 16) >> 5;
    if (midpoint == 0 or aidr < 0x39b0_0000_0000_0000 or aidr > 0x3c9f_ffff_ffff_ff80) {
        root = hypotdHard(normalized_x, normalized_y);
    }

    const result_bits: u64 = @intCast(
        @as(i64, @intCast(@as(u64, @bitCast(root)))) - offset,
    );
    if (result_bits >= infinity_bits_d) return @bitCast(infinity_bits_d);
    return @bitCast(result_bits);
}
