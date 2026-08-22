// Correctly rounded binary64 sine for nearest-even rounding.
//
// Copyright (c) 2022-2025 Paul Zimmermann and Tom Hubrecht.
//
// This code is derived from CORE-MATH src/binary64/sin/sin.c at revision
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

const data = @import("scalar_sind_data.zig");
const Dint = data.Dint;

const sign_mask: u64 = 0x8000_0000_0000_0000;
const magnitude_mask: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits: u64 = infinity_bits | quiet_nan_bit;

const zero: Dint = .{ .hi = 0, .lo = 0, .ex = -1076, .sgn = 0 };
const magic: Dint = .{ .hi = 0x8000_0000_0000_0000, .lo = 0, .ex = -10, .sgn = 0 };

const DoubleDouble = struct {
    high: f64,
    low: f64,
};

const FastResult = struct {
    high: f64,
    low: f64,
    error_bound: f64,
};

/// Return correctly rounded sine with the caller's nearest-even mode.
///
/// The caller must set the floating-point environment to nearest-even before
/// calling this function.
pub fn sind(value: f64) f64 {
    @setFloatMode(.strict);

    const bits: u64 = @bitCast(value);
    const magnitude = bits & magnitude_mask;

    if (magnitude > infinity_bits) return @bitCast(bits | quiet_nan_bit);
    if (magnitude == infinity_bits) return @bitCast(canonical_nan_bits);

    // In this interval, sin(x) rounds to x in nearest-even mode, including
    // every binary64 subnormal and both signed zeros.
    if (magnitude <= 0x3e57_1374_4912_3ef6) return value;

    const fast = sinFast(value);
    const left = fast.high + (fast.low - fast.error_bound);
    const right = fast.high + (fast.low + fast.error_bound);
    if (left == right) return left;

    return sinAccurate(value);
}

fn dintBits(value: Dint) u128 {
    return (@as(u128, value.hi) << 64) | value.lo;
}

fn setDintBits(value: *Dint, bits: u128) void {
    value.hi = @truncate(bits >> 64);
    value.lo = @truncate(bits);
}

fn dintZero(value: Dint) bool {
    return value.hi == 0;
}

fn compareDintAbs(a: Dint, b: Dint) i2 {
    if (dintZero(a)) return if (dintZero(b)) 0 else -1;
    if (dintZero(b)) return 1;
    if (a.ex < b.ex) return -1;
    if (a.ex > b.ex) return 1;
    const a_bits = dintBits(a);
    const b_bits = dintBits(b);
    if (a_bits < b_bits) return -1;
    if (a_bits > b_bits) return 1;
    return 0;
}

fn addDint(result: *Dint, a_input: *const Dint, b_input: *const Dint) void {
    var a = a_input.*;
    var b = b_input.*;

    if (dintBits(a) == 0) {
        result.* = b;
        return;
    }

    const ordering = compareDintAbs(a, b);
    if (ordering == 0) {
        if (a.sgn != b.sgn) {
            result.* = zero;
        } else {
            result.* = a;
            result.ex += 1;
        }
        return;
    }
    if (ordering < 0) {
        const temporary = a;
        a = b;
        b = temporary;
    }

    const a_bits = dintBits(a);
    var b_bits = dintBits(b);
    const exponent_shift: u64 = @bitCast(a.ex - b.ex);
    if (exponent_shift != 0) {
        b_bits = if (exponent_shift < 128)
            b_bits >> @as(u7, @intCast(exponent_shift))
        else
            0;
    }

    var result_bits: u128 = undefined;
    result.ex = a.ex;
    if (a.sgn != b.sgn) {
        result_bits = a_bits - b_bits;
        var leading_zero_count: u8 = @intCast(@clz(result_bits));
        if (leading_zero_count != 0) {
            const shift: u7 = @intCast(leading_zero_count);
            if (exponent_shift == 1) {
                result_bits = (a_bits << shift) -% (dintBits(b) << (shift - 1));
            } else {
                result_bits = (a_bits << shift) -% (b_bits << shift);
            }
            result.ex -= leading_zero_count;
            leading_zero_count = @intCast(@clz(@as(u64, @truncate(result_bits >> 64))));
        }
        result_bits <<= @as(u7, @intCast(leading_zero_count));
        result.ex -= leading_zero_count;
    } else {
        result_bits = a_bits +% b_bits;
        if (result_bits < a_bits) {
            result_bits = (@as(u128, 1) << 127) | (result_bits >> 1);
            result.ex += 1;
        }
    }

    result.sgn = a.sgn;
    setDintBits(result, result_bits);
}

fn mulDint(result: *Dint, a_input: *const Dint, b_input: *const Dint) void {
    const a = a_input.*;
    const b = b_input.*;
    const middle_1 = @as(u128, a.hi) * b.lo;
    const middle_2 = @as(u128, a.lo) * b.hi;
    var result_bits = @as(u128, a.hi) * b.hi;
    result_bits +%= (middle_1 >> 64) + (middle_2 >> 64);

    const exponent_carry: u1 = @truncate(result_bits >> 127);
    result_bits <<= @as(u1, 1 - exponent_carry);
    result.ex = a.ex + b.ex + exponent_carry - 1;
    result.sgn = a.sgn ^ b.sgn;
    setDintBits(result, result_bits);
}

fn mulDint21(result: *Dint, a_input: *const Dint, b_input: *const Dint) void {
    const a = a_input.*;
    const b = b_input.*;
    const high_product = @as(u128, a.hi) * b.hi;
    const low_product = @as(u128, a.lo) * b.hi;
    var result_bits = high_product +% (low_product >> 64);

    const exponent_carry: u1 = @truncate(result_bits >> 127);
    result_bits <<= @as(u1, 1 - exponent_carry);
    result.ex = a.ex + b.ex + exponent_carry - 1;
    result.sgn = a.sgn ^ b.sgn;
    setDintBits(result, result_bits);
}

fn dintFromDouble(value: f64) Dint {
    const bits: u64 = @bitCast(value);
    const exponent_field = (bits >> 52) & 0x7ff;
    var exponent: i64 = @as(i64, @intCast(exponent_field)) - 0x3fe;
    var significand = (bits & 0x000f_ffff_ffff_ffff) +
        (if (exponent_field != 0) @as(u64, 1) << 52 else 0);
    const leading_zero_count: u7 = @intCast(@clz(significand));
    significand <<= @as(u6, @intCast(leading_zero_count));
    if (leading_zero_count > 11) exponent -= leading_zero_count - 12;
    return .{
        .hi = significand,
        .lo = 0,
        .ex = exponent,
        .sgn = @truncate(bits >> 63),
    };
}

fn roundRightEven(value: u128, shift_count: i64) u128 {
    if (shift_count <= 0) return value << @as(u7, @intCast(-shift_count));
    if (shift_count > 128) return 0;
    if (shift_count == 128) {
        const round_bit = (value >> 127) != 0;
        const sticky = (value & ((@as(u128, 1) << 127) - 1)) != 0;
        return @intFromBool(round_bit and sticky);
    }

    const shift: u7 = @intCast(shift_count);
    var rounded = value >> shift;
    const round_bit = ((value >> (shift - 1)) & 1) != 0;
    const sticky_mask = if (shift == 1)
        0
    else
        (@as(u128, 1) << (shift - 1)) - 1;
    const sticky = (value & sticky_mask) != 0;
    if (round_bit and (sticky or (rounded & 1) != 0)) rounded += 1;
    return rounded;
}

fn dintToDouble(value: Dint) f64 {
    @setFloatMode(.strict);
    const sign_bits = @as(u64, value.sgn) << 63;
    if (dintZero(value)) return @bitCast(sign_bits);

    if (value.ex >= -1021) {
        var exponent = value.ex;
        var significand = value.hi >> 11;
        const round_bit = ((value.hi >> 10) & 1) != 0;
        const sticky = (value.hi & 0x3ff) != 0 or value.lo != 0;
        if (round_bit and (sticky or (significand & 1) != 0)) {
            significand += 1;
            if (significand == (@as(u64, 1) << 53)) {
                significand >>= 1;
                exponent += 1;
            }
        }
        if (exponent > 1024) return @bitCast(sign_bits | infinity_bits);
        const exponent_bits: u64 = @intCast(exponent + 1022);
        return @bitCast(sign_bits | (exponent_bits << 52) |
            (significand & 0x000f_ffff_ffff_ffff));
    }

    // Convert directly to units of 2^-1074. This is independent of the host
    // floating-point environment and implements midpoint-to-even explicitly.
    const shift_count = -value.ex - 946;
    const subnormal: u64 = @truncate(roundRightEven(dintBits(value), shift_count));
    return @bitCast(sign_bits | subnormal);
}

fn normalize(value: *Dint) void {
    if (value.hi != 0) {
        const leading_zero_count: u7 = @intCast(@clz(value.hi));
        if (leading_zero_count != 0) {
            const shift: u6 = @intCast(leading_zero_count);
            value.hi = (value.hi << shift) |
                (value.lo >> @as(u6, @intCast(@as(u7, 64) - shift)));
            value.lo <<= shift;
        }
        value.ex -= leading_zero_count;
    } else if (value.lo != 0) {
        const leading_zero_count: u7 = @intCast(@clz(value.lo));
        value.hi = value.lo << @as(u6, @intCast(leading_zero_count));
        value.lo = 0;
        value.ex -= 64 + leading_zero_count;
    }
}

fn reduce(value: *Dint) void {
    var exponent = value.ex;
    if (exponent <= 1) {
        var product = @as(u128, value.hi) * data.t[1];
        const tiny: u64 = @truncate(product);
        value.lo = @truncate(product >> 64);
        product = @as(u128, value.hi) * data.t[0];
        const product_low: u64 = @truncate(product);
        const old_low = value.lo;
        value.lo +%= product_low;
        value.hi = @as(u64, @truncate(product >> 64)) +% @intFromBool(value.lo < old_low);
        exponent = value.ex;
        normalize(value);
        const normalization_shift: u7 = @intCast(exponent - value.ex);
        if (normalization_shift != 0) {
            value.lo |= tiny >> @as(u6, @intCast(@as(u7, 64) - normalization_shift));
        }
        return;
    }

    const table_index: usize = if (exponent < 127)
        0
    else
        @intCast(@divTrunc(exponent - 127 + 63, 64));
    var limbs: [5]u64 = undefined;

    var product = @as(u128, value.hi) * data.t[table_index + 3];
    limbs[0] = @truncate(product);
    limbs[1] = @truncate(product >> 64);
    product = @as(u128, value.hi) * data.t[table_index + 2];
    var addend: u64 = @truncate(product);
    const carry_1 = limbs[1] +% addend < limbs[1];
    limbs[1] +%= addend;
    limbs[2] = @as(u64, @truncate(product >> 64)) +% @intFromBool(carry_1);
    product = @as(u128, value.hi) * data.t[table_index + 1];
    addend = @truncate(product);
    const carry_2 = limbs[2] +% addend < limbs[2];
    limbs[2] +%= addend;
    limbs[3] = @as(u64, @truncate(product >> 64)) +% @intFromBool(carry_2);
    product = @as(u128, value.hi) * data.t[table_index];
    addend = @truncate(product);
    const carry_3 = limbs[3] +% addend < limbs[3];
    limbs[3] +%= addend;
    limbs[4] = @as(u64, @truncate(product >> 64)) +% @intFromBool(carry_3);

    const fraction_shift = exponent - @as(i64, @intCast(table_index)) * 64;
    var tiny: u64 = undefined;
    if (fraction_shift < 64) {
        const shift: u6 = @intCast(fraction_shift);
        value.hi = (limbs[4] << shift) |
            (limbs[3] >> @as(u6, @intCast(@as(u7, 64) - shift)));
        value.lo = (limbs[3] << shift) |
            (limbs[2] >> @as(u6, @intCast(@as(u7, 64) - shift)));
        tiny = (limbs[2] << shift) |
            (limbs[1] >> @as(u6, @intCast(@as(u7, 64) - shift)));
    } else if (fraction_shift == 64) {
        value.hi = limbs[3];
        value.lo = limbs[2];
        tiny = limbs[1];
    } else {
        const shift: u6 = @intCast(fraction_shift - 64);
        product = @as(u128, value.hi) * data.t[table_index + 4];
        const upper_product: u64 = @truncate(product >> 64);
        const carry = limbs[0] +% upper_product < limbs[0];
        limbs[0] +%= upper_product;
        if (carry) {
            limbs[1] +%= 1;
            if (limbs[1] == 0) {
                limbs[2] +%= 1;
                if (limbs[2] == 0) {
                    limbs[3] +%= 1;
                    if (limbs[3] == 0) limbs[4] +%= 1;
                }
            }
        }
        value.hi = (limbs[3] << shift) |
            (limbs[2] >> @as(u6, @intCast(@as(u7, 64) - shift)));
        value.lo = (limbs[2] << shift) |
            (limbs[1] >> @as(u6, @intCast(@as(u7, 64) - shift)));
        tiny = (limbs[1] << shift) |
            (limbs[0] >> @as(u6, @intCast(@as(u7, 64) - shift)));
    }

    value.ex = 0;
    normalize(value);
    if (value.ex < 0) value.lo |= tiny >> @as(u6, @intCast(64 + value.ex));
}

fn reduce2(value: *Dint) usize {
    if (value.ex <= -11) return 0;
    const shift: u6 = @intCast(53 - value.ex);
    const index: usize = @intCast(value.hi >> shift);
    value.hi &= (@as(u64, 1) << shift) - 1;
    normalize(value);
    return index;
}

fn setDoubleDouble(high: *f64, low: *f64, high_limb_input: u64, low_limb_input: u64) void {
    @setFloatMode(.strict);
    var high_limb = high_limb_input;
    var low_limb = low_limb_input;
    if (high_limb != 0) {
        const leading_zero_count: u7 = @intCast(@clz(high_limb));
        if (leading_zero_count != 0) {
            const shift: u6 = @intCast(leading_zero_count);
            high_limb = (high_limb << shift) |
                (low_limb >> @as(u6, @intCast(@as(u7, 64) - shift)));
            low_limb <<= shift;
        }
        const exponent_field: u64 = 0x3fe - @as(u64, leading_zero_count);
        high.* = @bitCast((exponent_field << 52) | ((high_limb << 1) >> 12));
        low_limb = (high_limb << 53) | (low_limb >> 11);
        if (low_limb != 0) {
            const low_zero_count: u7 = @intCast(@clz(low_limb));
            low_limb <<= @as(u6, @intCast(low_zero_count));
            const low_exponent: u64 = exponent_field - 53 - low_zero_count;
            low.* = @bitCast((low_exponent << 52) | ((low_limb << 1) >> 12));
        } else {
            low.* = 0.0;
        }
    } else if (low_limb != 0) {
        const leading_zero_count: u7 = @intCast(@clz(low_limb));
        const exponent_field: u64 = 0x3fe - 64 - @as(u64, leading_zero_count);
        const first_shift: u7 = leading_zero_count + 1;
        low_limb = if (first_shift == 64) 0 else low_limb << @as(u6, @intCast(first_shift));
        high.* = @bitCast((exponent_field << 52) | (low_limb >> 12));
        low_limb <<= 52;
        if (low_limb != 0) {
            const low_zero_count: u7 = @intCast(@clz(low_limb));
            const second_shift: u7 = low_zero_count + 1;
            low_limb = if (second_shift == 64) 0 else low_limb << @as(u6, @intCast(second_shift));
            const low_exponent: u64 = exponent_field - 64 - low_zero_count;
            low.* = @bitCast((low_exponent << 52) | (low_limb >> 12));
        } else {
            low.* = 0.0;
        }
    } else {
        high.* = 0.0;
        low.* = 0.0;
    }
}

fn reduceFast(high: *f64, low: *f64, value: f64, error_bound: *f64) usize {
    @setFloatMode(.strict);
    if (value <= 0x1.921fb54442d17p+2) {
        const product = aMul(0x1.45f306dc9c883p-3, value);
        high.* = product.high;
        low.* = @mulAdd(f64, -0x1.6b01ec5417056p-57, value, product.low);
        error_bound.* = 0x1.d9p-105 * high.*;
    } else {
        const bits: u64 = @bitCast(value);
        var exponent: i64 = @intCast((bits >> 52) & 0x7ff);
        const significand = (@as(u64, 1) << 52) | (bits & 0x000f_ffff_ffff_ffff);
        var limbs: [3]u64 = undefined;
        var product: u128 = undefined;

        if (exponent <= 1074) {
            product = @as(u128, significand) * data.t[1];
            limbs[0] = @truncate(product);
            limbs[1] = @truncate(product >> 64);
            product = @as(u128, significand) * data.t[0];
            const addend: u64 = @truncate(product);
            const carry = limbs[1] +% addend < limbs[1];
            limbs[1] +%= addend;
            limbs[2] = @as(u64, @truncate(product >> 64)) +% @intFromBool(carry);
            exponent = 1075 - exponent;
        } else {
            const table_index: usize = @intCast(@divTrunc(exponent - 1138 + 63, 64));
            product = @as(u128, significand) * data.t[table_index + 2];
            limbs[0] = @truncate(product);
            limbs[1] = @truncate(product >> 64);
            product = @as(u128, significand) * data.t[table_index + 1];
            var addend: u64 = @truncate(product);
            var carry = limbs[1] +% addend < limbs[1];
            limbs[1] +%= addend;
            limbs[2] = @as(u64, @truncate(product >> 64)) +% @intFromBool(carry);
            product = @as(u128, significand) * data.t[table_index];
            addend = @truncate(product);
            carry = limbs[2] +% addend < limbs[2];
            limbs[2] +%= addend;
            exponent = 1139 + @as(i64, @intCast(table_index << 6)) - exponent;
        }

        if (exponent == 64) {
            limbs[0] = limbs[1];
            limbs[1] = limbs[2];
        } else {
            const shift: u6 = @intCast(exponent);
            limbs[0] = (limbs[1] << @as(u6, @intCast(@as(u7, 64) - shift))) |
                (limbs[0] >> shift);
            limbs[1] = (limbs[2] << @as(u6, @intCast(@as(u7, 64) - shift))) |
                (limbs[1] >> shift);
        }
        setDoubleDouble(high, low, limbs[1], limbs[0]);
        error_bound.* = 0x1.01p-76;
    }

    const index: usize = @intFromFloat(high.* * 0x1p11);
    high.* = @mulAdd(f64, @as(f64, @floatFromInt(index)), -0x1p-11, high.*);
    return index;
}

fn aMul(a: f64, b: f64) DoubleDouble {
    @setFloatMode(.strict);
    const high = a * b;
    return .{ .high = high, .low = @mulAdd(f64, a, b, -high) };
}

fn sMul(a: f64, b_high: f64, b_low: f64) DoubleDouble {
    @setFloatMode(.strict);
    const product = aMul(a, b_high);
    return .{
        .high = product.high,
        .low = @mulAdd(f64, a, b_low, product.low),
    };
}

fn dMul(a_high: f64, a_low: f64, b_high: f64, b_low: f64) DoubleDouble {
    @setFloatMode(.strict);
    const high_product = aMul(a_high, b_high);
    const middle = @mulAdd(f64, a_low, b_high, high_product.low);
    return .{
        .high = high_product.high,
        .low = @mulAdd(f64, a_high, b_low, middle),
    };
}

fn fastTwoSum(a: f64, b: f64) DoubleDouble {
    @setFloatMode(.strict);
    const high = a + b;
    const error_value = high - a;
    return .{ .high = high, .low = b - error_value };
}

fn evalSinFast(x_high: f64, x_low: f64, square_high: f64, square_low: f64) DoubleDouble {
    @setFloatMode(.strict);
    var high = data.ps_fast[4];
    high = @mulAdd(f64, high, square_high, data.ps_fast[3]);
    high = @mulAdd(f64, high, square_high, data.ps_fast[2]);
    var product = sMul(high, square_high, square_low);
    const sum = fastTwoSum(data.ps_fast[0], product.high);
    product.high = sum.high;
    product.low += data.ps_fast[1] + sum.low;
    return dMul(product.high, product.low, x_high, x_low);
}

fn evalCosFast(square_high: f64, square_low: f64) DoubleDouble {
    @setFloatMode(.strict);
    var high = data.pc_fast[4];
    high = @mulAdd(f64, high, square_high, data.pc_fast[3]);
    high = @mulAdd(f64, high, square_high, data.pc_fast[2]);
    var product = sMul(high, square_high, square_low);
    const sum = fastTwoSum(data.pc_fast[0], product.high);
    product.high = sum.high;
    product.low += data.pc_fast[1] + sum.low;
    return product;
}

fn sinFast(value: f64) FastResult {
    @setFloatMode(.strict);
    var negative = value < 0.0;
    var use_sine = true;
    const absolute_value = if (negative) -value else value;

    var high: f64 = undefined;
    var low: f64 = undefined;
    var reduction_error: f64 = undefined;
    var index = reduceFast(&high, &low, absolute_value, &reduction_error);

    negative = negative != ((index >> 10) != 0);
    index &= 0x3ff;
    use_sine = use_sine != ((index >> 9) != 0);
    index &= 0x1ff;

    if ((index & 0x100) != 0) {
        use_sine = !use_sine;
        index = 0x1ff - index;
        high = 0x1p-11 - high;
        low = -low;
    }

    high -= data.sc[index][0];
    var square = aMul(high, high);
    square.low = @mulAdd(f64, high + high, low, square.low);
    const sine = evalSinFast(high, low, square.high, square.low);
    const cosine = evalCosFast(square.high, square.low);
    const sign: f64 = if (negative) -1.0 else 1.0;

    if (use_sine) {
        const sine_scaled = sMul(sign * data.sc[index][2], sine.high, sine.low);
        const cosine_scaled = sMul(sign * data.sc[index][1], cosine.high, cosine.low);
        const sum = fastTwoSum(cosine_scaled.high, sine_scaled.high);
        return .{
            .high = sum.high,
            .low = sum.low + (sine_scaled.low + cosine_scaled.low),
            .error_bound = 0x1.55p-69 + reduction_error,
        };
    }

    const cosine_scaled = sMul(sign * data.sc[index][2], cosine.high, cosine.low);
    const sine_scaled = sMul(sign * data.sc[index][1], sine.high, sine.low);
    const sum = fastTwoSum(cosine_scaled.high, -sine_scaled.high);
    return .{
        .high = sum.high,
        .low = sum.low + (cosine_scaled.low - sine_scaled.low),
        .error_bound = 0x1.81p-69 + reduction_error,
    };
}

fn evalSinAccurate(x: *const Dint, square: *const Dint) Dint {
    var result: Dint = undefined;
    mulDint21(&result, square, &data.ps[5]);
    addDint(&result, &result, &data.ps[4]);
    mulDint(&result, &result, square);
    addDint(&result, &result, &data.ps[3]);
    mulDint(&result, &result, square);
    addDint(&result, &result, &data.ps[2]);
    mulDint(&result, &result, square);
    addDint(&result, &result, &data.ps[1]);
    mulDint(&result, &result, square);
    addDint(&result, &result, &data.ps[0]);
    mulDint(&result, &result, x);
    return result;
}

fn evalCosAccurate(square: *const Dint) Dint {
    var result: Dint = undefined;
    mulDint21(&result, square, &data.pc[5]);
    addDint(&result, &result, &data.pc[4]);
    mulDint(&result, &result, square);
    addDint(&result, &result, &data.pc[3]);
    mulDint(&result, &result, square);
    addDint(&result, &result, &data.pc[2]);
    mulDint(&result, &result, square);
    addDint(&result, &result, &data.pc[1]);
    mulDint(&result, &result, square);
    addDint(&result, &result, &data.pc[0]);
    return result;
}

fn sinAccurate(value: f64) f64 {
    @setFloatMode(.strict);
    const input_bits: u64 = @bitCast(value);
    const absolute_bits = input_bits & magnitude_mask;
    const absolute_value: f64 = @bitCast(absolute_bits);
    var reduced = dintFromDouble(absolute_value);
    reduce(&reduced);

    var negative = (input_bits & sign_mask) != 0;
    var use_sine = true;
    var index = reduce2(&reduced);

    if ((index & 0x400) != 0) {
        negative = !negative;
        index &= 0x3ff;
    }
    if ((index & 0x200) != 0) {
        use_sine = false;
        index &= 0x1ff;
    }
    if ((index & 0x100) != 0) {
        use_sine = !use_sine;
        reduced.sgn = 1;
        addDint(&reduced, &magic, &reduced);
        index = 0x1ff - index;
    }

    var square: Dint = undefined;
    mulDint(&square, &reduced, &reduced);
    var cosine = evalCosAccurate(&square);
    var sine = evalSinAccurate(&reduced, &square);

    if (use_sine) {
        mulDint(&cosine, &data.s[index], &cosine);
        mulDint(&sine, &data.c[index], &sine);
    } else {
        mulDint(&cosine, &data.c[index], &cosine);
        mulDint(&sine, &data.s[index], &sine);
        sine.sgn = 1 - sine.sgn;
    }
    addDint(&cosine, &cosine, &sine);

    const error_units: u64 = 41;
    const lower_low = cosine.lo -% error_units;
    const lower_high = cosine.hi -% @intFromBool(lower_low > cosine.lo);
    const upper_low = cosine.lo +% error_units;
    const upper_high = cosine.hi +% @intFromBool(upper_low < cosine.lo);
    if ((lower_high >> 10) != (upper_high >> 10)) {
        if (absolute_bits == 0x3ebe_0000_0000_01c2) {
            return @bitCast(@as(u64, 0x3ebd_ffff_ffff_f02e) |
                (input_bits & sign_mask));
        }
        if (absolute_bits == 0x7506_ac5b_262c_a1ff) {
            return @bitCast(@as(u64, 0x3ff0_0000_0000_0000) |
                (input_bits & sign_mask));
        }
    }

    if (negative) cosine.sgn = 1 - cosine.sgn;
    return dintToDouble(cosine);
}

// The cosine path is derived from CORE-MATH src/binary64/cos/cos.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and uses the MIT License.
fn cosFast(value: f64) FastResult {
    @setFloatMode(.strict);
    var negative = false;
    var use_cosine = true;

    var high: f64 = undefined;
    var low: f64 = undefined;
    var reduction_error: f64 = undefined;
    var index = reduceFast(&high, &low, value, &reduction_error);

    negative = negative != ((index >> 10) != 0);
    index &= 0x3ff;
    use_cosine = use_cosine != ((index >> 9) != 0);
    negative = negative != ((index >> 9) != 0);
    index &= 0x1ff;

    if ((index & 0x100) != 0) {
        use_cosine = !use_cosine;
        index = 0x1ff - index;
        high = 0x1p-11 - high;
        low = -low;
    }

    high -= data.sc[index][0];
    var square = aMul(high, high);
    square.low = @mulAdd(f64, high + high, low, square.low);
    const sine = evalSinFast(high, low, square.high, square.low);
    const cosine = evalCosFast(square.high, square.low);
    const sign: f64 = if (negative) -1.0 else 1.0;

    if (!use_cosine) {
        const sine_scaled = sMul(sign * data.sc[index][2], sine.high, sine.low);
        const cosine_scaled = sMul(sign * data.sc[index][1], cosine.high, cosine.low);
        const sum = fastTwoSum(cosine_scaled.high, sine_scaled.high);
        return .{
            .high = sum.high,
            .low = sum.low + (sine_scaled.low + cosine_scaled.low),
            .error_bound = 0x1.55p-69 + reduction_error,
        };
    }

    const cosine_scaled = sMul(sign * data.sc[index][2], cosine.high, cosine.low);
    const sine_scaled = sMul(sign * data.sc[index][1], sine.high, sine.low);
    const sum = fastTwoSum(cosine_scaled.high, -sine_scaled.high);
    return .{
        .high = sum.high,
        .low = sum.low + (cosine_scaled.low - sine_scaled.low),
        .error_bound = 0x1.81p-69 + reduction_error,
    };
}

fn cosAccurate(value: f64) f64 {
    @setFloatMode(.strict);
    const absolute_bits: u64 = @bitCast(value);
    var reduced = dintFromDouble(value);
    reduce(&reduced);

    var negative = false;
    var use_cosine = true;
    var index = reduce2(&reduced);

    if ((index & 0x400) != 0) {
        negative = true;
        index &= 0x3ff;
    }
    if ((index & 0x200) != 0) {
        negative = !negative;
        use_cosine = false;
        index &= 0x1ff;
    }
    if ((index & 0x100) != 0) {
        use_cosine = !use_cosine;
        reduced.sgn = 1;
        addDint(&reduced, &magic, &reduced);
        index = 0x1ff - index;
    }

    var square: Dint = undefined;
    mulDint(&square, &reduced, &reduced);
    var cosine = evalCosAccurate(&square);
    var sine = evalSinAccurate(&reduced, &square);

    if (!use_cosine) {
        mulDint(&cosine, &data.s[index], &cosine);
        mulDint(&sine, &data.c[index], &sine);
    } else {
        mulDint(&cosine, &data.c[index], &cosine);
        mulDint(&sine, &data.s[index], &sine);
        sine.sgn = 1 - sine.sgn;
    }
    addDint(&cosine, &cosine, &sine);

    const error_units: u64 = 41;
    const lower_low = cosine.lo -% error_units;
    const lower_high = cosine.hi -% @intFromBool(lower_low > cosine.lo);
    const upper_low = cosine.lo +% error_units;
    const upper_high = cosine.hi +% @intFromBool(upper_low < cosine.lo);
    if ((lower_high >> 10) != (upper_high >> 10)) {
        if (absolute_bits == 0x3e88_0000_0000_0009) {
            return @bitCast(@as(u64, 0x3fef_ffff_ffff_ff70));
        }
        if (absolute_bits == 0x3e98_0000_0000_0024) {
            return @bitCast(@as(u64, 0x3fef_ffff_ffff_fdc0));
        }
        if (absolute_bits == 0x3ea8_0000_0000_0090) {
            return @bitCast(@as(u64, 0x3fef_ffff_ffff_f700));
        }
        if (absolute_bits == 0x3eb2_0000_0000_00f3) {
            return @bitCast(@as(u64, 0x3fef_ffff_ffff_ebc0));
        }
        if (absolute_bits == 0x3eb8_0000_0000_0240) {
            return @bitCast(@as(u64, 0x3fef_ffff_ffff_dc00));
        }
    }

    if (negative) cosine.sgn = 1 - cosine.sgn;
    return dintToDouble(cosine);
}

/// Return correctly rounded cosine with the caller's nearest-even mode.
///
/// The caller must set the floating-point environment to nearest-even before
/// calling this function.
pub fn cosd(value: f64) f64 {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    const magnitude = bits & magnitude_mask;

    if (magnitude > infinity_bits) return @bitCast(bits | quiet_nan_bit);
    if (magnitude == infinity_bits) return @bitCast(canonical_nan_bits);

    // For this interval, cos(x) rounds to one in nearest-even mode.
    if (magnitude <= 0x3e46_a09e_667f_3bcc) {
        const absolute_value: f64 = @bitCast(magnitude);
        return @mulAdd(f64, absolute_value, -0x1p-28, 1.0);
    }

    const absolute_value: f64 = @bitCast(magnitude);
    const fast = cosFast(absolute_value);
    const left = fast.high + (fast.low - fast.error_bound);
    const right = fast.high + (fast.low + fast.error_bound);
    if (left == right) return left;

    return cosAccurate(absolute_value);
}
