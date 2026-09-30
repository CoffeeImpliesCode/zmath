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

    if (magnitude <= 0x3fc0_0000_0000_0000) {
        // For |x| <= 1/8, the degree-eleven Taylor tail is below
        // |x|^13/13!. Evaluation and truncation fit within u*|x^3|,
        // where u=2^-53. FastTwoSum retains the correction exactly.
        const square = value * value;
        const cubic = value * square;
        var polynomial: f64 = -0x1.ae64567f544e4p-26;
        polynomial = @mulAdd(f64, polynomial, square, 0x1.71de3a556c734p-19);
        polynomial = @mulAdd(f64, polynomial, square, -0x1.a01a01a01a01ap-13);
        polynomial = @mulAdd(f64, polynomial, square, 0x1.1111111111111p-7);
        polynomial = @mulAdd(f64, polynomial, square, -0x1.5555555555555p-3);
        const approximation = fastTwoSum(value, cubic * polynomial);
        const error_bound = @abs(cubic) * 0x1p-53;
        const left = approximation.high + (approximation.low - error_bound);
        const right = approximation.high + (approximation.low + error_bound);
        if (left == right) return left;
        return sinAccurate(value);
    }

    const absolute_value: f64 = @bitCast(magnitude);
    const fast = sinFast(absolute_value, (bits & sign_mask) != 0);
    const left = fast.high + (fast.low - fast.error_bound);
    const right = fast.high + (fast.low + fast.error_bound);
    if (left != right) {
        @branchHint(.cold);
        return sinAccurate(value);
    }
    return left;
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
    // `dintToDouble` is the only caller and reaches this line under
    // `value.ex < -1021`, so `-value.ex - 946 >= 76` and this branch is empty
    // today. `@abs` keeps it total anyway: it takes the magnitude without the
    // `-x` overflow at `minInt(i64)`, leaving only the `u7` narrowing, which no
    // left shift of 128 or more can satisfy.
    if (shift_count <= 0) return value << @as(u7, @intCast(@abs(shift_count)));
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
    // Keep the leading 53 fraction bits in high and the next 53 in low.
    // Each integer conversion is exact; the discarded suffix is below
    // 2^-106 turns.
    const high_integer: u64 = high_limb_input >> 11;
    const low_integer: u64 = ((high_limb_input & 0x7ff) << 42) | (low_limb_input >> 22);
    high.* = @as(f64, @floatFromInt(high_integer)) * 0x1p-53;
    low.* = @as(f64, @floatFromInt(low_integer)) * 0x1p-106;
}

// The reducer guarantees 0 <= value < 1. Decode floor(value * 2^11)
// directly so ReleaseSafe does not need a checked float-to-integer conversion.
inline fn reductionIndex(value: f64) u16 {
    const bits: u64 = @bitCast(value);
    const exponent: u11 = @truncate(bits >> 52);
    if (exponent < 0x3f4) return 0;
    if (exponent >= 0x3ff) unreachable;

    const significand = (@as(u64, 1) << 52) | (bits & 0x000f_ffff_ffff_ffff);
    const shift: u6 = @intCast(0x428 - exponent);
    return @truncate(significand >> shift);
}

// Multiplication by either +1 or -1 changes only the sign bit for finite
// binary64 values. Express that exact operation without an FP multiply.
inline fn negateIf(value: f64, negative: bool) f64 {
    const bits: u64 = @bitCast(value);
    const sign = @as(u64, @intFromBool(negative)) << 63;
    return @bitCast(bits ^ sign);
}

inline fn reduceFast(
    comptime tangent_precision: bool,
    high: *f64,
    low: *f64,
    value: f64,
    error_bound: *f64,
) u16 {
    @setFloatMode(.strict);
    if (value <= 0x1.921fb54442d17p+2) {
        // Here product_high is below one, so index extraction and residual
        // subtraction can precede the split-tail accumulation.
        const product_high = 0x1.45f306dc9c883p-3 * value;
        var index: u16 = undefined;
        var product_low: f64 = undefined;
        if (tangent_precision) {
            index = reductionIndex(product_high);
            product_low = @mulAdd(
                f64,
                0x1.45f306dc9c883p-3,
                value,
                -product_high,
            );
            high.* = @mulAdd(
                f64,
                @as(f64, @floatFromInt(index)),
                -0x1p-11,
                product_high,
            );
        } else {
            // The sine and cosine callers enter here with value > 1/8, so
            // 2^-6 < product_high < 1. Removing the encoded fraction bits
            // forms the exact floor grid without integer conversion.
            const product_bits: u64 = @bitCast(product_high);
            const product_exponent: u11 = @truncate(product_bits >> 52);
            // The shift lies in [42, 47], so truncation loses no bits.
            const fraction_bit_count: u6 = @truncate(0x428 -% product_exponent);
            const significand = (@as(u64, 1) << 52) |
                (product_bits & 0x000f_ffff_ffff_ffff);
            const grid_units = significand >> fraction_bit_count;
            const grid_prefix = product_bits >> fraction_bit_count;
            index = @truncate(grid_units & 0x7ff);
            const grid_high: f64 = @bitCast(grid_prefix << fraction_bit_count);
            high.* = product_high - grid_high;
            product_low = @mulAdd(
                f64,
                0x1.45f306dc9c883p-3,
                value,
                -product_high,
            );
        }
        low.* = @mulAdd(f64, -0x1.6b01ec5417056p-57, value, product_low);
        // The reducer error is in turns. Since |d sin(2*pi*x)/dx| and
        // |d cos(2*pi*x)/dx| are below 2*pi < 8, scale it to an output
        // bound for the sine and cosine callers.
        if (!tangent_precision) error_bound.* = 0x1.d9p-102 * product_high;
        return index;
    } else if (value <= 0x1p+26) {
        // For 2*pi < value <= 2^26, the split reciprocal gives a turn
        // error below 2^-80.768.
        const product_high = 0x1.45f306dc9c883p-3 * value;
        var provisional_index: u16 = undefined;
        var grid_high: f64 = undefined;
        if (tangent_precision) {
            // Retain tangent's branchless compensation of the rounded grid.
            const grid_bias: f64 = 0x1.8p41;
            const biased = product_high + grid_bias;
            const biased_bits: u64 = @bitCast(biased);
            const bias_bits: u64 = @bitCast(grid_bias);
            var grid_units: u64 = biased_bits -% bias_bits;
            const rounded_grid = biased - grid_bias;
            const overshoot = rounded_grid > product_high;
            grid_units -%= @intFromBool(overshoot);
            const correction: f64 = @bitCast(
                @as(u64, @intFromBool(overshoot)) * 0x3f40_0000_0000_0000,
            );
            grid_high = rounded_grid - correction;
            provisional_index = @truncate(grid_units & 0x7ff);
        } else {
            // Here 1 <= product_high < 2^24. Removing the encoded fraction
            // bits below 2^-11 forms the same floor grid directly. The implicit
            // leading bit is above the low 11 grid bits, which hold the index.
            const product_bits: u64 = @bitCast(product_high);
            const product_exponent: u11 = @truncate(product_bits >> 52);
            // The shift lies in [18, 41], so truncation loses no bits.
            const fraction_bit_count: u6 = @truncate(0x428 -% product_exponent);
            const grid_prefix = product_bits >> fraction_bit_count;
            provisional_index = @truncate(grid_prefix & 0x7ff);
            grid_high = @bitCast(grid_prefix << fraction_bit_count);
        }
        const provisional_high = product_high - grid_high;
        const product_low = @mulAdd(
            f64,
            0x1.45f306dc9c883p-3,
            value,
            -product_high,
        );
        const fractional_low = @mulAdd(
            f64,
            -0x1.6b01ec5417056p-57,
            value,
            product_low,
        );
        // The rounded multiplication residual is at most 2^-30. The
        // split-tail contribution keeps |fractional_low| below 2^-29.
        // Thus a 2^-27 margin prevents changes to the table index and
        // crossings of either neighboring grid boundary.
        const index_margin = 0x1p-27;
        if (!tangent_precision and value <= 0x1p+6 and
            provisional_high >= index_margin and
            provisional_high <= 0x1p-11 - index_margin)
        {
            high.* = provisional_high;
            low.* = fractional_low;
            error_bound.* = 0x1p-67;
        } else {
            // The reflection step needs a coarse high-part lattice so
            // 2^-11-high is exact. Adding one is exact here; normalize
            // on its binade, then remove it without losing the tail.
            const reduced = fastTwoSum(1.0 + provisional_high, fractional_low);
            high.* = reduced.high - 1.0;
            low.* = reduced.low;
            error_bound.* = if (tangent_precision) 0x1p-80 else 0x1p-77;
        }
        return provisional_index;
    } else {
        const bits: u64 = @bitCast(value);
        var exponent: i64 = @intCast((bits >> 52) & 0x7ff);
        const significand = (@as(u64, 1) << 52) | (bits & 0x000f_ffff_ffff_ffff);
        var limbs: [3]u64 = undefined;
        var product: u128 = undefined;

        // Tangent keeps one more 64-bit limb because it needs relative reduction
        // accuracy near zeros of both sine and cosine.
        if (exponent <= 1074) {
            if (tangent_precision) {
                const tail_product = @as(u128, significand) * data.t[2];
                product = @as(u128, significand) * data.t[1];
                const product_low: u64 = @truncate(product);
                const tail_high: u64 = @truncate(tail_product >> 64);
                limbs[0] = product_low +% tail_high;
                limbs[1] = @as(u64, @truncate(product >> 64)) +%
                    @intFromBool(limbs[0] < product_low);
            } else {
                product = @as(u128, significand) * data.t[1];
                limbs[0] = @truncate(product);
                limbs[1] = @truncate(product >> 64);
            }
            product = @as(u128, significand) * data.t[0];
            const addend: u64 = @truncate(product);
            const carry = limbs[1] +% addend < limbs[1];
            limbs[1] +%= addend;
            limbs[2] = @as(u64, @truncate(product >> 64)) +% @intFromBool(carry);
            exponent = 1075 - exponent;
        } else {
            const table_index: usize = @intCast(@divTrunc(exponent - 1138 + 63, 64));
            if (tangent_precision) {
                const tail_product = @as(u128, significand) * data.t[table_index + 3];
                product = @as(u128, significand) * data.t[table_index + 2];
                const product_low: u64 = @truncate(product);
                const tail_high: u64 = @truncate(tail_product >> 64);
                limbs[0] = product_low +% tail_high;
                limbs[1] = @as(u64, @truncate(product >> 64)) +%
                    @intFromBool(limbs[0] < product_low);
            } else {
                product = @as(u128, significand) * data.t[table_index + 2];
                limbs[0] = @truncate(product);
                limbs[1] = @truncate(product >> 64);
            }
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
        // The large reducer error is in turns; convert it to an absolute
        // sine/cosine output bound with 2*pi < 8.
        if (!tangent_precision) error_bound.* = 0x1.01p-73;
    }

    const index = reductionIndex(high.*);
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

const FastTail = struct {
    linear: DoubleDouble,
    sine: f64,
    cosine: f64,
};

inline fn evalTailFast(high: f64, low: f64) FastTail {
    @setFloatMode(.strict);
    // |high| < 2^-11 + 2^-24 and |low| < 2^-49. Retain the linear
    // sine term as a pair and include both first-order low corrections.
    const square = high * high;
    const cubic = square * high;
    const cross = high * low;
    var sine = cubic * @mulAdd(f64, 0x1.466bc6775aae2p+6, square, -0x1.4abbce625be53p+5);
    sine = @mulAdd(f64, -0x1.f019b59389d7cp+6, high * cross, sine);
    var cosine = square * @mulAdd(
        f64,
        @mulAdd(f64, -0x1.55d3c7e3cbffap+6, square, 0x1.03c1f081b5ac4p+6),
        square,
        -0x1.3bd3cc9be45dep+4,
    );
    cosine = @mulAdd(f64, -0x1.3bd3cc9be45dep+5, cross, cosine);
    var linear = aMul(0x1.921fb54442d18p+2, high);
    linear.low = @mulAdd(f64, 0x1.1a62633145c07p-52, high, linear.low);
    linear.low = @mulAdd(f64, 0x1.921fb54442d18p+2, low, linear.low);
    return .{ .linear = linear, .sine = sine, .cosine = cosine };
}

inline fn combineTailFast(
    cosine_coefficient: f64,
    sine_coefficient: f64,
    tail: FastTail,
    negative: bool,
    reduction_error: f64,
) FastResult {
    @setFloatMode(.strict);
    const scaled = aMul(sine_coefficient, tail.linear.high);
    var correction = @mulAdd(f64, sine_coefficient, tail.linear.low + tail.sine, scaled.low);
    correction = @mulAdd(f64, cosine_coefficient, tail.cosine, correction);
    const sum = fastTwoSum(cosine_coefficient, scaled.high);
    // The degree-five sine and degree-six cosine tails, coefficient
    // errors, table errors and compensated evaluation fit within 2^-65.
    return .{
        .high = negateIf(sum.high, negative),
        .low = negateIf(sum.low + correction, negative),
        .error_bound = 0x1p-65 + reduction_error,
    };
}

fn sinFast(absolute_value: f64, input_negative: bool) FastResult {
    @setFloatMode(.strict);
    var negative = input_negative;
    var use_sine = true;

    var high: f64 = undefined;
    var low: f64 = undefined;
    var reduction_error: f64 = undefined;
    var index = reduceFast(false, &high, &low, absolute_value, &reduction_error);

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

    const table_index: u8 = @truncate(index);

    high -= data.sc[table_index][0];
    const tail = evalTailFast(high, low);
    if (use_sine) return combineTailFast(
        data.sc[table_index][1],
        data.sc[table_index][2],
        tail,
        negative,
        reduction_error,
    );
    return combineTailFast(data.sc[table_index][2], -data.sc[table_index][1], tail, negative, reduction_error);
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

noinline fn sinAccurate(value: f64) f64 {
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
    var index = reduceFast(false, &high, &low, value, &reduction_error);

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

    const table_index: u8 = @truncate(index);

    high -= data.sc[table_index][0];
    const tail = evalTailFast(high, low);
    if (!use_cosine) return combineTailFast(
        data.sc[table_index][1],
        data.sc[table_index][2],
        tail,
        negative,
        reduction_error,
    );
    return combineTailFast(data.sc[table_index][2], -data.sc[table_index][1], tail, negative, reduction_error);
}

noinline fn cosAccurate(value: f64) f64 {
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
    if (magnitude <= 0x3e46_a09e_667f_3bcc) return 1.0;
    if (magnitude <= 0x3fc0_0000_0000_0000) return cosSmall(value);

    const absolute_value: f64 = @bitCast(magnitude);
    const fast = cosFast(absolute_value);
    const left = fast.high + (fast.low - fast.error_bound);
    const right = fast.high + (fast.low + fast.error_bound);
    if (left != right) {
        @branchHint(.cold);
        return cosAccurate(absolute_value);
    }
    return left;
}

inline fn cosSmall(value: f64) f64 {
    @setFloatMode(.strict);
    // Above the one-rounding cutoff, x^2 >= u. On |x| <= 1/8, the
    // degree-twelve tail and evaluation error fit within 2u*x^4.
    // Include rounding of the retained low sum separately (u=2^-53).
    const square = aMul(value, value);
    const fourth = square.high * square.high;
    var polynomial: f64 = 0x1.1eed8eff8d898p-29;
    polynomial = @mulAdd(f64, polynomial, square.high, -0x1.27e4fb7789f5cp-22);
    polynomial = @mulAdd(f64, polynomial, square.high, 0x1.a01a01a01a01ap-16);
    polynomial = @mulAdd(f64, polynomial, square.high, -0x1.6c16c16c16c17p-10);
    polynomial = @mulAdd(f64, polynomial, square.high, 0x1.5555555555555p-5);
    const correction = @mulAdd(f64, fourth, polynomial, -0.5 * square.low);
    const sum = fastTwoSum(1.0, -0.5 * square.high);
    const low = sum.low + correction;
    const error_bound = (@abs(sum.low) + fourth) * 0x1p-52;
    const left = sum.high + (low - error_bound);
    const right = sum.high + (low + error_bound);
    if (left == right) return left;
    return cosAccurate(@abs(value));
}

fn tangentRatio(
    index: u8,
    high: f64,
    low: f64,
) struct { numerator: DoubleDouble, denominator: DoubleDouble } {
    @setFloatMode(.strict);
    // Only the nonlinear term uses the rounded argument. Its evaluation
    // error is bounded by 16*u times that term, u=2^-53. The paired
    // linear term retains the reduction tail and the low part of 2*pi.
    const argument = high + low;
    const square = argument * argument;
    const cubic = square * argument;
    var polynomial = @mulAdd(f64, 0x1.45fff9b48e95ep+14, square, 0x1.466bc6775aae2p+10);
    polynomial = @mulAdd(f64, polynomial, square, 0x1.4abbce625be53p+6);
    const nonlinear = cubic * polynomial;
    var linear = aMul(0x1.921fb54442d18p+2, high);
    linear.low = @mulAdd(f64, 0x1.1a62633145c07p-52, high, linear.low);
    linear.low = @mulAdd(f64, 0x1.921fb54442d18p+2, low, linear.low);
    const tail = fastTwoSum(linear.high, linear.low + nonlinear);

    // tan(a+b) = (tan(a)+tan(b))/(1-tan(a)*tan(b)).
    // Both sums are well-conditioned on this reduced quadrant.
    const row = data.tan_fast[index];
    // Row zero is exact zero. Every other row is at least tan(2*pi/2048);
    // the tail has no larger exponent. FastTwoSum is exact in both cases.
    const numerator = fastTwoSum(row[0], tail.high);
    const cross = dMul(row[0], row[1], tail.high, tail.low);
    const denominator = fastTwoSum(1.0, -cross.high);
    return .{
        .numerator = .{ .high = numerator.high, .low = numerator.low + (row[1] + tail.low) },
        .denominator = .{ .high = denominator.high, .low = denominator.low - cross.low },
    };
}

fn fastDiv(
    numerator_high: f64,
    numerator_low: f64,
    denominator_high: f64,
    denominator_low: f64,
) DoubleDouble {
    @setFloatMode(.strict);
    const reciprocal = 1.0 / denominator_high;
    const high = numerator_high * reciprocal;
    const residual_high = @mulAdd(f64, denominator_high, -high, numerator_high);
    const residual_low = @mulAdd(f64, denominator_low, -high, numerator_low);
    return .{
        .high = high,
        .low = reciprocal * (residual_high + residual_low),
    };
}

pub fn tanFast(high: *f64, low: *f64, value: f64) f64 {
    @setFloatMode(.strict);
    const input_bits: u64 = @bitCast(value);
    var negative = (input_bits & sign_mask) != 0;
    const absolute_value: f64 = @bitCast(input_bits & magnitude_mask);
    var is_tan = true;
    var reduction_error_bound: f64 = 0.0;
    var index = reduceFast(
        true,
        high,
        low,
        absolute_value,
        &reduction_error_bound,
    );

    index &= 0x3ff;
    is_tan = is_tan != ((index >> 9) != 0);
    negative = negative != ((index >> 9) != 0);
    index &= 0x1ff;

    const reflected = (index & 0x100) != 0;
    is_tan = is_tan != reflected;
    const reflection_mask = @as(u16, 0) -% @as(u16, @intFromBool(reflected));
    index = (index ^ reflection_mask) & 0xff;
    const offset_bits = @as(u64, 0x3f40_0000_0000_0000) & (@as(u64, 0) -% @as(u64, @intFromBool(reflected)));
    high.* = @as(f64, @bitCast(offset_bits)) + negateIf(high.*, reflected);
    low.* = negateIf(low.*, reflected);

    const table_index: u8 = @truncate(index);
    // The medium reducer uses a conservative 2^-80 absolute turn bound.
    // Excluding a 2^-15 neighborhood keeps its tangent contribution below
    // 2^-63 relative; the older reducers retain their tighter guard.
    const minimum_remainder: f64 = if (reduction_error_bound == 0.0)
        0x1p-37
    else
        0x1p-15;
    if (index == 0 and high.* < minimum_remainder) {
        @branchHint(.cold);
        return 1.0;
    }

    const ratio = tangentRatio(table_index, high.*, low.*);
    const quotient = if (is_tan)
        fastDiv(ratio.numerator.high, ratio.numerator.low, ratio.denominator.high, ratio.denominator.low)
    else
        fastDiv(ratio.denominator.high, ratio.denominator.low, ratio.numerator.high, ratio.numerator.low);
    high.* = negateIf(quotient.high, negative);
    low.* = negateIf(quotient.low, negative);

    // The tangent polynomial, table pairs, and compensated quotient fit
    // within 2^-65 relative. Fixed-fraction reconstruction adds less than
    // 2^-69; include reflection and low-part slack, and the medium bound.
    const relative_error_bound: f64 = if (reduction_error_bound == 0.0)
        0x1.2p-65
    else
        0x1p-62;
    return @abs(high.*) * relative_error_bound;
}

// The following Dint quotient and finite tangent path are derived from
// CORE-MATH src/binary64/tan/tan.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa (MIT License).

fn divideDint(result: *Dint, numerator: *const Dint, denominator: *const Dint) void {
    const numerator_bits = dintBits(numerator.*);
    const denominator_bits = dintBits(denominator.*);
    const integer_part = numerator_bits / denominator_bits;
    var remainder = numerator_bits % denominator_bits;
    var quotient = integer_part << 127;
    var bit_index: i32 = 126;
    const half_denominator = (denominator_bits >> 1) + (denominator_bits & 1);
    while (bit_index >= 0) : (bit_index -= 1) {
        if (remainder >= half_denominator) {
            remainder = (remainder << 1) -% denominator_bits;
            quotient |= @as(u128, 1) << @as(u7, @intCast(bit_index));
        } else {
            remainder <<= 1;
        }
    }
    const round_up = if ((denominator_bits & 1) == 0)
        remainder > half_denominator or
            (remainder == half_denominator and (quotient & 1) != 0)
    else
        remainder >= half_denominator;
    const quotient_overflow = round_up and quotient == ~@as(u128, 0);
    if (round_up and !quotient_overflow) quotient += 1;

    var exponent = numerator.ex - denominator.ex;
    if (quotient_overflow) {
        quotient = @as(u128, 1) << 127;
        exponent += 2;
    } else if ((quotient >> 127) != 0) {
        exponent += 1;
    } else {
        quotient <<= 1;
    }
    if (quotient == 0) {
        result.* = zero;
        result.sgn = numerator.sgn ^ denominator.sgn;
        return;
    }
    result.* = .{
        .hi = @truncate(quotient >> 64),
        .lo = @truncate(quotient),
        .ex = exponent,
        .sgn = numerator.sgn ^ denominator.sgn,
    };
}

/// Return the direct CORE-MATH accurate tangent path.
pub noinline fn tanAccurate(value: f64) f64 {
    @setFloatMode(.strict);
    const input_bits: u64 = @bitCast(value);
    const absolute_bits = input_bits & magnitude_mask;
    const absolute_value: f64 = @bitCast(absolute_bits);
    var reduced = dintFromDouble(absolute_value);
    reduce(&reduced);

    var is_tan = true;
    var negative = (input_bits & sign_mask) != 0;
    var index = reduce2(&reduced);
    index &= 0x3ff;

    if ((index & 0x200) != 0) {
        is_tan = false;
        negative = !negative;
        index &= 0x1ff;
    }
    if ((index & 0x100) != 0) {
        is_tan = !is_tan;
        reduced.sgn = 1;
        addDint(&reduced, &magic, &reduced);
        index = 0x1ff - index;
    }

    var square: Dint = undefined;
    mulDint(&square, &reduced, &reduced);
    const cosine_tail = evalCosAccurate(&square);
    const sine_tail = evalSinAccurate(&reduced, &square);
    var sine_result: Dint = undefined;
    var cosine_result: Dint = undefined;
    var sine_from_cosine: Dint = undefined;
    var cosine_from_sine: Dint = undefined;
    mulDint(&sine_from_cosine, &data.s[index], &cosine_tail);
    mulDint(&cosine_from_sine, &data.c[index], &sine_tail);
    addDint(&sine_result, &sine_from_cosine, &cosine_from_sine);
    mulDint(&cosine_result, &data.c[index], &cosine_tail);
    mulDint(&cosine_from_sine, &data.s[index], &sine_tail);
    cosine_from_sine.sgn = 1 - cosine_from_sine.sgn;
    addDint(&cosine_result, &cosine_result, &cosine_from_sine);
    var numerator: Dint = undefined;
    if (is_tan) {
        divideDint(&numerator, &sine_result, &cosine_result);
    } else {
        divideDint(&numerator, &cosine_result, &sine_result);
    }
    if (negative) numerator.sgn = 1 - numerator.sgn;
    return dintToDouble(numerator);
}

test "fast cosine encloses a wrapped full-turn remainder" {
    const testing = @import("std").testing;
    const result = cosFast(0x1.921fb54442d18p+25);
    // MPFR-512 rounded to binary128; its rounding error is below 2^-113.
    const reference: f128 = 0xf.fffffffffffffd91062f4eb457p-4;
    try testing.expectApproxEqAbs(
        reference,
        @as(f128, result.high) + @as(f128, result.low),
        @as(f128, result.error_bound),
    );
}

test "sine preserves rounding after a reflected medium reduction" {
    const testing = @import("std").testing;
    const input: f64 = @bitCast(@as(u64, 0x413e_a544_0398_0ced));
    try testing.expectEqual(
        @as(u64, 0xbfe5_e4c7_0940_cf5e),
        @as(u64, @bitCast(sind(input))),
    );
}
