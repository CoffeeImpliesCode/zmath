const std = @import("std");
const common = @import("scalar_c23_common.zig");

/// The bounded C23 power families use at most 68,800 exact bits for a
/// binary64 minimum subnormal raised to the 64th power.  The extra capacity is
/// deliberate: it also covers one-bit rounding shifts and the compoundn
/// expansion of `1 + x` before exponentiation.
pub const limb_count: usize = 2048;
pub const capacity_bits: usize = limb_count * 64;
pub const max_degree: i64 = 64;
pub const max_exact_bits: usize = 68_800;

comptime {
    if (capacity_bits < 131_072) @compileError("exact power limb capacity is too small");
    if (max_exact_bits + 2 > capacity_bits) @compileError("exact power bound exceeds limb capacity");
}

fn Arithmetic(comptime storage_limbs: usize) type {
    return struct {
        const storage_bits = storage_limbs * 64;

        const Big = struct {
            limbs: [storage_limbs]u64,
            len: usize,

            fn zero() Big {
                return .{
                    .limbs = [_]u64{0} ** storage_limbs,
                    .len = 0,
                };
            }

            fn one() Big {
                var result = zero();
                result.limbs[0] = 1;
                result.len = 1;
                return result;
            }

            fn fromU64(value: u64) Big {
                var result = zero();
                if (value != 0) {
                    result.limbs[0] = value;
                    result.len = 1;
                }
                return result;
            }

            fn normalize(self: *Big) void {
                while (self.len != 0 and self.limbs[self.len - 1] == 0) {
                    self.len -= 1;
                }
            }

            fn isZero(self: *const Big) bool {
                return self.len == 0;
            }

            fn bitLength(self: *const Big) usize {
                if (self.len == 0) return 0;
                const top = self.limbs[self.len - 1];
                return (self.len - 1) * 64 + (64 - @as(usize, @clz(top)));
            }

            fn setBit(self: *Big, bit_index: usize) void {
                std.debug.assert(bit_index < storage_bits);
                const limb_index = bit_index / 64;
                const bit_in_limb: u6 = @intCast(bit_index & 63);
                self.limbs[limb_index] |= @as(u64, 1) << bit_in_limb;
                if (self.len <= limb_index) self.len = limb_index + 1;
            }

            fn lowU64(self: *const Big) u64 {
                if (self.len == 0) return 0;
                return self.limbs[0];
            }

            fn compare(left: *const Big, right: *const Big) i8 {
                if (left.len < right.len) return -1;
                if (left.len > right.len) return 1;
                var index = left.len;
                while (index != 0) {
                    index -= 1;
                    if (left.limbs[index] < right.limbs[index]) return -1;
                    if (left.limbs[index] > right.limbs[index]) return 1;
                }
                return 0;
            }

            fn copy(source: *const Big, destination: *Big) void {
                destination.* = source.*;
            }

            /// Computes `source * 2^shift` into a distinct fixed-limb value.
            fn shifted(source: *const Big, shift: usize, destination: *Big) void {
                destination.* = zero();
                if (source.len == 0) return;

                const whole_limbs = shift / 64;
                const bit_shift: usize = shift & 63;
                const shift_amount: u6 = @intCast(bit_shift);
                const reverse_shift: u6 = @intCast(if (bit_shift == 0) 0 else 64 - bit_shift);
                std.debug.assert(source.len + whole_limbs <= storage_limbs);
                if (bit_shift != 0) std.debug.assert(source.len + whole_limbs < storage_limbs);

                var index: usize = 0;
                while (index < source.len) : (index += 1) {
                    const destination_index = index + whole_limbs;
                    const limb = source.limbs[index];
                    destination.limbs[destination_index] |= limb << shift_amount;
                    if (bit_shift != 0) {
                        destination.limbs[destination_index + 1] |= limb >> reverse_shift;
                    }
                }
                destination.len = source.len + whole_limbs;
                if (bit_shift != 0 and source.limbs[source.len - 1] >> reverse_shift != 0) {
                    destination.len += 1;
                }
                destination.normalize();
            }

            fn add(left: *const Big, right: *const Big, destination: *Big) void {
                destination.* = zero();
                const count = if (left.len > right.len) left.len else right.len;
                std.debug.assert(count <= storage_limbs);

                var carry: u64 = 0;
                var index: usize = 0;
                while (index < count) : (index += 1) {
                    const left_limb = if (index < left.len) left.limbs[index] else 0;
                    const right_limb = if (index < right.len) right.limbs[index] else 0;
                    const first = @addWithOverflow(left_limb, right_limb);
                    const second = @addWithOverflow(first[0], carry);
                    destination.limbs[index] = second[0];
                    carry = if (first[1] != 0 or second[1] != 0) 1 else 0;
                }
                destination.len = count;
                if (carry != 0) {
                    std.debug.assert(destination.len < storage_limbs);
                    destination.limbs[destination.len] = carry;
                    destination.len += 1;
                }
            }

            fn decrement(self: *Big) void {
                std.debug.assert(self.len != 0);
                var index: usize = 0;
                while (index < self.len) : (index += 1) {
                    const result = @subWithOverflow(self.limbs[index], 1);
                    self.limbs[index] = result[0];
                    if (result[1] == 0) break;
                }
                self.normalize();
            }

            fn increment(self: *Big) void {
                var index: usize = 0;
                while (index < self.len) : (index += 1) {
                    const result = @addWithOverflow(self.limbs[index], 1);
                    self.limbs[index] = result[0];
                    if (result[1] == 0) return;
                }
                std.debug.assert(self.len < storage_limbs);
                self.limbs[self.len] = 1;
                self.len += 1;
            }

            fn subtract(left: *const Big, right: *const Big, destination: *Big) void {
                std.debug.assert(compare(left, right) >= 0);
                destination.* = zero();

                var borrow: u64 = 0;
                var index: usize = 0;
                while (index < left.len) : (index += 1) {
                    const left_limb = left.limbs[index];
                    const right_limb = if (index < right.len) right.limbs[index] else 0;
                    const first = @subWithOverflow(left_limb, right_limb);
                    const second = @subWithOverflow(first[0], borrow);
                    destination.limbs[index] = second[0];
                    borrow = if (first[1] != 0 or second[1] != 0) 1 else 0;
                }
                std.debug.assert(borrow == 0);
                destination.len = left.len;
                destination.normalize();
            }

            /// Fixed-limb schoolbook multiplication.  The `u128` accumulator is only
            /// an intra-limb carry; all persistent state remains the stack `u64` array.
            fn multiply(left: *const Big, right: *const Big, destination: *Big) void {
                destination.* = zero();
                if (left.len == 0 or right.len == 0) return;
                std.debug.assert(left.len + right.len <= storage_limbs);

                var left_index: usize = 0;
                while (left_index < left.len) : (left_index += 1) {
                    var carry: u128 = 0;
                    var right_index: usize = 0;
                    while (right_index < right.len) : (right_index += 1) {
                        const destination_index = left_index + right_index;
                        const product = @as(u128, left.limbs[left_index]) *
                            @as(u128, right.limbs[right_index]);
                        const total = product + @as(u128, destination.limbs[destination_index]) + carry;
                        destination.limbs[destination_index] = @truncate(total);
                        carry = total >> 64;
                    }

                    var destination_index = left_index + right.len;
                    while (carry != 0) {
                        std.debug.assert(destination_index < storage_limbs);
                        const total = @as(u128, destination.limbs[destination_index]) + carry;
                        destination.limbs[destination_index] = @truncate(total);
                        carry = total >> 64;
                        destination_index += 1;
                    }
                }
                destination.len = left.len + right.len;
                destination.normalize();
            }

            fn power(base: *const Big, exponent: u64, destination: *Big) void {
                std.debug.assert(exponent <= @as(u64, @intCast(max_degree)));
                if (exponent != 0) {
                    std.debug.assert(base.bitLength() <= storage_bits / @as(usize, @intCast(exponent)));
                }
                var result = one();
                // SAFETY: copy initializes factor before its first read.
                var factor: Big = undefined;
                copy(base, &factor);
                var remaining = exponent;
                while (remaining != 0) {
                    if ((remaining & 1) != 0) {
                        // SAFETY: multiply initializes product before it is copied.
                        var product: Big = undefined;
                        multiply(&result, &factor, &product);
                        result = product;
                    }
                    remaining >>= 1;
                    if (remaining != 0) {
                        // SAFETY: multiply initializes square before it is copied.
                        var square: Big = undefined;
                        multiply(&factor, &factor, &square);
                        factor = square;
                    }
                }
                destination.* = result;
            }

            /// Divides `numerator` by nonzero `denominator`, retaining the exact
            /// remainder. Work scales with quotient bits, not the numerator's leading bits.
            fn divideRemainder(
                numerator: *const Big,
                denominator: *const Big,
                quotient: *Big,
                remainder: *Big,
            ) void {
                std.debug.assert(denominator.len != 0);
                quotient.* = zero();
                copy(numerator, remainder);
                const numerator_bits = numerator.bitLength();
                const denominator_bits = denominator.bitLength();
                if (numerator_bits < denominator_bits) return;

                var bit_index = numerator_bits - denominator_bits + 1;
                while (bit_index != 0) {
                    bit_index -= 1;
                    var divisor: Big = undefined;
                    shifted(denominator, bit_index, &divisor);
                    if (compare(remainder, &divisor) >= 0) {
                        var difference: Big = undefined;
                        subtract(remainder, &divisor, &difference);
                        copy(&difference, remainder);
                        quotient.setBit(bit_index);
                    }
                }
                std.debug.assert(compare(remainder, denominator) < 0);
            }
        };

        const FormatInfo = struct {
            fraction_bits: u8,
            precision_bits: u8,
            exponent_bias: i64,
            minimum_exponent: i64,
            maximum_exponent: i64,
            exponent_field_mask: u64,
            fraction_mask: u64,
            sign_mask: u64,
            infinity_bits: u64,
            maximum_finite_bits: u64,
        };

        fn formatInfo(comptime Format: type) FormatInfo {
            if (Format == common.H) {
                return .{
                    .fraction_bits = 10,
                    .precision_bits = 11,
                    .exponent_bias = 15,
                    .minimum_exponent = -14,
                    .maximum_exponent = 15,
                    .exponent_field_mask = 0x7c00,
                    .fraction_mask = 0x03ff,
                    .sign_mask = 0x8000,
                    .infinity_bits = 0x7c00,
                    .maximum_finite_bits = 0x7bff,
                };
            }
            if (Format == common.F) {
                return .{
                    .fraction_bits = 23,
                    .precision_bits = 24,
                    .exponent_bias = 127,
                    .minimum_exponent = -126,
                    .maximum_exponent = 127,
                    .exponent_field_mask = 0x7f80_0000,
                    .fraction_mask = 0x007f_ffff,
                    .sign_mask = 0x8000_0000,
                    .infinity_bits = 0x7f80_0000,
                    .maximum_finite_bits = 0x7f7f_ffff,
                };
            }
            if (Format == common.D) {
                return .{
                    .fraction_bits = 52,
                    .precision_bits = 53,
                    .exponent_bias = 1023,
                    .minimum_exponent = -1022,
                    .maximum_exponent = 1023,
                    .exponent_field_mask = 0x7ff0_0000_0000_0000,
                    .fraction_mask = 0x000f_ffff_ffff_ffff,
                    .sign_mask = 0x8000_0000_0000_0000,
                    .infinity_bits = 0x7ff0_0000_0000_0000,
                    .maximum_finite_bits = 0x7fef_ffff_ffff_ffff,
                };
            }
            @compileError("unsupported exact power format");
        }

        const Decoded = struct {
            significand: Big,
            exponent: i64,
        };

        fn decodeFinite(comptime Format: type, bits: Format.Bits) Decoded {
            const info = formatInfo(Format);
            const fraction_shift: u6 = @intCast(info.fraction_bits);
            const raw_bits: u64 = @intCast(bits);
            const exponent_field = (raw_bits & info.exponent_field_mask) >> fraction_shift;
            const fraction = raw_bits & info.fraction_mask;

            if (exponent_field == 0) {
                return .{
                    .significand = Big.fromU64(fraction),
                    .exponent = info.minimum_exponent - @as(i64, info.fraction_bits),
                };
            }

            const hidden_bit = @as(u64, 1) << fraction_shift;
            return .{
                .significand = Big.fromU64(hidden_bit | fraction),
                .exponent = @as(i64, @intCast(exponent_field)) - info.exponent_bias -
                    @as(i64, info.fraction_bits),
            };
        }

        fn magnitude(exponent: i64) u64 {
            const bits: u64 = @bitCast(exponent);
            if (exponent < 0) return 0 -% bits;
            return bits;
        }

        fn compareScaled(
            left: *const Big,
            left_shift: i64,
            right: *const Big,
            right_shift: i64,
        ) i8 {
            if (left.len == 0 and right.len == 0) return 0;
            if (left.len == 0) return -1;
            if (right.len == 0) return 1;

            const left_top = @as(i64, @intCast(left.bitLength())) + left_shift;
            const right_top = @as(i64, @intCast(right.bitLength())) + right_shift;
            if (left_top < right_top) return -1;
            if (left_top > right_top) return 1;

            if (left_shift > right_shift) {
                // SAFETY: `Big.shifted` overwrites all of `shifted_left` before comparison.
                var shifted_left: Big = undefined;
                Big.shifted(left, @intCast(left_shift - right_shift), &shifted_left);
                return Big.compare(&shifted_left, right);
            }
            if (right_shift > left_shift) {
                // SAFETY: `Big.shifted` overwrites all of `shifted_right` before comparison.
                var shifted_right: Big = undefined;
                Big.shifted(right, @intCast(right_shift - left_shift), &shifted_right);
                return Big.compare(left, &shifted_right);
            }
            return Big.compare(left, right);
        }

        const Rational = struct {
            numerator: Big,
            denominator: Big,
            shift: i64,
        };

        fn rationalFromPower(decoded: *const Decoded, exponent: i64, rational: *Rational) void {
            const degree = magnitude(exponent);
            if (degree != 0) {
                std.debug.assert(decoded.significand.bitLength() <=
                    storage_bits / @as(usize, @intCast(degree)));
            }
            // SAFETY: `Big.power` always overwrites `powered` before the sign branch reads it.
            var powered: Big = undefined;
            Big.power(&decoded.significand, degree, &powered);

            if (exponent >= 0) {
                rational.* = .{
                    .numerator = powered,
                    .denominator = Big.one(),
                    .shift = decoded.exponent * @as(i64, @intCast(degree)),
                };
            } else {
                rational.* = .{
                    .numerator = Big.one(),
                    .denominator = powered,
                    .shift = -(decoded.exponent * @as(i64, @intCast(degree))),
                };
            }
        }

        fn roundScaledQuotient(rational: *const Rational, scale: i64) u64 {
            // SAFETY: Each `total_shift` branch overwrites `numerator` before the next read.
            var numerator: Big = undefined;
            // SAFETY: Each `total_shift` branch overwrites `denominator` before the next read.
            var denominator: Big = undefined;
            const total_shift = rational.shift + scale;
            if (total_shift >= 0) {
                Big.shifted(&rational.numerator, @intCast(total_shift), &numerator);
                Big.copy(&rational.denominator, &denominator);
            } else {
                Big.copy(&rational.numerator, &numerator);
                Big.shifted(&rational.denominator, @intCast(-total_shift), &denominator);
            }
            std.debug.assert(denominator.bitLength() < storage_bits);

            // SAFETY: `Big.divideRemainder` initializes `quotient` before its first read.
            var quotient: Big = undefined;
            // SAFETY: `Big.divideRemainder` initializes `remainder` before its first read.
            var remainder: Big = undefined;
            Big.divideRemainder(&numerator, &denominator, &quotient, &remainder);

            // RN-even proof: `quotient` is floor(value / ulp), and the exact
            // remainder is below, above, or exactly halfway according to 2*r
            // compared with the divisor.  At equality, increment only an odd low bit.
            // SAFETY: `Big.shifted` overwrites `doubled_remainder` before comparison.
            var doubled_remainder: Big = undefined;
            Big.shifted(&remainder, 1, &doubled_remainder);
            const midpoint_relation = Big.compare(&doubled_remainder, &denominator);
            if (midpoint_relation > 0 or
                (midpoint_relation == 0 and (quotient.lowU64() & 1) != 0))
            {
                quotient.increment();
            }
            std.debug.assert(quotient.bitLength() <= 64);
            return quotient.lowU64();
        }

        fn roundRational(comptime Format: type, rational: *const Rational, negative: bool) Format.Float {
            const info = formatInfo(Format);
            const fraction_shift: u6 = @intCast(info.fraction_bits);
            const precision_shift: u6 = @intCast(info.precision_bits);
            const hidden_shift: u6 = @intCast(info.precision_bits - 1);
            std.debug.assert(!rational.numerator.isZero());
            std.debug.assert(!rational.denominator.isZero());

            const rough_exponent = @as(i64, @intCast(rational.numerator.bitLength())) -
                @as(i64, @intCast(rational.denominator.bitLength())) + rational.shift;
            var exponent = rough_exponent;
            if (compareScaled(
                &rational.numerator,
                rational.shift,
                &rational.denominator,
                exponent,
            ) < 0) exponent -= 1;

            var magnitude_bits: u64 = 0;
            if (exponent < info.minimum_exponent) {
                const scale = @as(i64, info.precision_bits - 1) - info.minimum_exponent;
                const subnormal_significand = roundScaledQuotient(rational, scale);
                const hidden_bit = @as(u64, 1) << hidden_shift;
                if (subnormal_significand == 0) {
                    magnitude_bits = 0;
                } else if (subnormal_significand >= hidden_bit) {
                    magnitude_bits = @as(u64, 1) << fraction_shift;
                } else {
                    magnitude_bits = subnormal_significand;
                }
            } else {
                var significand = roundScaledQuotient(
                    rational,
                    @as(i64, info.precision_bits - 1) - exponent,
                );
                const carry_bit = @as(u64, 1) << precision_shift;
                const hidden_bit = @as(u64, 1) << hidden_shift;
                if (significand >= carry_bit) {
                    significand >>= 1;
                    exponent += 1;
                }

                if (exponent > info.maximum_exponent) {
                    magnitude_bits = info.infinity_bits;
                } else {
                    std.debug.assert(significand >= hidden_bit);
                    const exponent_field: u64 = @intCast(exponent + info.exponent_bias);
                    const fraction = significand - hidden_bit;
                    magnitude_bits = (exponent_field << fraction_shift) | fraction;
                }
            }

            if (negative) magnitude_bits |= info.sign_mask;
            return Format.fromBits(@intCast(magnitude_bits));
        }

        fn makeCompoundBase(comptime Format: type, value_bits: Format.Bits) Decoded {
            const info = formatInfo(Format);
            const raw_bits: u64 = @intCast(value_bits);
            const negative = (raw_bits & info.sign_mask) != 0;
            const decoded = decodeFinite(Format, value_bits);
            if (decoded.significand.isZero()) return .{ .significand = Big.one(), .exponent = 0 };

            if (decoded.exponent >= 0) {
                // SAFETY: `Big.shifted` overwrites `integer_part` before `Big.add` reads it.
                var integer_part: Big = undefined;
                Big.shifted(&decoded.significand, @intCast(decoded.exponent), &integer_part);
                const one = Big.one();
                // SAFETY: `Big.add` overwrites all of `result` before the returned copy reads it.
                var result: Big = undefined;
                Big.add(&integer_part, &one, &result);
                return .{ .significand = result, .exponent = 0 };
            }

            const denominator_shift: usize = @intCast(-decoded.exponent);
            var unit: Big = Big.zero();
            unit.setBit(denominator_shift);
            // SAFETY: Both `negative` branches overwrite `result` before the returned copy reads it.
            var result: Big = undefined;
            if (negative) {
                std.debug.assert(Big.compare(&unit, &decoded.significand) > 0);
                Big.subtract(&unit, &decoded.significand, &result);
            } else {
                Big.add(&unit, &decoded.significand, &result);
            }
            return .{ .significand = result, .exponent = decoded.exponent };
        }

        fn midpoint(lower_bits: u64, upper_bits: u64, info: FormatInfo, result: *Decoded) void {
            const lower_decoded = decodeFiniteBits(lower_bits, info);
            const upper_decoded = decodeFiniteBits(upper_bits, info);

            if (lower_decoded.significand.isZero()) {
                result.* = .{
                    .significand = upper_decoded.significand,
                    .exponent = upper_decoded.exponent - 1,
                };
                return;
            }

            const minimum_exponent = if (lower_decoded.exponent < upper_decoded.exponent)
                lower_decoded.exponent
            else
                upper_decoded.exponent;
            // SAFETY: The first `Big.shifted` call overwrites `lower_scaled` before `Big.add`.
            var lower_scaled: Big = undefined;
            // SAFETY: The second `Big.shifted` call overwrites `upper_scaled` before `Big.add`.
            var upper_scaled: Big = undefined;
            Big.shifted(
                &lower_decoded.significand,
                @intCast(lower_decoded.exponent - minimum_exponent),
                &lower_scaled,
            );
            Big.shifted(
                &upper_decoded.significand,
                @intCast(upper_decoded.exponent - minimum_exponent),
                &upper_scaled,
            );
            // SAFETY: `Big.add` overwrites all of `sum` before the result copies it.
            var sum: Big = undefined;
            Big.add(&lower_scaled, &upper_scaled, &sum);
            result.* = .{ .significand = sum, .exponent = minimum_exponent - 1 };
        }

        fn decodeFiniteBits(raw_bits: u64, info: FormatInfo) Decoded {
            const fraction_shift: u6 = @intCast(info.fraction_bits);
            const exponent_field = (raw_bits & info.exponent_field_mask) >> fraction_shift;
            const fraction = raw_bits & info.fraction_mask;
            if (exponent_field == 0) {
                return .{
                    .significand = Big.fromU64(fraction),
                    .exponent = info.minimum_exponent - @as(i64, info.fraction_bits),
                };
            }
            return .{
                .significand = Big.fromU64((@as(u64, 1) << fraction_shift) | fraction),
                .exponent = @as(i64, @intCast(exponent_field)) - info.exponent_bias -
                    @as(i64, info.fraction_bits),
            };
        }

        fn candidateRelation(
            candidate_bits: u64,
            input: *const Decoded,
            degree: u64,
            reciprocal: bool,
            info: FormatInfo,
        ) i8 {
            const candidate = decodeFiniteBits(candidate_bits, info);
            if (candidate.significand.isZero()) return -1;

            // SAFETY: `Big.power` overwrites `candidate_power` before either relation path reads it.
            var candidate_power: Big = undefined;
            Big.power(&candidate.significand, degree, &candidate_power);
            const candidate_shift = candidate.exponent * @as(i64, @intCast(degree));
            if (!reciprocal) {
                return compareScaled(
                    &candidate_power,
                    candidate_shift,
                    &input.significand,
                    input.exponent,
                );
            }
            // SAFETY: `Big.multiply` overwrites all of `product` before comparison.
            var product: Big = undefined;
            Big.multiply(&candidate_power, &input.significand, &product);
            const one = Big.one();
            return compareScaled(&product, candidate_shift + input.exponent, &one, 0);
        }

        fn midpointRelation(
            lower_bits: u64,
            upper_bits: u64,
            input: *const Decoded,
            degree: u64,
            reciprocal: bool,
            info: FormatInfo,
        ) i8 {
            // SAFETY: `midpoint` overwrites `middle` on every path before the next statement reads it.
            var middle: Decoded = undefined;
            midpoint(lower_bits, upper_bits, info, &middle);
            // SAFETY: `Big.power` overwrites `middle_power` before either relation path reads it.
            var middle_power: Big = undefined;
            Big.power(&middle.significand, degree, &middle_power);
            const middle_shift = middle.exponent * @as(i64, @intCast(degree));
            if (!reciprocal) {
                return compareScaled(&middle_power, middle_shift, &input.significand, input.exponent);
            }
            // SAFETY: `Big.multiply` overwrites all of `product` before comparison.
            var product: Big = undefined;
            Big.multiply(&middle_power, &input.significand, &product);
            const one = Big.one();
            return compareScaled(&product, middle_shift + input.exponent, &one, 0);
        }

        fn overflowThresholdRelation(input: *const Decoded, degree: u64, info: FormatInfo) i8 {
            // RN-even overflow begins at max_finite + half its final ulp.  Written
            // exactly, that midpoint is `(2^(p+1)-1) * 2^(emax-p)`.
            var threshold_numerator = Big.zero();
            threshold_numerator.setBit(@intCast(info.precision_bits + 1));
            threshold_numerator.decrement();
            // SAFETY: `Big.power` overwrites `threshold_power` before multiplication.
            var threshold_power: Big = undefined;
            Big.power(&threshold_numerator, degree, &threshold_power);
            const threshold_shift =
                (info.maximum_exponent - @as(i64, info.precision_bits)) * @as(i64, @intCast(degree));

            // SAFETY: `Big.multiply` overwrites all of `product` before comparison.
            var product: Big = undefined;
            Big.multiply(&threshold_power, &input.significand, &product);
            const one = Big.one();
            return compareScaled(&product, threshold_shift + input.exponent, &one, 0);
        }

        fn rootFiniteBits(comptime Format: type, value_bits: Format.Bits, exponent: i64) u64 {
            const info = formatInfo(Format);
            const input = decodeFinite(Format, value_bits);
            std.debug.assert(!input.significand.isZero());
            const degree = magnitude(exponent);
            const reciprocal = exponent < 0;

            var lower: u64 = 0;
            var upper_bound = info.maximum_finite_bits;
            while (lower < upper_bound) {
                const span = upper_bound - lower;
                const middle = lower + span / 2 + (span & 1);
                const relation = candidateRelation(middle, &input, degree, reciprocal, info);
                if (relation <= 0) {
                    lower = middle;
                } else {
                    upper_bound = middle - 1;
                }
            }

            const lower_relation = candidateRelation(lower, &input, degree, reciprocal, info);
            if (lower_relation == 0) return lower;
            if (lower == info.maximum_finite_bits) {
                if (reciprocal and overflowThresholdRelation(&input, degree, info) <= 0) {
                    return info.infinity_bits;
                }
                return lower;
            }

            const upper = lower + 1;
            const middle_relation = midpointRelation(lower, upper, &input, degree, reciprocal, info);
            if (middle_relation < 0) return upper;
            if (middle_relation > 0) return lower;
            // Exact midpoint: the candidate with an even encoded significand wins.
            if ((lower & 1) == 0) return lower;
            return upper;
        }

        /// Exact, one-round pown for a finite positive base and `|exponent| <= 64`.
        fn pownFinite(comptime Format: type, base: Format.Float, exponent: i64) Format.Float {
            std.debug.assert(exponent >= -max_degree);
            std.debug.assert(exponent <= max_degree);
            const decoded = decodeFinite(Format, Format.toBits(base));
            std.debug.assert(!decoded.significand.isZero());
            // SAFETY: `rationalFromPower` overwrites `rational` on both sign paths before return.
            var rational: Rational = undefined;
            rationalFromPower(&decoded, exponent, &rational);
            return roundRational(Format, &rational, false);
        }

        /// Exact, one-round compoundn for finite `x > -1` and `|exponent| <= 64`.
        fn compoundFinite(comptime Format: type, value: Format.Float, exponent: i64) Format.Float {
            std.debug.assert(exponent >= -max_degree);
            std.debug.assert(exponent <= max_degree);
            const base = makeCompoundBase(Format, Format.toBits(value));
            // SAFETY: `rationalFromPower` overwrites `rational` on both sign paths before return.
            var rational: Rational = undefined;
            rationalFromPower(&base, exponent, &rational);
            return roundRational(Format, &rational, false);
        }

        /// Exact candidate-power rootn for finite nonzero positive values and
        /// `|exponent| <= 64`.  Candidate powers and midpoint powers are compared as
        /// integers; no rounded root, reciprocal, or target-format composition is used.
        fn rootFinite(comptime Format: type, value: Format.Float, exponent: i64) Format.Float {
            std.debug.assert(exponent != 0);
            std.debug.assert(exponent >= -max_degree);
            std.debug.assert(exponent <= max_degree);
            const value_bits = Format.toBits(value);
            const result_bits = rootFiniteBits(Format, value_bits, exponent);
            return Format.fromBits(@intCast(result_bits));
        }
    };
}

fn storageLimbs(comptime Format: type) usize {
    const limits: struct { limbs: usize, base_bits: usize, precision: usize } =
        if (Format == common.H)
            .{ .limbs = 32, .base_bits = 25, .precision = 11 }
        else if (Format == common.F)
            .{ .limbs = 256, .base_bits = 150, .precision = 24 }
        else if (Format == common.D)
            .{ .limbs = limb_count, .base_bits = 1075, .precision = 53 }
        else
            @compileError("unsupported exact power format");
    // Exact (1+x)^n needs base_bits*|n| bits, plus final rounding alignment.
    comptime {
        const required = limits.base_bits * @as(usize, @intCast(max_degree)) + limits.precision + 2;
        if (limits.limbs * 64 < required) @compileError("exact power storage bound exceeded");
    }
    return limits.limbs;
}

/// Rounds an integer power of a finite positive base once to the target format.
pub fn pownFinite(comptime Format: type, base: Format.Float, exponent: i64) Format.Float {
    return Arithmetic(storageLimbs(Format)).pownFinite(Format, base, exponent);
}

/// Rounds exact (1+x)^n once for finite x > -1 and |n| <= max_degree.
pub fn compoundFinite(comptime Format: type, value: Format.Float, exponent: i64) Format.Float {
    return Arithmetic(storageLimbs(Format)).compoundFinite(Format, value, exponent);
}

/// Rounds a finite positive root by exact candidate and midpoint comparisons.
pub fn rootFinite(comptime Format: type, value: Format.Float, exponent: i64) Format.Float {
    return Arithmetic(storageLimbs(Format)).rootFinite(Format, value, exponent);
}

test "exact division retains a small quotient and remainder across many limbs" {
    const Big = Arithmetic(limb_count).Big;
    var denominator = Big.one();
    denominator.setBit(1024);
    var numerator = Big.fromU64(49);
    numerator.setBit(1029);
    numerator.setBit(1026);
    numerator.setBit(1024);
    var quotient: Big = undefined;
    var remainder: Big = undefined;
    Big.divideRemainder(&numerator, &denominator, &quotient, &remainder);
    const expected_quotient = Big.fromU64(37);
    const expected_remainder = Big.fromU64(12);
    try std.testing.expectEqual(@as(i8, 0), Big.compare(&quotient, &expected_quotient));
    try std.testing.expectEqual(@as(i8, 0), Big.compare(&remainder, &expected_remainder));
    Big.divideRemainder(&expected_remainder, &denominator, &quotient, &remainder);
    try std.testing.expect(quotient.isZero());
    try std.testing.expectEqual(@as(i8, 0), Big.compare(&remainder, &expected_remainder));
}
