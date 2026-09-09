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
            // Exactly limbs[0..len] are initialized; limbs[len..] are never read.
            limbs: [storage_limbs]u64,
            len: usize,

            fn zero(destination: *Big) void {
                destination.len = 0;
            }

            fn one(destination: *Big) void {
                destination.limbs[0] = 1;
                destination.len = 1;
            }

            fn fromU64(value: u64, destination: *Big) void {
                destination.len = @intFromBool(value != 0);
                if (value != 0) destination.limbs[0] = value;
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
                if (self.len <= limb_index) {
                    @memset(self.limbs[self.len .. limb_index + 1], 0);
                    self.len = limb_index + 1;
                }
                self.limbs[limb_index] |= @as(u64, 1) << bit_in_limb;
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
                std.debug.assert(source != destination);
                @memcpy(destination.limbs[0..source.len], source.limbs[0..source.len]);
                destination.len = source.len;
            }

            /// Computes `source * 2^shift` into a distinct fixed-limb value.
            fn shifted(source: *const Big, shift: usize, destination: *Big) void {
                std.debug.assert(source != destination);
                destination.len = 0;
                if (source.len == 0) return;

                const whole_limbs = shift / 64;
                const bit_shift = shift & 63;
                const output_len = source.len + whole_limbs + @intFromBool(bit_shift != 0);
                std.debug.assert(output_len <= storage_limbs);
                if (storage_limbs != 32) {
                    // The larger formats retain their measured active-prefix shift loop.
                    const shift_amount: u6 = @intCast(bit_shift);
                    const reverse_shift: u6 = @intCast(if (bit_shift == 0) 0 else 64 - bit_shift);
                    @memset(destination.limbs[0..output_len], 0);
                    var index: usize = 0;
                    while (index < source.len) : (index += 1) {
                        const destination_index = index + whole_limbs;
                        const limb = source.limbs[index];
                        destination.limbs[destination_index] |= limb << shift_amount;
                        if (bit_shift != 0) {
                            destination.limbs[destination_index + 1] |= limb >> reverse_shift;
                        }
                    }
                    destination.len = output_len;
                    destination.normalize();
                    return;
                }

                @memset(destination.limbs[0..whole_limbs], 0);

                if (bit_shift == 0) {
                    @memcpy(
                        destination.limbs[whole_limbs..output_len],
                        source.limbs[0..source.len],
                    );
                    destination.len = output_len;
                    return;
                }

                const shift_amount: u6 = @intCast(bit_shift);
                const reverse_shift: u6 = @intCast(64 - bit_shift);
                var carry: u64 = 0;
                for (source.limbs[0..source.len], 0..) |limb, index| {
                    destination.limbs[whole_limbs + index] = (limb << shift_amount) | carry;
                    carry = limb >> reverse_shift;
                }
                destination.len = source.len + whole_limbs;
                if (carry != 0) {
                    destination.limbs[destination.len] = carry;
                    destination.len += 1;
                }
            }

            fn add(left: *const Big, right: *const Big, destination: *Big) void {
                std.debug.assert(left != destination);
                std.debug.assert(right != destination);
                destination.len = 0;
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
                std.debug.assert(left != destination);
                std.debug.assert(right != destination);
                std.debug.assert(compare(left, right) >= 0);
                destination.len = 0;

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
                std.debug.assert(left != destination);
                std.debug.assert(right != destination);
                destination.len = 0;
                if (left.len == 0 or right.len == 0) return;
                std.debug.assert(left.len + right.len <= storage_limbs);

                var left_index: usize = 0;
                if (storage_limbs == 32) {
                    // The first row initializes its digits without reading destination.
                    var carry: u128 = 0;
                    var right_index: usize = 0;
                    while (right_index < right.len) : (right_index += 1) {
                        const total = @as(u128, left.limbs[0]) *
                            @as(u128, right.limbs[right_index]) + carry;
                        destination.limbs[right_index] = @truncate(total);
                        carry = total >> 64;
                    }
                    destination.limbs[right.len] = @intCast(carry);
                    left_index = 1;
                } else {
                    // Larger formats retain the measured active-prefix clearing path.
                    @memset(destination.limbs[0 .. left.len + right.len], 0);
                }

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
                    if (storage_limbs == 32) {
                        destination.limbs[left_index + right.len] = @intCast(carry);
                    } else {
                        var destination_index = left_index + right.len;
                        while (carry != 0) {
                            std.debug.assert(destination_index < storage_limbs);
                            const total = @as(u128, destination.limbs[destination_index]) + carry;
                            destination.limbs[destination_index] = @truncate(total);
                            carry = total >> 64;
                            destination_index += 1;
                        }
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
                var result: Big = undefined;
                one(&result);
                // SAFETY: copy initializes factor before its first read.
                var factor: Big = undefined;
                copy(base, &factor);
                var remaining = exponent;
                while (remaining != 0) {
                    if ((remaining & 1) != 0) {
                        // SAFETY: multiply initializes the active product before copying it.
                        var product: Big = undefined;
                        multiply(&result, &factor, &product);
                        copy(&product, &result);
                    }
                    remaining >>= 1;
                    if (remaining != 0) {
                        // SAFETY: multiply initializes the active square before copying it.
                        var square: Big = undefined;
                        multiply(&factor, &factor, &square);
                        copy(&square, &factor);
                    }
                }
                copy(&result, destination);
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
                quotient.len = 0;
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

        fn decodeFinite(comptime Format: type, bits: Format.Bits, result: *Decoded) void {
            decodeFiniteBits(@intCast(bits), formatInfo(Format), result);
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
                // SAFETY: `Big.shifted` initializes the active limbs before comparison.
                var shifted_left: Big = undefined;
                Big.shifted(left, @intCast(left_shift - right_shift), &shifted_left);
                return Big.compare(&shifted_left, right);
            }
            if (right_shift > left_shift) {
                // SAFETY: `Big.shifted` initializes the active limbs before comparison.
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
            // SAFETY: `Big.power` initializes the active limbs before the sign branch.
            var powered: Big = undefined;
            Big.power(&decoded.significand, degree, &powered);
            const power_shift = decoded.exponent * @as(i64, @intCast(degree));
            if (exponent >= 0) {
                Big.copy(&powered, &rational.numerator);
                rational.denominator.limbs[0] = 1;
                rational.denominator.len = 1;
                rational.shift = power_shift;
            } else {
                rational.numerator.limbs[0] = 1;
                rational.numerator.len = 1;
                Big.copy(&powered, &rational.denominator);
                rational.shift = -power_shift;
            }
        }

        fn roundScaledQuotient(rational: *const Rational, scale: i64) u64 {
            // SAFETY: Both shift branches initialize the active numerator limbs.
            var numerator: Big = undefined;
            // SAFETY: Both shift branches initialize the active denominator limbs.
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
            // SAFETY: `Big.shifted` initializes the active limbs before comparison.
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

        fn makeCompoundBase(
            comptime Format: type,
            value_bits: Format.Bits,
            result: *Decoded,
        ) void {
            const info = formatInfo(Format);
            const raw_bits: u64 = @intCast(value_bits);
            const negative = (raw_bits & info.sign_mask) != 0;
            var decoded: Decoded = undefined;
            decodeFinite(Format, value_bits, &decoded);
            if (decoded.significand.isZero()) {
                result.significand.limbs[0] = 1;
                result.significand.len = 1;
                result.exponent = 0;
                return;
            }

            if (decoded.exponent >= 0) {
                // SAFETY: `Big.shifted` initializes the active integer part.
                var integer_part: Big = undefined;
                Big.shifted(&decoded.significand, @intCast(decoded.exponent), &integer_part);
                var one: Big = undefined;
                Big.one(&one);
                Big.add(&integer_part, &one, &result.significand);
                result.exponent = 0;
                return;
            }

            const denominator_shift: usize = @intCast(-decoded.exponent);
            var unit: Big = undefined;
            Big.zero(&unit);
            unit.setBit(denominator_shift);
            if (negative) {
                std.debug.assert(Big.compare(&unit, &decoded.significand) > 0);
                Big.subtract(&unit, &decoded.significand, &result.significand);
            } else {
                Big.add(&unit, &decoded.significand, &result.significand);
            }
            result.exponent = decoded.exponent;
        }

        fn midpoint(lower_bits: u64, upper_bits: u64, info: FormatInfo, result: *Decoded) void {
            var lower_decoded: Decoded = undefined;
            var upper_decoded: Decoded = undefined;
            decodeFiniteBits(lower_bits, info, &lower_decoded);
            decodeFiniteBits(upper_bits, info, &upper_decoded);

            if (lower_decoded.significand.isZero()) {
                Big.copy(&upper_decoded.significand, &result.significand);
                result.exponent = upper_decoded.exponent - 1;
                return;
            }

            const minimum_exponent = if (lower_decoded.exponent < upper_decoded.exponent)
                lower_decoded.exponent
            else
                upper_decoded.exponent;
            // SAFETY: Each shift initializes its active limbs before the addition.
            var lower_scaled: Big = undefined;
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
            Big.add(&lower_scaled, &upper_scaled, &result.significand);
            result.exponent = minimum_exponent - 1;
        }

        fn decodeFiniteBits(raw_bits: u64, info: FormatInfo, result: *Decoded) void {
            const fraction_shift: u6 = @intCast(info.fraction_bits);
            const exponent_field = (raw_bits & info.exponent_field_mask) >> fraction_shift;
            const fraction = raw_bits & info.fraction_mask;
            const significand = if (exponent_field == 0)
                fraction
            else
                (@as(u64, 1) << fraction_shift) | fraction;
            result.significand.len = @intFromBool(significand != 0);
            if (significand != 0) result.significand.limbs[0] = significand;
            result.exponent = if (exponent_field == 0)
                info.minimum_exponent - @as(i64, info.fraction_bits)
            else
                @as(i64, @intCast(exponent_field)) - info.exponent_bias -
                    @as(i64, info.fraction_bits);
        }

        fn candidateRelation(
            candidate_bits: u64,
            input: *const Decoded,
            degree: u64,
            reciprocal: bool,
            info: FormatInfo,
        ) i8 {
            var candidate: Decoded = undefined;
            decodeFiniteBits(candidate_bits, info, &candidate);
            if (candidate.significand.isZero()) return -1;

            // SAFETY: `Big.power` initializes the active limbs before either relation path.
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
            // SAFETY: `Big.multiply` initializes the active limbs before comparison.
            var product: Big = undefined;
            Big.multiply(&candidate_power, &input.significand, &product);
            var one: Big = undefined;
            Big.one(&one);
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
            // SAFETY: `midpoint` initializes the active limbs and exponent on both paths.
            var middle: Decoded = undefined;
            midpoint(lower_bits, upper_bits, info, &middle);
            // SAFETY: `Big.power` initializes the active limbs before either relation path.
            var middle_power: Big = undefined;
            Big.power(&middle.significand, degree, &middle_power);
            const middle_shift = middle.exponent * @as(i64, @intCast(degree));
            if (!reciprocal) {
                return compareScaled(&middle_power, middle_shift, &input.significand, input.exponent);
            }
            // SAFETY: `Big.multiply` initializes the active limbs before comparison.
            var product: Big = undefined;
            Big.multiply(&middle_power, &input.significand, &product);
            var one: Big = undefined;
            Big.one(&one);
            return compareScaled(&product, middle_shift + input.exponent, &one, 0);
        }

        fn overflowThresholdRelation(input: *const Decoded, degree: u64, info: FormatInfo) i8 {
            // RN-even overflow begins at max_finite + half its final ulp.  Written
            // exactly, that midpoint is `(2^(p+1)-1) * 2^(emax-p)`.
            var threshold_numerator: Big = undefined;
            Big.zero(&threshold_numerator);
            threshold_numerator.setBit(@intCast(info.precision_bits + 1));
            threshold_numerator.decrement();
            // SAFETY: `Big.power` initializes the active limbs before multiplication.
            var threshold_power: Big = undefined;
            Big.power(&threshold_numerator, degree, &threshold_power);
            const threshold_shift =
                (info.maximum_exponent - @as(i64, info.precision_bits)) * @as(i64, @intCast(degree));

            // SAFETY: `Big.multiply` initializes the active limbs before comparison.
            var product: Big = undefined;
            Big.multiply(&threshold_power, &input.significand, &product);
            var one: Big = undefined;
            Big.one(&one);
            return compareScaled(&product, threshold_shift + input.exponent, &one, 0);
        }

        fn powerOfTwoBits(exponent: i64, info: FormatInfo) u64 {
            const subnormal_exponent = info.minimum_exponent - @as(i64, info.fraction_bits);
            if (exponent < subnormal_exponent) return 0;
            if (exponent > info.maximum_exponent) return info.infinity_bits;
            if (exponent < info.minimum_exponent) {
                const bit_index: u6 = @intCast(exponent - subnormal_exponent);
                return @as(u64, 1) << bit_index;
            }
            const exponent_field: u64 = @intCast(exponent + info.exponent_bias);
            const fraction_shift: u6 = @intCast(info.fraction_bits);
            return exponent_field << fraction_shift;
        }

        /// Proposes a search start from a rounded estimate and certifies it with
        /// exact candidate comparisons.  The estimate never decides a result: a
        /// proposal counts only when `C(lower) <= 0 < C(lower + 1)` holds exactly,
        /// which is the same invariant the binary search maintains.  Every other
        /// case, including a proposal outside the search window, returns false and
        /// leaves the exact search to decide.
        fn seedLowerBound(
            comptime Format: type,
            value_bits: Format.Bits,
            degree: u64,
            reciprocal: bool,
            input: *const Decoded,
            lower_bound: u64,
            upper_bound: u64,
            lower: *u64,
            lower_relation: *i8,
        ) bool {
            @setFloatMode(.strict);
            const info = formatInfo(Format);
            const input_magnitude: f64 = @bitCast(@as(u64, value_bits) & ~info.sign_mask);
            const sign: f64 = if (reciprocal) -1.0 else 1.0;
            const step: f64 = sign / @as(f64, @floatFromInt(degree));
            const estimate: f64 = std.math.pow(f64, input_magnitude, step);
            const estimate_bits = @as(u64, @bitCast(estimate));
            if (estimate_bits < lower_bound or estimate_bits > upper_bound) return false;

            if (estimate_bits > lower_bound) {
                const below = candidateRelation(estimate_bits - 1, input, degree, reciprocal, info);
                if (below > 0) return false;
                const above = candidateRelation(estimate_bits, input, degree, reciprocal, info);
                if (above <= 0) return false;
                lower.* = estimate_bits - 1;
                lower_relation.* = below;
                return true;
            }

            const at = candidateRelation(estimate_bits, input, degree, reciprocal, info);
            if (at > 0) return false;
            lower.* = estimate_bits;
            lower_relation.* = at;
            return true;
        }

        fn rootFiniteBits(comptime Format: type, value_bits: Format.Bits, exponent: i64) u64 {
            const info = formatInfo(Format);
            var input: Decoded = undefined;
            decodeFinite(Format, value_bits, &input);
            std.debug.assert(!input.significand.isZero());
            const degree = magnitude(exponent);
            const reciprocal = exponent < 0;

            // If 2^E <= x < 2^(E+1), positive roots lie in
            // [2^floor(E/d), 2^(floor(E/d)+1)). Reciprocal roots lie in
            // (2^floor((-E-1)/d), 2^(floor(-E/d)+1)).
            // Clamp only the encoded endpoints; exact candidate and midpoint
            // comparisons still decide the rounded answer and overflow.
            const input_floor = input.exponent +
                @as(i64, @intCast(input.significand.bitLength())) - 1;
            const divisor: i64 = @intCast(degree);
            const lower_exponent = @divFloor(
                if (reciprocal) -input_floor - 1 else input_floor,
                divisor,
            );
            const upper_exponent = @divFloor(
                if (reciprocal) -input_floor else input_floor,
                divisor,
            ) + 1;
            const upper_bits = powerOfTwoBits(upper_exponent, info);
            std.debug.assert(upper_bits != 0);
            var lower = @min(powerOfTwoBits(lower_exponent, info), info.maximum_finite_bits);
            const upper_bound = @min(upper_bits - 1, info.maximum_finite_bits);
            std.debug.assert(lower <= upper_bound);

            // A rounded estimate can only start the search; exact comparisons
            // still decide the rounded answer and overflow.  Degrees below two
            // and every rejected proposal use the exact binary search.
            var lower_relation: i8 = undefined;
            var seeded = false;
            if (Format == common.D and degree >= 2) {
                seeded = seedLowerBound(
                    Format,
                    value_bits,
                    degree,
                    reciprocal,
                    &input,
                    lower,
                    upper_bound,
                    &lower,
                    &lower_relation,
                );
            }
            if (!seeded) {
                var search_upper = upper_bound;
                while (lower < search_upper) {
                    const span = search_upper - lower;
                    const middle = lower + span / 2 + (span & 1);
                    const relation = candidateRelation(middle, &input, degree, reciprocal, info);
                    if (relation <= 0) {
                        lower = middle;
                    } else {
                        search_upper = middle - 1;
                    }
                }
                lower_relation = candidateRelation(lower, &input, degree, reciprocal, info);
            }

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
            var decoded: Decoded = undefined;
            decodeFinite(Format, Format.toBits(base), &decoded);
            std.debug.assert(!decoded.significand.isZero());
            // SAFETY: Both sign paths initialize the rational's active limbs and shift.
            var rational: Rational = undefined;
            rationalFromPower(&decoded, exponent, &rational);
            return roundRational(Format, &rational, false);
        }

        /// Exact, one-round compoundn for finite `x > -1` and `|exponent| <= 64`.
        fn compoundFinite(comptime Format: type, value: Format.Float, exponent: i64) Format.Float {
            std.debug.assert(exponent >= -max_degree);
            std.debug.assert(exponent <= max_degree);
            var base: Decoded = undefined;
            makeCompoundBase(Format, Format.toBits(value), &base);
            // SAFETY: Both sign paths initialize the rational's active limbs and shift.
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
    var denominator: Big = undefined;
    Big.one(&denominator);
    denominator.setBit(1024);
    var numerator: Big = undefined;
    Big.fromU64(49, &numerator);
    numerator.setBit(1029);
    numerator.setBit(1026);
    numerator.setBit(1024);
    var quotient: Big = undefined;
    var remainder: Big = undefined;
    Big.divideRemainder(&numerator, &denominator, &quotient, &remainder);
    var expected_quotient: Big = undefined;
    var expected_remainder: Big = undefined;
    Big.fromU64(37, &expected_quotient);
    Big.fromU64(12, &expected_remainder);
    try std.testing.expectEqual(@as(i8, 0), Big.compare(&quotient, &expected_quotient));
    try std.testing.expectEqual(@as(i8, 0), Big.compare(&remainder, &expected_remainder));
    Big.divideRemainder(&expected_remainder, &denominator, &quotient, &remainder);
    try std.testing.expect(quotient.isZero());
    try std.testing.expectEqual(@as(i8, 0), Big.compare(&remainder, &expected_remainder));
}

test "shifted exact integer initializes gaps and carry limbs" {
    const Big = Arithmetic(32).Big;
    var source: Big = undefined;
    Big.fromU64(0x8000_0000_0000_0001, &source);
    source.setBit(64);

    for ([_]usize{ 0, 1, 63, 64, 65, 1024 }) |shift| {
        var expected: Big = undefined;
        Big.zero(&expected);
        for ([_]usize{ 0, 63, 64 }) |bit_index| {
            expected.setBit(bit_index + shift);
        }

        var actual: Big = undefined;
        Big.one(&actual);
        actual.setBit(1900);
        Big.shifted(&source, shift, &actual);
        try std.testing.expectEqual(@as(i8, 0), Big.compare(&expected, &actual));
    }

    var zero: Big = undefined;
    Big.zero(&zero);
    var actual: Big = undefined;
    Big.one(&actual);
    Big.shifted(&zero, 65, &actual);
    try std.testing.expect(actual.isZero());
}
