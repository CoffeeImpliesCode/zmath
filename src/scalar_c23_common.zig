const std = @import("std");

pub const Class = enum {
    zero,
    subnormal,
    normal,
    infinity,
    quiet_nan,
    signaling_nan,
};

pub const H = IeeeFormat(f16, u16, .{
    .sign_mask = 0x8000,
    .exponent_mask = 0x7c00,
    .fraction_mask = 0x03ff,
    .quiet_nan_bit = 0x0200,
});

pub const F = IeeeFormat(f32, u32, .{
    .sign_mask = 0x8000_0000,
    .exponent_mask = 0x7f80_0000,
    .fraction_mask = 0x007f_ffff,
    .quiet_nan_bit = 0x0040_0000,
});

pub const D = IeeeFormat(f64, u64, .{
    .sign_mask = 0x8000_0000_0000_0000,
    .exponent_mask = 0x7ff0_0000_0000_0000,
    .fraction_mask = 0x000f_ffff_ffff_ffff,
    .quiet_nan_bit = 0x0008_0000_0000_0000,
});

fn IeeeFormat(
    comptime FloatType: type,
    comptime BitsType: type,
    comptime options: struct {
        sign_mask: BitsType,
        exponent_mask: BitsType,
        fraction_mask: BitsType,
        quiet_nan_bit: BitsType,
    },
) type {
    return struct {
        pub const Float = FloatType;
        pub const Bits = BitsType;

        const sign_mask: Bits = options.sign_mask;
        const exponent_mask: Bits = options.exponent_mask;
        const fraction_mask: Bits = options.fraction_mask;
        const quiet_nan_bit: Bits = options.quiet_nan_bit;
        const magnitude_mask: Bits = ~sign_mask;
        const infinity_bits: Bits = exponent_mask;
        const canonical_nan_bits: Bits = infinity_bits | quiet_nan_bit;

        pub fn toBits(value: Float) Bits {
            return @bitCast(value);
        }

        pub fn fromBits(bits: Bits) Float {
            return @bitCast(bits);
        }

        pub fn classifyBits(bits: Bits) Class {
            const exponent = bits & exponent_mask;
            const fraction = bits & fraction_mask;

            if (exponent == 0) {
                if (fraction == 0) return .zero;
                return .subnormal;
            }
            if (exponent != exponent_mask) return .normal;
            if (fraction == 0) return .infinity;
            if (fraction & quiet_nan_bit != 0) return .quiet_nan;
            return .signaling_nan;
        }

        pub fn isNaNBits(bits: Bits) bool {
            return bits & magnitude_mask > infinity_bits;
        }

        pub fn quietNaNBits(bits: Bits) Bits {
            std.debug.assert(isNaNBits(bits));
            return bits | quiet_nan_bit;
        }

        pub fn canonicalNaN() Float {
            return fromBits(canonical_nan_bits);
        }

        pub fn isZeroBits(bits: Bits) bool {
            return bits & magnitude_mask == 0;
        }

        pub fn isPositiveZeroBits(bits: Bits) bool {
            return bits == 0;
        }

        pub fn isNegativeZeroBits(bits: Bits) bool {
            return bits == sign_mask;
        }

        pub fn isInfinityBits(bits: Bits) bool {
            return bits & magnitude_mask == infinity_bits;
        }

        pub fn isPositiveInfinityBits(bits: Bits) bool {
            return bits == infinity_bits;
        }

        pub fn isNegativeInfinityBits(bits: Bits) bool {
            return bits == sign_mask | infinity_bits;
        }

        /// Returns a monotonic key for a non-NaN IEEE value.
        pub fn finiteOrderKey(bits: Bits) Bits {
            std.debug.assert(!isNaNBits(bits));
            if (bits & sign_mask != 0) return ~bits;
            return bits | sign_mask;
        }

        pub fn firstNaN2(first: Float, second: Float) ?Float {
            const first_bits = toBits(first);
            if (isNaNBits(first_bits)) return fromBits(quietNaNBits(first_bits));

            const second_bits = toBits(second);
            if (isNaNBits(second_bits)) return fromBits(quietNaNBits(second_bits));
            return null;
        }

        pub fn firstNaN3(first: Float, second: Float, third: Float) ?Float {
            if (firstNaN2(first, second)) |value| return value;

            const third_bits = toBits(third);
            if (isNaNBits(third_bits)) return fromBits(quietNaNBits(third_bits));
            return null;
        }
    };
}
