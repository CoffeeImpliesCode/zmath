//! Optional IEC 60559 magnitude-order and NaN-payload operations.
//!
//! The diagnostic payload is the fraction field below the quiet-NaN bit. The
//! quiet bit is not part of the payload. This gives 9 payload bits for `f16`,
//! 22 for `f32`, and 51 for `f64`.
//!
//! `setPayload*` accepts either floating encoding of zero and positive integral
//! payloads through 511, 4_194_303, and 2_251_799_813_685_247 respectively.
//! `setPayloadSig*` uses the same encoding but starts at one because a zero
//! fraction with the quiet bit clear encodes infinity. Invalid inputs produce
//! positive zero with `ok = false`; valid results are sign-positive NaNs.

const std = @import("std");
const common = @import("scalar_c23_common.zig");
const c23_types = @import("scalar_c23_types.zig");

const H = OptionalFormat(common.H, .{
    .sign_mask = 0x8000,
    .exponent_mask = 0x7c00,
    .fraction_mask = 0x03ff,
    .quiet_nan_bit = 0x0200,
});

const F = OptionalFormat(common.F, .{
    .sign_mask = 0x8000_0000,
    .exponent_mask = 0x7f80_0000,
    .fraction_mask = 0x007f_ffff,
    .quiet_nan_bit = 0x0040_0000,
});

const D = OptionalFormat(common.D, .{
    .sign_mask = 0x8000_0000_0000_0000,
    .exponent_mask = 0x7ff0_0000_0000_0000,
    .fraction_mask = 0x000f_ffff_ffff_ffff,
    .quiet_nan_bit = 0x0008_0000_0000_0000,
});

/// Returns true when the magnitude of `left` is at most the magnitude of
/// `right`, comparing the sign-cleared encodings directly.
///
/// The comparison is a total order over the magnitudes and deliberately does
/// not distinguish `-0` from `+0`; use `ieee.totalOrderh` when the sign must
/// take part.
pub fn totalOrderMagh(left: f16, right: f16) bool {
    return H.totalOrderMag(left, right);
}

/// Returns true when the magnitude of `left` is at most the magnitude of
/// `right`, at binary32 width, with the same total-order and signed-zero rule
/// as `totalOrderMagh`.
pub fn totalOrderMagf(left: f32, right: f32) bool {
    return F.totalOrderMag(left, right);
}

/// Returns true when the magnitude of `left` is at most the magnitude of
/// `right`, at binary64 width, with the same total-order and signed-zero rule
/// as `totalOrderMagh`.
pub fn totalOrderMagd(left: f64, right: f64) bool {
    return D.totalOrderMag(left, right);
}

/// Returns the diagnostic payload of a NaN `value` as a binary16 value, and
/// `-1.0` for every non-NaN input.
///
/// The 9-bit payload is exactly representable, so the result is exact.
/// `-1.0` is a sentinel for "the operand was not a NaN"; no payload is
/// negative, so a caller can always tell the two apart.
pub fn getPayloadh(value: f16) f16 {
    return H.getPayload(value);
}

/// Returns the diagnostic payload of a NaN `value` as a binary32 value, and
/// `-1.0` for every non-NaN input, as `getPayloadh` does at binary32 width.
pub fn getPayloadf(value: f32) f32 {
    return F.getPayload(value);
}

/// Returns the diagnostic payload of a NaN `value` as a binary64 value, and
/// `-1.0` for every non-NaN input, as `getPayloadh` does at binary64 width.
pub fn getPayloadd(value: f64) f64 {
    return D.getPayload(value);
}

/// Builds a positive quiet NaN carrying the diagnostic payload `payload` and
/// returns it with `ok` true.
///
/// A `payload` that is neither a signed zero nor an in-range positive integral
/// value returns `+0.0` with `ok` false; that `0.0` is a sentinel and not a
/// payload of zero, so a caller must read `ok` before using the result.
pub fn setPayloadh(payload: f16) c23_types.PayloadResult(f16) {
    return H.setPayload(payload, false);
}

/// Builds a positive quiet NaN from a binary32 `payload`, returning `+0.0`
/// with `ok` false on the same rejected inputs as `setPayloadh`.
pub fn setPayloadf(payload: f32) c23_types.PayloadResult(f32) {
    return F.setPayload(payload, false);
}

/// Builds a positive quiet NaN from a binary64 `payload`, returning `+0.0`
/// with `ok` false on the same rejected inputs as `setPayloadh`.
pub fn setPayloadd(payload: f64) c23_types.PayloadResult(f64) {
    return D.setPayload(payload, false);
}

/// Builds a positive SIGNALING NaN carrying the diagnostic payload `payload`,
/// and returns it with `ok` true.
///
/// The encoding is the same as `setPayloadh` with the quiet bit left clear, so
/// a zero payload is not representable: `payload` of `+0.0` or `-0.0` returns
/// `+0.0` with `ok` false, exactly as any other rejected input does.
pub fn setPayloadSigh(payload: f16) c23_types.PayloadResult(f16) {
    return H.setPayload(payload, true);
}

/// Builds a positive signaling binary32 NaN from `payload`, returning `+0.0`
/// with `ok` false on the same rejected inputs as `setPayloadSigh`.
pub fn setPayloadSigf(payload: f32) c23_types.PayloadResult(f32) {
    return F.setPayload(payload, true);
}

/// Builds a positive signaling binary64 NaN from `payload`, returning `+0.0`
/// with `ok` false on the same rejected inputs as `setPayloadSigh`.
pub fn setPayloadSigd(payload: f64) c23_types.PayloadResult(f64) {
    return D.setPayload(payload, true);
}

fn OptionalFormat(
    comptime Format: type,
    comptime options: struct {
        sign_mask: Format.Bits,
        exponent_mask: Format.Bits,
        fraction_mask: Format.Bits,
        quiet_nan_bit: Format.Bits,
    },
) type {
    return struct {
        const Float = Format.Float;
        const Bits = Format.Bits;
        const Shift = std.math.Log2Int(Bits);

        const fraction_bits: comptime_int = @ctz(options.quiet_nan_bit) + 1;
        const payload_bits: comptime_int = fraction_bits - 1;
        const payload_mask: Bits = options.quiet_nan_bit - 1;
        const exponent_max: u16 = @intCast(options.exponent_mask >> fraction_bits);
        const exponent_bias: u16 = exponent_max >> 1;

        comptime {
            std.debug.assert(options.fraction_mask == (options.quiet_nan_bit << 1) - 1);
            std.debug.assert(options.exponent_mask & options.fraction_mask == 0);
            std.debug.assert(options.sign_mask & options.exponent_mask == 0);
        }

        fn totalOrderMag(left: Float, right: Float) bool {
            const left_magnitude = Format.toBits(left) & ~options.sign_mask;
            const right_magnitude = Format.toBits(right) & ~options.sign_mask;

            // Positive IEEE binary encodings have the required magnitude order:
            // numbers, infinity, signaling NaNs, then quiet NaNs. Within each
            // NaN class, the diagnostic payload increases with the encoding.
            return left_magnitude <= right_magnitude;
        }

        fn getPayload(value: Float) Float {
            const bits = Format.toBits(value);
            if (!Format.isNaNBits(bits)) return -1.0;

            return @floatFromInt(bits & payload_mask);
        }

        fn setPayload(
            payload: Float,
            comptime signaling: bool,
        ) c23_types.PayloadResult(Float) {
            const payload_integer = payloadInteger(payload) orelse return invalidPayload();
            if (signaling and payload_integer == 0) return invalidPayload();

            const nan_bits = options.exponent_mask |
                payload_integer |
                if (signaling) 0 else options.quiet_nan_bit;
            return .{
                .value = Format.fromBits(nan_bits),
                .ok = true,
            };
        }

        fn payloadInteger(value: Float) ?Bits {
            const bits = Format.toBits(value);
            const magnitude = bits & ~options.sign_mask;

            // Both signed zeros represent the admissible integer payload zero.
            if (magnitude == 0) return 0;
            if (bits & options.sign_mask != 0) return null;

            const exponent_field: u16 = @intCast(
                (bits & options.exponent_mask) >> fraction_bits,
            );
            if (exponent_field == 0) return null;
            if (exponent_field == exponent_max) return null;
            if (exponent_field < exponent_bias) return null;

            const exponent = exponent_field - exponent_bias;
            if (exponent >= payload_bits) return null;

            const shift: Shift = @intCast(fraction_bits - exponent);
            const significand = (@as(Bits, 1) << fraction_bits) |
                (bits & options.fraction_mask);
            const discarded_mask = (@as(Bits, 1) << shift) - 1;
            if (significand & discarded_mask != 0) return null;

            const integer = significand >> shift;
            if (integer > payload_mask) return null;
            return integer;
        }

        fn invalidPayload() c23_types.PayloadResult(Float) {
            return .{
                .value = Format.fromBits(0),
                .ok = false,
            };
        }
    };
}
