//! Exact-bit tests for the optional IEC 60559 magnitude and payload families.
//!
//! Payload values are nonnegative integral floating values. The quiet family
//! accepts zero through the format's quiet-bit-minus-one mask. The signaling
//! family accepts one through that mask, because signaling payload zero would
//! have the infinity encoding. Both signed zeros represent payload zero.

const std = @import("std");
const optional = @import("ieee_c23_optional.zig");

fn expectHBits(expected: u16, actual: f16) !void {
    try std.testing.expectEqual(expected, @as(u16, @bitCast(actual)));
}

fn expectFBits(expected: u32, actual: f32) !void {
    try std.testing.expectEqual(expected, @as(u32, @bitCast(actual)));
}

fn expectDBits(expected: u64, actual: f64) !void {
    try std.testing.expectEqual(expected, @as(u64, @bitCast(actual)));
}

fn expectHPayload(expected: u16, actual: f16) !void {
    try std.testing.expectEqual(@as(f16, @floatFromInt(expected)), actual);
}

fn expectFPayload(expected: u32, actual: f32) !void {
    try std.testing.expectEqual(@as(f64, @floatFromInt(expected)), @as(f64, actual));
}

fn expectDPayload(expected: u64, actual: f64) !void {
    try std.testing.expectEqual(@as(f64, @floatFromInt(expected)), actual);
}

test "totalOrderMagh uses magnitude bits for values zeros infinities and NaNs" {
    const positive_one: f16 = @bitCast(@as(u16, 0x3c00));
    const negative_one: f16 = @bitCast(@as(u16, 0xbc00));
    const positive_two: f16 = @bitCast(@as(u16, 0x4000));
    const positive_zero: f16 = @bitCast(@as(u16, 0x0000));
    const negative_zero: f16 = @bitCast(@as(u16, 0x8000));
    const positive_infinity: f16 = @bitCast(@as(u16, 0x7c00));
    const negative_infinity: f16 = @bitCast(@as(u16, 0xfc00));
    const signaling_low: f16 = @bitCast(@as(u16, 0x7c01));
    const signaling_high: f16 = @bitCast(@as(u16, 0x7dff));
    const quiet_low: f16 = @bitCast(@as(u16, 0x7e01));
    const quiet_high: f16 = @bitCast(@as(u16, 0x7fff));

    try std.testing.expect(optional.totalOrderMagh(negative_one, positive_one));
    try std.testing.expect(optional.totalOrderMagh(positive_one, negative_one));
    try std.testing.expect(!optional.totalOrderMagh(positive_two, positive_one));
    try std.testing.expect(optional.totalOrderMagh(negative_zero, positive_zero));
    try std.testing.expect(optional.totalOrderMagh(positive_zero, negative_zero));
    try std.testing.expect(optional.totalOrderMagh(negative_infinity, positive_infinity));
    try std.testing.expect(optional.totalOrderMagh(signaling_low, signaling_high));
    try std.testing.expect(optional.totalOrderMagh(signaling_high, quiet_low));
    try std.testing.expect(optional.totalOrderMagh(quiet_low, quiet_high));
    try std.testing.expect(!optional.totalOrderMagh(quiet_high, quiet_low));
}

test "totalOrderMagf uses magnitude bits for values zeros infinities and NaNs" {
    const positive_one: f32 = @bitCast(@as(u32, 0x3f80_0000));
    const negative_one: f32 = @bitCast(@as(u32, 0xbf80_0000));
    const positive_two: f32 = @bitCast(@as(u32, 0x4000_0000));
    const positive_zero: f32 = @bitCast(@as(u32, 0x0000_0000));
    const negative_zero: f32 = @bitCast(@as(u32, 0x8000_0000));
    const positive_infinity: f32 = @bitCast(@as(u32, 0x7f80_0000));
    const negative_infinity: f32 = @bitCast(@as(u32, 0xff80_0000));
    const signaling_low: f32 = @bitCast(@as(u32, 0x7f80_0001));
    const signaling_high: f32 = @bitCast(@as(u32, 0x7fbf_ffff));
    const quiet_low: f32 = @bitCast(@as(u32, 0x7fc0_0001));
    const quiet_high: f32 = @bitCast(@as(u32, 0x7fff_ffff));

    try std.testing.expect(optional.totalOrderMagf(negative_one, positive_one));
    try std.testing.expect(optional.totalOrderMagf(positive_one, negative_one));
    try std.testing.expect(!optional.totalOrderMagf(positive_two, positive_one));
    try std.testing.expect(optional.totalOrderMagf(negative_zero, positive_zero));
    try std.testing.expect(optional.totalOrderMagf(positive_zero, negative_zero));
    try std.testing.expect(optional.totalOrderMagf(negative_infinity, positive_infinity));
    try std.testing.expect(optional.totalOrderMagf(signaling_low, signaling_high));
    try std.testing.expect(optional.totalOrderMagf(signaling_high, quiet_low));
    try std.testing.expect(optional.totalOrderMagf(quiet_low, quiet_high));
    try std.testing.expect(!optional.totalOrderMagf(quiet_high, quiet_low));
}

test "totalOrderMagd uses magnitude bits for values zeros infinities and NaNs" {
    const positive_one: f64 = @bitCast(@as(u64, 0x3ff0_0000_0000_0000));
    const negative_one: f64 = @bitCast(@as(u64, 0xbff0_0000_0000_0000));
    const positive_two: f64 = @bitCast(@as(u64, 0x4000_0000_0000_0000));
    const positive_zero: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0000));
    const negative_zero: f64 = @bitCast(@as(u64, 0x8000_0000_0000_0000));
    const positive_infinity: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0000));
    const negative_infinity: f64 = @bitCast(@as(u64, 0xfff0_0000_0000_0000));
    const signaling_low: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0001));
    const signaling_high: f64 = @bitCast(@as(u64, 0x7ff7_ffff_ffff_ffff));
    const quiet_low: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_0001));
    const quiet_high: f64 = @bitCast(@as(u64, 0x7fff_ffff_ffff_ffff));

    try std.testing.expect(optional.totalOrderMagd(negative_one, positive_one));
    try std.testing.expect(optional.totalOrderMagd(positive_one, negative_one));
    try std.testing.expect(!optional.totalOrderMagd(positive_two, positive_one));
    try std.testing.expect(optional.totalOrderMagd(negative_zero, positive_zero));
    try std.testing.expect(optional.totalOrderMagd(positive_zero, negative_zero));
    try std.testing.expect(optional.totalOrderMagd(negative_infinity, positive_infinity));
    try std.testing.expect(optional.totalOrderMagd(signaling_low, signaling_high));
    try std.testing.expect(optional.totalOrderMagd(signaling_high, quiet_low));
    try std.testing.expect(optional.totalOrderMagd(quiet_low, quiet_high));
    try std.testing.expect(!optional.totalOrderMagd(quiet_high, quiet_low));
}

test "getPayloadh extracts q and s NaN diagnostic fields" {
    try expectHPayload(0x0055, optional.getPayloadh(@bitCast(@as(u16, 0x7e55))));
    try expectHPayload(0x0155, optional.getPayloadh(@bitCast(@as(u16, 0x7d55))));
    try expectHPayload(0x0155, optional.getPayloadh(@bitCast(@as(u16, 0xfd55))));
    try std.testing.expectEqual(@as(f16, -1.0), optional.getPayloadh(@bitCast(@as(u16, 0x0000))));
    try std.testing.expectEqual(@as(f16, -1.0), optional.getPayloadh(@bitCast(@as(u16, 0x7c00))));
}

test "getPayloadf extracts q and s NaN diagnostic fields" {
    try expectFPayload(0x0001_2345, optional.getPayloadf(@bitCast(@as(u32, 0x7fc1_2345))));
    try expectFPayload(0x003f_ffff, optional.getPayloadf(@bitCast(@as(u32, 0x7fbf_ffff))));
    try expectFPayload(0x003f_ffff, optional.getPayloadf(@bitCast(@as(u32, 0xffbf_ffff))));
    try std.testing.expectEqual(@as(f32, -1.0), optional.getPayloadf(@bitCast(@as(u32, 0x8000_0000))));
    try std.testing.expectEqual(@as(f32, -1.0), optional.getPayloadf(@bitCast(@as(u32, 0x7f80_0000))));
}

test "getPayloadd extracts q and s NaN diagnostic fields" {
    try expectDPayload(
        0x0001_2345_6789_abcd,
        optional.getPayloadd(@bitCast(@as(u64, 0x7ff9_2345_6789_abcd))),
    );
    try expectDPayload(
        0x0007_ffff_ffff_ffff,
        optional.getPayloadd(@bitCast(@as(u64, 0x7ff7_ffff_ffff_ffff))),
    );
    try expectDPayload(
        0x0007_ffff_ffff_ffff,
        optional.getPayloadd(@bitCast(@as(u64, 0xfff7_ffff_ffff_ffff))),
    );
    try std.testing.expectEqual(
        @as(f64, -1.0),
        optional.getPayloadd(@bitCast(@as(u64, 0x8000_0000_0000_0000))),
    );
    try std.testing.expectEqual(
        @as(f64, -1.0),
        optional.getPayloadd(@bitCast(@as(u64, 0x7ff0_0000_0000_0000))),
    );
}

test "setPayloadh constructs exact q and s NaNs at both payload limits" {
    const quiet_zero = optional.setPayloadh(@bitCast(@as(u16, 0x0000)));
    try std.testing.expect(quiet_zero.ok);
    try expectHBits(0x7e00, quiet_zero.value);

    const quiet_max = optional.setPayloadh(@floatFromInt(@as(u16, 511)));
    try std.testing.expect(quiet_max.ok);
    try expectHBits(0x7fff, quiet_max.value);

    const signaling_min = optional.setPayloadSigh(@floatFromInt(@as(u16, 1)));
    try std.testing.expect(signaling_min.ok);
    try expectHBits(0x7c01, signaling_min.value);

    const signaling_max = optional.setPayloadSigh(@floatFromInt(@as(u16, 511)));
    try std.testing.expect(signaling_max.ok);
    try expectHBits(0x7dff, signaling_max.value);

    const signed_zero = optional.setPayloadh(@bitCast(@as(u16, 0x8000)));
    try std.testing.expect(signed_zero.ok);
    try expectHBits(0x7e00, signed_zero.value);
}

test "setPayloadf constructs exact q and s NaNs at both payload limits" {
    const quiet_zero = optional.setPayloadf(@bitCast(@as(u32, 0x0000_0000)));
    try std.testing.expect(quiet_zero.ok);
    try expectFBits(0x7fc0_0000, quiet_zero.value);

    const quiet_max = optional.setPayloadf(@floatFromInt(@as(u32, 0x003f_ffff)));
    try std.testing.expect(quiet_max.ok);
    try expectFBits(0x7fff_ffff, quiet_max.value);

    const signaling_min = optional.setPayloadSigf(@floatFromInt(@as(u32, 1)));
    try std.testing.expect(signaling_min.ok);
    try expectFBits(0x7f80_0001, signaling_min.value);

    const signaling_max = optional.setPayloadSigf(@floatFromInt(@as(u32, 0x003f_ffff)));
    try std.testing.expect(signaling_max.ok);
    try expectFBits(0x7fbf_ffff, signaling_max.value);

    const signed_zero = optional.setPayloadf(@bitCast(@as(u32, 0x8000_0000)));
    try std.testing.expect(signed_zero.ok);
    try expectFBits(0x7fc0_0000, signed_zero.value);
}

test "setPayloadd constructs exact q and s NaNs at both payload limits" {
    const quiet_zero = optional.setPayloadd(@bitCast(@as(u64, 0x0000_0000_0000_0000)));
    try std.testing.expect(quiet_zero.ok);
    try expectDBits(0x7ff8_0000_0000_0000, quiet_zero.value);

    const quiet_max = optional.setPayloadd(
        @floatFromInt(@as(u64, 0x0007_ffff_ffff_ffff)),
    );
    try std.testing.expect(quiet_max.ok);
    try expectDBits(0x7fff_ffff_ffff_ffff, quiet_max.value);

    const signaling_min = optional.setPayloadSigd(@floatFromInt(@as(u64, 1)));
    try std.testing.expect(signaling_min.ok);
    try expectDBits(0x7ff0_0000_0000_0001, signaling_min.value);

    const signaling_max = optional.setPayloadSigd(
        @floatFromInt(@as(u64, 0x0007_ffff_ffff_ffff)),
    );
    try std.testing.expect(signaling_max.ok);
    try expectDBits(0x7ff7_ffff_ffff_ffff, signaling_max.value);

    const signed_zero = optional.setPayloadd(@bitCast(@as(u64, 0x8000_0000_0000_0000)));
    try std.testing.expect(signed_zero.ok);
    try expectDBits(0x7ff8_0000_0000_0000, signed_zero.value);
}

fn expectInvalidH(payload: f16) !void {
    const quiet = optional.setPayloadh(payload);
    try std.testing.expect(!quiet.ok);
    try expectHBits(0, quiet.value);

    const signaling = optional.setPayloadSigh(payload);
    try std.testing.expect(!signaling.ok);
    try expectHBits(0, signaling.value);
}

fn expectInvalidF(payload: f32) !void {
    const quiet = optional.setPayloadf(payload);
    try std.testing.expect(!quiet.ok);
    try expectFBits(0, quiet.value);

    const signaling = optional.setPayloadSigf(payload);
    try std.testing.expect(!signaling.ok);
    try expectFBits(0, signaling.value);
}

fn expectInvalidD(payload: f64) !void {
    const quiet = optional.setPayloadd(payload);
    try std.testing.expect(!quiet.ok);
    try expectDBits(0, quiet.value);

    const signaling = optional.setPayloadSigd(payload);
    try std.testing.expect(!signaling.ok);
    try expectDBits(0, signaling.value);
}

test "setPayloadh rejects fractional negative infinite NaN and out-of-range values" {
    try expectInvalidH(1.5);
    try expectInvalidH(-1.0);
    try expectInvalidH(@bitCast(@as(u16, 0x7c00)));
    try expectInvalidH(@bitCast(@as(u16, 0x7e01)));
    try expectInvalidH(@bitCast(@as(u16, 0x6000)));

    const zero = optional.setPayloadSigh(@bitCast(@as(u16, 0)));
    try std.testing.expect(!zero.ok);
    try expectHBits(0, zero.value);
}

test "setPayloadf rejects fractional negative infinite NaN and out-of-range values" {
    try expectInvalidF(1.5);
    try expectInvalidF(-1.0);
    try expectInvalidF(@bitCast(@as(u32, 0x7f80_0000)));
    try expectInvalidF(@bitCast(@as(u32, 0x7fc0_0001)));
    try expectInvalidF(@bitCast(@as(u32, 0x4b80_0000)));

    const zero = optional.setPayloadSigf(@bitCast(@as(u32, 0)));
    try std.testing.expect(!zero.ok);
    try expectFBits(0, zero.value);
}

test "setPayloadd rejects fractional negative infinite NaN and out-of-range values" {
    try expectInvalidD(1.5);
    try expectInvalidD(-1.0);
    try expectInvalidD(@bitCast(@as(u64, 0x7ff0_0000_0000_0000)));
    try expectInvalidD(@bitCast(@as(u64, 0x7ff8_0000_0000_0001)));
    try expectInvalidD(@bitCast(@as(u64, 0x4320_0000_0000_0000)));

    const zero = optional.setPayloadSigd(@bitCast(@as(u64, 0)));
    try std.testing.expect(!zero.ok);
    try expectDBits(0, zero.value);
}
