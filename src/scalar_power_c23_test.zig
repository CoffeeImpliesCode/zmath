const std = @import("std");
const common = @import("scalar_c23_common.zig");
const c23 = @import("scalar_c23_types.zig");
const power = @import("scalar_power_c23.zig");

fn expectPownH(base_bits: u16, exponent: i64, expected_bits: u16) !void {
    const result = power.pownh(common.H.fromBits(base_bits), exponent);
    try std.testing.expect(result.ok);
    try std.testing.expectEqual(expected_bits, common.H.toBits(result.value));
}

fn expectPownF(base_bits: u32, exponent: i64, expected_bits: u32) !void {
    const result = power.pownf(common.F.fromBits(base_bits), exponent);
    try std.testing.expect(result.ok);
    try std.testing.expectEqual(expected_bits, common.F.toBits(result.value));
}

fn expectPownD(base_bits: u64, exponent: i64, expected_bits: u64) !void {
    const result = power.pownd(common.D.fromBits(base_bits), exponent);
    try std.testing.expect(result.ok);
    try std.testing.expectEqual(expected_bits, common.D.toBits(result.value));
}

fn expectCompoundH(value_bits: u16, exponent: i64, expected_bits: u16) !void {
    const result = power.compoundnh(common.H.fromBits(value_bits), exponent);
    try std.testing.expect(result.ok);
    try std.testing.expectEqual(expected_bits, common.H.toBits(result.value));
}

fn expectCompoundF(value_bits: u32, exponent: i64, expected_bits: u32) !void {
    const result = power.compoundnf(common.F.fromBits(value_bits), exponent);
    try std.testing.expect(result.ok);
    try std.testing.expectEqual(expected_bits, common.F.toBits(result.value));
}

fn expectCompoundD(value_bits: u64, exponent: i64, expected_bits: u64) !void {
    const result = power.compoundnd(common.D.fromBits(value_bits), exponent);
    try std.testing.expect(result.ok);
    try std.testing.expectEqual(expected_bits, common.D.toBits(result.value));
}

fn expectRootH(value_bits: u16, exponent: i64, expected_bits: u16) !void {
    const result = power.rootnh(common.H.fromBits(value_bits), exponent);
    try std.testing.expect(result.ok);
    try std.testing.expectEqual(expected_bits, common.H.toBits(result.value));
}

fn expectRootF(value_bits: u32, exponent: i64, expected_bits: u32) !void {
    const result = power.rootnf(common.F.fromBits(value_bits), exponent);
    try std.testing.expect(result.ok);
    try std.testing.expectEqual(expected_bits, common.F.toBits(result.value));
}

fn expectRootD(value_bits: u64, exponent: i64, expected_bits: u64) !void {
    const result = power.rootnd(common.D.fromBits(value_bits), exponent);
    try std.testing.expect(result.ok);
    try std.testing.expectEqual(expected_bits, common.D.toBits(result.value));
}

test "C23 power families expose PowerResult and scalar powr signatures" {
    comptime {
        std.debug.assert(@TypeOf(power.pownh) == fn (f16, i64) c23.PowerResult(f16));
        std.debug.assert(@TypeOf(power.pownf) == fn (f32, i64) c23.PowerResult(f32));
        std.debug.assert(@TypeOf(power.pownd) == fn (f64, i64) c23.PowerResult(f64));
        std.debug.assert(@TypeOf(power.powrh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(power.powrf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(power.powrd) == fn (f64, f64) f64);
        std.debug.assert(@TypeOf(power.compoundnh) == fn (f16, i64) c23.PowerResult(f16));
        std.debug.assert(@TypeOf(power.compoundnf) == fn (f32, i64) c23.PowerResult(f32));
        std.debug.assert(@TypeOf(power.compoundnd) == fn (f64, i64) c23.PowerResult(f64));
        std.debug.assert(@TypeOf(power.rootnh) == fn (f16, i64) c23.PowerResult(f16));
        std.debug.assert(@TypeOf(power.rootnf) == fn (f32, i64) c23.PowerResult(f32));
        std.debug.assert(@TypeOf(power.rootnd) == fn (f64, i64) c23.PowerResult(f64));
    }
}

test "bounded pown rounds exact powers and preserves special precedence" {
    try expectPownH(0x3c00, 0, 0x3c00);
    try expectPownH(0x4000, 3, 0x4800);
    try expectPownH(0xc000, 3, 0xc800);
    try expectPownH(0x4000, -3, 0x3000);
    try expectPownH(0x8000, 3, 0x8000);
    try expectPownH(0x8000, -2, 0x7c00);
    try expectPownH(0x7c00, -3, 0x0000);
    try expectPownH(0xfc00, 3, 0xfc00);

    try expectPownF(0x4000_0000, 3, 0x4100_0000);
    try expectPownF(0xc000_0000, 3, 0xc100_0000);
    try expectPownF(0x4000_0000, -3, 0x3e00_0000);
    try expectPownF(0x8000_0000, -2, 0x7f80_0000);
    try expectPownF(0x7f80_0000, -3, 0x0000_0000);

    try expectPownD(0x4000_0000_0000_0000, 3, 0x4020_0000_0000_0000);
    try expectPownD(0xc000_0000_0000_0000, 3, 0xc020_0000_0000_0000);
    try expectPownD(0x4000_0000_0000_0000, -3, 0x3fc0_0000_0000_0000);
    try expectPownD(0x8000_0000_0000_0000, -2, 0x7ff0_0000_0000_0000);
    try expectPownD(0x7ff0_0000_0000_0000, -3, 0x0000_0000_0000_0000);
}

test "pown and rootn apply ordered NaN and identity rules" {
    const pown_h_snan = power.pownh(common.H.fromBits(0x7d55), 0);
    try std.testing.expect(pown_h_snan.ok);
    try std.testing.expectEqual(@as(u16, 0x7f55), common.H.toBits(pown_h_snan.value));
    const pown_h_qnan = power.pownh(common.H.fromBits(0x7e55), 0);
    try std.testing.expectEqual(@as(u16, 0x3c00), common.H.toBits(pown_h_qnan.value));

    const pown_f_snan = power.pownf(common.F.fromBits(0x7f80_0123), 4);
    try std.testing.expectEqual(@as(u32, 0x7fc0_0123), common.F.toBits(pown_f_snan.value));
    const pown_f_snan_one = power.pownf(common.F.fromBits(0x7f80_0123), 1);
    try std.testing.expect(pown_f_snan_one.ok);
    try std.testing.expectEqual(@as(u32, 0x7fc0_0123), common.F.toBits(pown_f_snan_one.value));
    const pown_d_snan = power.pownd(common.D.fromBits(0x7ff0_0000_0000_0123), 4);
    try std.testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0123), common.D.toBits(pown_d_snan.value));

    const root_h_zero = power.rootnh(common.H.fromBits(0x7e55), 0);
    try std.testing.expect(root_h_zero.ok);
    try std.testing.expectEqual(@as(u16, 0x7e00), common.H.toBits(root_h_zero.value));
    const root_f_zero = power.rootnf(common.F.fromBits(0x7fc0_0123), 0);
    try std.testing.expectEqual(@as(u32, 0x7fc0_0000), common.F.toBits(root_f_zero.value));
    const root_d_zero = power.rootnd(common.D.fromBits(0x7ff0_0000_0000_0123), 0);
    try std.testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), common.D.toBits(root_d_zero.value));
    const root_h_snan_one = power.rootnh(common.H.fromBits(0x7d55), 1);
    try std.testing.expect(root_h_snan_one.ok);
    try std.testing.expectEqual(@as(u16, 0x7f55), common.H.toBits(root_h_snan_one.value));
}

test "bounded pown covers subnormal bases reciprocal powers halfway cases and 64" {
    try expectPownH(0x0001, 1, 0x0001);
    try expectPownH(0x0001, 2, 0x0000);
    try expectPownH(0x4000, 64, 0x7c00);
    try expectPownH(0x3800, 64, 0x0000);
    try expectPownH(0x3c01, 2, 0x3c02);

    try expectPownF(0x0000_0001, 1, 0x0000_0001);
    try expectPownF(0x0000_0001, 2, 0x0000_0000);
    try expectPownF(0x4000_0000, 64, 0x5f80_0000);
    try expectPownF(0x3f00_0000, 64, 0x1f80_0000);
    try expectPownF(0x3f80_0001, 2, 0x3f80_0002);

    try expectPownD(0x0000_0000_0000_0001, 1, 0x0000_0000_0000_0001);
    try expectPownD(0x0000_0000_0000_0001, 2, 0x0000_0000_0000_0000);
    try expectPownD(0x4000_0000_0000_0000, 64, 0x43f0_0000_0000_0000);
    try expectPownD(0x3fe0_0000_0000_0000, 64, 0x3bf0_0000_0000_0000);
    try expectPownD(0x3ff0_0000_0000_0001, 2, 0x3ff0_0000_0000_0002);
}

test "finite power identities preserve signs and round squares only once" {
    try expectPownH(0xbc01, 1, 0xbc01);
    try expectPownF(0xbf80_0001, 1, 0xbf80_0001);
    try expectPownD(0xbff0_0000_0000_0001, 1, 0xbff0_0000_0000_0001);
    try expectPownH(0xbc00, -63, 0xbc00);
    try expectPownF(0x3f80_0000, -64, 0x3f80_0000);
    try expectPownD(0xbff0_0000_0000_0000, 64, 0x3ff0_0000_0000_0000);
    try expectPownH(0x9400, 2, 0x0010);
    try expectPownF(0x1f80_0000, 2, 0x0020_0000);
    try expectPownD(0x1ff0_0000_0000_0000, 2, 0x0004_0000_0000_0000);
    try expectPownH(0xfbff, 2, 0x7c00);
    try expectPownF(0xff7f_ffff, 2, 0x7f80_0000);
    try expectPownD(0xffef_ffff_ffff_ffff, 2, 0x7ff0_0000_0000_0000);
}

test "bounded pown rejects both range boundaries with canonical quiet NaN" {
    const h_low = power.pownh(common.H.fromBits(0x3c00), -65);
    const h_high = power.pownh(common.H.fromBits(0x3c00), 65);
    try std.testing.expect(!h_low.ok);
    try std.testing.expect(!h_high.ok);
    try std.testing.expectEqual(@as(u16, 0x7e00), common.H.toBits(h_low.value));
    try std.testing.expectEqual(@as(u16, 0x7e00), common.H.toBits(h_high.value));

    const f_low = power.pownf(1.0, -65);
    const f_high = power.pownf(1.0, 65);
    try std.testing.expect(!f_low.ok);
    try std.testing.expect(!f_high.ok);
    try std.testing.expectEqual(@as(u32, 0x7fc0_0000), common.F.toBits(f_low.value));
    try std.testing.expectEqual(@as(u32, 0x7fc0_0000), common.F.toBits(f_high.value));

    const d_low = power.pownd(1.0, -65);
    const d_high = power.pownd(1.0, 65);
    try std.testing.expect(!d_low.ok);
    try std.testing.expect(!d_high.ok);
    try std.testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), common.D.toBits(d_low.value));
    try std.testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), common.D.toBits(d_high.value));
}

test "compoundn expands one plus x exactly before rounding" {
    try expectCompoundH(0x1000, 2, 0x3c01);
    try expectCompoundH(0x3c00, 3, 0x4800);
    try expectCompoundH(0x3c00, -3, 0x3000);
    try expectCompoundH(0xbc00, 3, 0x0000);
    try expectCompoundH(0xbc00, -3, 0x7c00);
    try expectCompoundH(0xbbff, 2, 0x0004);
    try expectCompoundH(0x7bff, -1, 0x0100);
    try expectCompoundH(0x7c00, 3, 0x7c00);
    try expectCompoundH(0x7c00, -3, 0x0000);

    try expectCompoundF(0x3380_0000, 2, 0x3f80_0001);
    try expectCompoundF(0x3f80_0000, 64, 0x5f80_0000);
    try expectCompoundF(0xbf80_0000, 3, 0x0000_0000);

    try expectCompoundD(0x3ca0_0000_0000_0000, 2, 0x3ff0_0000_0000_0001);
    try expectCompoundD(0x3ff0_0000_0000_0000, -64, 0x3bf0_0000_0000_0000);
}

test "tiny compound bases round to one but the degree-limit boundary does not" {
    try expectCompoundH(0x0020, 64, 0x3c00);
    try expectCompoundH(0x8020, 64, 0x3c00);
    try expectCompoundH(0x0020, -64, 0x3c00);
    try expectCompoundH(0x8020, -64, 0x3c00);
    try expectCompoundH(0x8080, 64, 0x3bff);
    try expectCompoundH(0x8080, -64, 0x3c01);
    try expectCompoundF(0x2f80_0000, 64, 0x3f80_0000);
    try expectCompoundF(0xaf80_0000, 64, 0x3f80_0000);
    try expectCompoundF(0x2f80_0000, -64, 0x3f80_0000);
    try expectCompoundF(0xaf80_0000, -64, 0x3f80_0000);
    try expectCompoundF(0xb080_0000, 64, 0x3f7f_ffff);
    try expectCompoundF(0xb080_0000, -64, 0x3f80_0001);
    try expectCompoundD(0x3c20_0000_0000_0000, 64, 0x3ff0_0000_0000_0000);
    try expectCompoundD(0xbc20_0000_0000_0000, 64, 0x3ff0_0000_0000_0000);
    try expectCompoundD(0x3c20_0000_0000_0000, -64, 0x3ff0_0000_0000_0000);
    try expectCompoundD(0xbc20_0000_0000_0000, -64, 0x3ff0_0000_0000_0000);
    try expectCompoundD(0xbc40_0000_0000_0000, 64, 0x3fef_ffff_ffff_ffff);
    try expectCompoundD(0xbc40_0000_0000_0000, -64, 0x3ff0_0000_0000_0001);
}

test "compoundn preserves domain and NaN precedence and rejects 65" {
    const h_domain = power.compoundnh(common.H.fromBits(0xc000), 1);
    try std.testing.expect(h_domain.ok);
    try std.testing.expectEqual(@as(u16, 0x7e00), common.H.toBits(h_domain.value));
    const h_snan = power.compoundnh(common.H.fromBits(0x7d55), 0);
    try std.testing.expect(h_snan.ok);
    try std.testing.expectEqual(@as(u16, 0x7f55), common.H.toBits(h_snan.value));
    const h_invalid = power.compoundnh(common.H.fromBits(0x0000), 65);
    try std.testing.expect(!h_invalid.ok);
    try std.testing.expectEqual(@as(u16, 0x7e00), common.H.toBits(h_invalid.value));
    const h_tiny_invalid = power.compoundnh(common.H.fromBits(0x0020), 65);
    try std.testing.expect(!h_tiny_invalid.ok);
    try std.testing.expectEqual(@as(u16, 0x7e00), common.H.toBits(h_tiny_invalid.value));

    const f_invalid = power.compoundnf(0.0, -65);
    try std.testing.expect(!f_invalid.ok);
    try std.testing.expectEqual(@as(u32, 0x7fc0_0000), common.F.toBits(f_invalid.value));
    const d_invalid = power.compoundnd(0.0, 65);
    try std.testing.expect(!d_invalid.ok);
    try std.testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), common.D.toBits(d_invalid.value));
}

test "rootn compares exact candidate powers for odd negatives reciprocal roots and ties" {
    try expectRootH(0x4800, 3, 0x4000);
    try expectRootH(0xc800, 3, 0xc000);
    try expectRootH(0x4400, 2, 0x4000);
    try expectRootH(0x4c00, -2, 0x3400);
    try expectRootH(0x8000, 3, 0x8000);
    try expectRootH(0x0000, -2, 0x7c00);
    try expectRootH(0x3c01, 2, 0x3c00);

    try expectRootF(0x4100_0000, 3, 0x4000_0000);
    try expectRootF(0xc100_0000, 3, 0xc000_0000);
    try expectRootF(0x4080_0000, 2, 0x4000_0000);
    try expectRootF(0x4180_0000, -2, 0x3e80_0000);
    try expectRootF(0x8000_0000, 3, 0x8000_0000);
    try expectRootF(0x3f80_0001, 2, 0x3f80_0000);

    try expectRootD(0x4020_0000_0000_0000, 3, 0x4000_0000_0000_0000);
    try expectRootD(0xc020_0000_0000_0000, 3, 0xc000_0000_0000_0000);
    try expectRootD(0x4010_0000_0000_0000, 2, 0x4000_0000_0000_0000);
    try expectRootD(0x4030_0000_0000_0000, -2, 0x3fd0_0000_0000_0000);
    try expectRootD(0x8000_0000_0000_0000, -3, 0xfff0_0000_0000_0000);
    try expectRootD(0x3ff0_0000_0000_0001, 2, 0x3ff0_0000_0000_0000);
    try expectRootH(0x6692, 3, 0x49f2);
    try expectRootF(0x77b8_902c, 3, 0x5236_3168);
    try expectRootD(0x2889_000d_ce03_d28f, 3, 0x3822_9100_14a0_721d);
}

test "binary64 rootn returns the certified estimate bracket for every degree" {
    // Each expected result comes from an exact half-ulp comparison of the
    // candidate power with the input, so a proposal that the exact comparisons
    // cannot certify must fall back to the search and still round once.
    try expectRootD(0x3ff0_0000_0000_0001, 64, 0x3ff0_0000_0000_0000);
    try expectRootD(0x3ff0_0000_0000_0002, 64, 0x3ff0_0000_0000_0000);
    try expectRootD(0x3fef_ffff_ffff_ffff, 64, 0x3ff0_0000_0000_0000);
    try expectRootD(0x3ff8_0000_0000_0000, 3, 0x3ff2_50bf_e1b0_82f5);
    try expectRootD(0x4000_0000_0000_0000, 64, 0x3ff0_2c9a_3e77_8061);
    try expectRootD(0x4000_0000_0000_0000, 2, 0x3ff6_a09e_667f_3bcd);
    try expectRootD(0x4008_0000_0000_0000, 3, 0x3ff7_1374_4912_3ef6);
    try expectRootD(0x4008_0000_0000_0000, -63, 0x3fef_7262_465b_163e);
    try expectRootD(0x4000_0000_0000_0000, -64, 0x3fef_a7c1_819e_90d8);
    try expectRootD(0x3ff8_0000_0000_0000, -3, 0x3feb_f45f_04ce_f0b9);
    try expectRootD(0x7fef_ffff_ffff_ffff, 64, 0x40f0_0000_0000_0000);
    try expectRootD(0x0000_0000_0000_0001, 64, 0x3ee2_9e9d_f51f_dee1);
    try expectRootD(0x0000_0000_0000_0001, -64, 0x40fb_7f76_f2fb_5e47);
    try expectRootD(0x0000_0000_0000_0001, -2, 0x6180_0000_0000_0000);
}

test "degree-one roots preserve finite bits and negative-one roots remain reciprocal" {
    try expectRootH(0x0001, 1, 0x0001);
    try expectRootH(0xbc01, 1, 0xbc01);
    try expectRootH(0xc000, -1, 0xb800);
    try expectRootF(0x0000_0001, 1, 0x0000_0001);
    try expectRootF(0xbf80_0001, 1, 0xbf80_0001);
    try expectRootF(0xc000_0000, -1, 0xbf00_0000);
    try expectRootD(0x0000_0000_0000_0001, 1, 0x0000_0000_0000_0001);
    try expectRootD(0xbff0_0000_0000_0001, 1, 0xbff0_0000_0000_0001);
    try expectRootD(0xc000_0000_0000_0000, -1, 0xbfe0_0000_0000_0000);
    try expectRootH(0xbc00, -63, 0xbc00);
    try expectRootF(0x3f80_0000, -64, 0x3f80_0000);
    try expectRootD(0xbff0_0000_0000_0000, 63, 0xbff0_0000_0000_0000);
}

test "reciprocal roots bracket the subnormal and overflow endpoints" {
    try expectRootH(0x7bff, -1, 0x0100);
    try expectRootH(0x0001, -1, 0x7c00);
    try expectRootF(0x7f7f_ffff, -1, 0x0020_0000);
    try expectRootF(0x0000_0001, -1, 0x7f80_0000);
    try expectRootD(0x7fef_ffff_ffff_ffff, -1, 0x0004_0000_0000_0000);
    try expectRootD(0x0000_0000_0000_0001, -1, 0x7ff0_0000_0000_0000);
}

test "rootn preserves exceptional values and exact degree bound" {
    try expectRootH(0xc000, 2, 0x7e00);
    try expectRootH(0x7c00, 3, 0x7c00);
    try expectRootH(0xfc00, -3, 0x8000);
    try expectRootH(0x3c00, 64, 0x3c00);
    try expectRootH(0x3c00, -64, 0x3c00);

    try expectRootF(0xc000_0000, 2, 0x7fc0_0000);
    try expectRootF(0x7f80_0000, -3, 0x0000_0000);
    try expectRootF(0x3f80_0000, 64, 0x3f80_0000);
    try expectRootF(0x3f80_0000, -64, 0x3f80_0000);

    try expectRootD(0xc000_0000_0000_0000, 2, 0x7ff8_0000_0000_0000);
    try expectRootD(0x7ff0_0000_0000_0000, -3, 0x0000_0000_0000_0000);
    try expectRootD(0x3ff0_0000_0000_0000, 64, 0x3ff0_0000_0000_0000);
    try expectRootD(0x3ff0_0000_0000_0000, -64, 0x3ff0_0000_0000_0000);

    const h_low = power.rootnh(common.H.fromBits(0x3c00), -65);
    const h_high = power.rootnh(common.H.fromBits(0x3c00), 65);
    try std.testing.expect(!h_low.ok);
    try std.testing.expect(!h_high.ok);
    try std.testing.expectEqual(@as(u16, 0x7e00), common.H.toBits(h_low.value));
    try std.testing.expectEqual(@as(u16, 0x7e00), common.H.toBits(h_high.value));

    const f_invalid = power.rootnf(1.0, 65);
    try std.testing.expect(!f_invalid.ok);
    try std.testing.expectEqual(@as(u32, 0x7fc0_0000), common.F.toBits(f_invalid.value));
    const d_invalid = power.rootnd(1.0, -65);
    try std.testing.expect(!d_invalid.ok);
    try std.testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), common.D.toBits(d_invalid.value));
}

test "binary16 power families preserve minimum subnormal boundaries at degree limits" {
    try expectPownH(0x0001, 64, 0x0000);
    try expectPownH(0x0001, -64, 0x7c00);
    try expectPownH(0x0002, 64, 0x0000);
    try expectPownH(0x0002, -64, 0x7c00);

    try expectCompoundH(0x0001, 64, 0x3c00);
    try expectCompoundH(0x0001, -64, 0x3c00);
    try expectCompoundH(0x0002, 64, 0x3c00);
    try expectCompoundH(0x0002, -64, 0x3c00);

    try expectRootH(0x0001, -64, 0x3d30);
    try expectRootH(0x0001, 64, 0x3a2b);
    try expectRootH(0x0002, -64, 0x3d22);
    try expectRootH(0x0002, 64, 0x3a3c);
}

test "binary32 power families preserve minimum subnormal boundaries at degree limits" {
    try expectPownF(0x0000_0001, 64, 0x0000_0000);
    try expectPownF(0x0000_0001, -64, 0x7f80_0000);
    try expectPownF(0x0000_0002, 64, 0x0000_0000);
    try expectPownF(0x0000_0002, -64, 0x7f80_0000);

    try expectCompoundF(0x0000_0001, 64, 0x3f80_0000);
    try expectCompoundF(0x0000_0001, -64, 0x3f80_0000);
    try expectCompoundF(0x0000_0002, 64, 0x3f80_0000);
    try expectCompoundF(0x0000_0002, -64, 0x3f80_0000);

    try expectRootF(0x0000_0001, -64, 0x40a0_b051);
    try expectRootF(0x0000_0001, 64, 0x3e4b_ec15);
    try expectRootF(0x0000_0002, -64, 0x409e_f532);
    try expectRootF(0x0000_0002, 64, 0x3e4e_248c);
}

test "binary64 power families preserve minimum subnormal boundaries at degree limits" {
    try expectPownD(0x0000_0000_0000_0001, 64, 0x0000_0000_0000_0000);
    try expectPownD(0x0000_0000_0000_0001, -64, 0x7ff0_0000_0000_0000);
    try expectPownD(0x0000_0000_0000_0002, 64, 0x0000_0000_0000_0000);
    try expectPownD(0x0000_0000_0000_0002, -64, 0x7ff0_0000_0000_0000);

    try expectCompoundD(0x0000_0000_0000_0001, 64, 0x3ff0_0000_0000_0000);
    try expectCompoundD(0x0000_0000_0000_0001, -64, 0x3ff0_0000_0000_0000);
    try expectCompoundD(0x0000_0000_0000_0002, 64, 0x3ff0_0000_0000_0000);
    try expectCompoundD(0x0000_0000_0000_0002, -64, 0x3ff0_0000_0000_0000);

    try expectRootD(0x0000_0000_0000_0001, -64, 0x40fb_7f76_f2fb_5e47);
    try expectRootD(0x0000_0000_0000_0001, 64, 0x3ee2_9e9d_f51f_dee1);
    try expectRootD(0x0000_0000_0000_0002, -64, 0x40fb_33a2_b84f_15fb);
    try expectRootD(0x0000_0000_0000_0002, 64, 0x3ee2_d285_a6e4_030b);
}

test "binary64 power fast paths reciprocal degree minus one rounds near one" {
    // 1/(1+2^-52) is above the double 1-2^-52 and below the midpoint
    // above it, while 1/(1-2^-53) is above the midpoint above one,
    // so neither case resolves a tie.
    try expectPownD(0x3ff0_0000_0000_0001, -1, 0x3fef_ffff_ffff_fffe);
    try expectPownD(0x3fef_ffff_ffff_ffff, -1, 0x3ff0_0000_0000_0001);
    try expectPownD(0xbfef_ffff_ffff_ffff, -1, 0xbff0_0000_0000_0001);
    try expectPownD(0xbff0_0000_0000_0000, -1, 0xbff0_0000_0000_0000);

    try expectRootD(0x3ff0_0000_0000_0001, -1, 0x3fef_ffff_ffff_fffe);
    try expectRootD(0x3fef_ffff_ffff_ffff, -1, 0x3ff0_0000_0000_0001);
    try expectRootD(0xbfef_ffff_ffff_ffff, -1, 0xbff0_0000_0000_0001);
}

test "binary64 power fast paths reciprocal degree minus one keeps range edges" {
    // The reciprocal of the smallest subnormal is 2^1074 and overflows
    // to infinity in both signs.  The largest subnormal reciprocal is
    // above 2^1022 + 2^970 but below the midpoint above it, and the
    // minimum normal reciprocal is the exact power 2^1022.
    try expectPownD(0x0000_0000_0000_0001, -1, 0x7ff0_0000_0000_0000);
    try expectPownD(0x8000_0000_0000_0001, -1, 0xfff0_0000_0000_0000);
    try expectPownD(0x000f_ffff_ffff_ffff, -1, 0x7fd0_0000_0000_0001);
    try expectPownD(0x0010_0000_0000_0000, -1, 0x7fd0_0000_0000_0000);
    // The largest finite value is (2 - 2^-52) * 2^1023, and its
    // reciprocal is above the subnormal 2^-1024 but rounds to it.
    try expectPownD(0x7fef_ffff_ffff_ffff, -1, 0x0004_0000_0000_0000);
    // 2^1022 is the largest base whose reciprocal is still normal, and it
    // lands exactly on the minimum normal.  The next base above it has the
    // largest subnormal reciprocal, since 2^104/(2^52+1) is 2^52-1 plus a
    // fraction too small to round up, and 2^1023 gives the exact 2^51 ulps.
    try expectPownD(0x7fd0_0000_0000_0000, -1, 0x0010_0000_0000_0000);
    try expectPownD(0x7fd0_0000_0000_0001, -1, 0x000f_ffff_ffff_ffff);
    try expectPownD(0x7fe0_0000_0000_0000, -1, 0x0008_0000_0000_0000);
    try expectRootD(0x7fd0_0000_0000_0000, -1, 0x0010_0000_0000_0000);
    try expectRootD(0x7fd0_0000_0000_0001, -1, 0x000f_ffff_ffff_ffff);
    try expectRootD(0x7fe0_0000_0000_0000, -1, 0x0008_0000_0000_0000);
    // The sign follows the base, and the negative base takes the same path.
    try expectPownD(0xffef_ffff_ffff_ffff, -1, 0x8004_0000_0000_0000);
    try expectRootD(0xffef_ffff_ffff_ffff, -1, 0x8004_0000_0000_0000);

    // Zero and infinity take their magnitude from the degree and their
    // sign from the base sign and the degree parity.
    try expectPownD(0x0000_0000_0000_0000, -1, 0x7ff0_0000_0000_0000);
    try expectPownD(0x8000_0000_0000_0000, -1, 0xfff0_0000_0000_0000);
    try expectPownD(0x7ff0_0000_0000_0000, -1, 0x0000_0000_0000_0000);
    try expectPownD(0xfff0_0000_0000_0000, -1, 0x8000_0000_0000_0000);
    try expectRootD(0x8000_0000_0000_0000, -1, 0xfff0_0000_0000_0000);
    try expectRootD(0x7ff0_0000_0000_0000, -1, 0x0000_0000_0000_0000);
    try expectRootD(0xfff0_0000_0000_0000, -1, 0x8000_0000_0000_0000);
}

test "binary64 power fast paths degree two sqrt rounds half ulp boundaries" {
    // sqrt(1+3*2^-52) sits just below the midpoint between 1+2^-52
    // and 1+2^-51, and sqrt(1-2^-53) sits just below the midpoint
    // below one, so neither resolves a tie.
    try expectRootD(0x3ff0_0000_0000_0003, 2, 0x3ff0_0000_0000_0001);
    try expectRootD(0x3fef_ffff_ffff_ffff, 2, 0x3fef_ffff_ffff_ffff);
    // The smallest subnormal root is the exact power 2^-537, and the
    // largest subnormal root is one ulp under the exact power 2^-511.
    try expectRootD(0x0000_0000_0000_0001, 2, 0x1e60_0000_0000_0000);
    try expectRootD(0x000f_ffff_ffff_ffff, 2, 0x1fff_ffff_ffff_ffff);
    // The root of the largest finite value lies just below the
    // midpoint 2^512 - 2^458 and rounds to the predecessor of 2^512.
    try expectRootD(0x7fef_ffff_ffff_ffff, 2, 0x5fef_ffff_ffff_ffff);
}

test "binary64 power fast paths degree one compound rounds one plus x once" {
    // Cancellation near -1 keeps every bit of the exact sum.
    try expectCompoundD(0xbfef_ffff_ffff_ffff, 1, 0x3ca0_0000_0000_0000);
    try expectCompoundD(0xbfef_ffff_ffff_fffe, 1, 0x3cb0_0000_0000_0000);
    // 1 + 2^-53 is the midpoint above one and takes the even one,
    // 1 + 2^-54 is only a quarter ulp up so it stays one, and
    // 1 - 2^-54 is the midpoint below one, which also takes one.
    try expectCompoundD(0x3ca0_0000_0000_0000, 1, 0x3ff0_0000_0000_0000);
    try expectCompoundD(0x3c90_0000_0000_0000, 1, 0x3ff0_0000_0000_0000);
    try expectCompoundD(0xbc90_0000_0000_0000, 1, 0x3ff0_0000_0000_0000);
    // 1 + 3*2^-53 is the midpoint between 1+2^-52 and 1+2^-51 and
    // takes the even mantissa 2, while 1 - 3*2^-54 is the midpoint
    // between 1-2^-53 and 1-2^-52 and takes the even 1-2^-52.
    try expectCompoundD(0x3cb8_0000_0000_0000, 1, 0x3ff0_0000_0000_0002);
    try expectCompoundD(0xbca8_0000_0000_0000, 1, 0x3fef_ffff_ffff_fffe);
    // 1 + (1 + 2^-52) is a midpoint above two and ties to two, and
    // negative zero adds nothing.
    try expectCompoundD(0x3ff0_0000_0000_0001, 1, 0x4000_0000_0000_0000);
    try expectCompoundD(0x8000_0000_0000_0000, 1, 0x3ff0_0000_0000_0000);
}

test "binary64 power fast paths NaN and domain rules stay ahead of the fast paths" {
    try expectPownD(0x7ff0_0000_0000_0123, -1, 0x7ff8_0000_0000_0123);
    try expectPownD(0xfff8_0000_0000_0456, -1, 0xfff8_0000_0000_0456);
    try expectRootD(0x7ff0_0000_0000_0123, -1, 0x7ff8_0000_0000_0123);
    try expectRootD(0x7ff0_0000_0000_0123, 2, 0x7ff8_0000_0000_0123);
    try expectRootD(0xfff8_0000_0000_0456, 2, 0xfff8_0000_0000_0456);
    try expectCompoundD(0x7ff0_0000_0000_0123, 1, 0x7ff8_0000_0000_0123);

    // A negative nonzero value with an even degree is a domain error.
    // Negative zero is exempt because it is not below zero, but
    // negative infinity is a negative nonzero value, so the domain
    // answer precedes the infinity rule.
    try expectRootD(0xbff0_0000_0000_0000, 2, 0x7ff8_0000_0000_0000);
    try expectRootD(0x8000_0000_0000_0000, 2, 0x0000_0000_0000_0000);
    try expectRootD(0xfff0_0000_0000_0000, 2, 0x7ff8_0000_0000_0000);

    // Every non-NaN value below -1 is outside the compoundn domain,
    // -1 has its own zero result, and positive infinity stays infinite.
    try expectCompoundD(0xc000_0000_0000_0000, 1, 0x7ff8_0000_0000_0000);
    try expectCompoundD(0xfff0_0000_0000_0000, 1, 0x7ff8_0000_0000_0000);
    try expectCompoundD(0xbff0_0000_0000_0000, 1, 0x0000_0000_0000_0000);
    try expectCompoundD(0x7ff0_0000_0000_0000, 1, 0x7ff0_0000_0000_0000);
}

test "binary64 power fast paths out of range degrees stay ahead of the fast paths" {
    const canonical: u64 = 0x7ff8_0000_0000_0000;

    // Degrees outside the accepted -64 through 64 range report the
    // canonical NaN sentinel with a false status.  Degree 65 is odd,
    // so an in-range root of -1 would answer -1 instead.
    const pown_snan = power.pownd(common.D.fromBits(0x7ff0_0000_0000_0123), -65);
    try std.testing.expect(!pown_snan.ok);
    try std.testing.expectEqual(canonical, common.D.toBits(pown_snan.value));
    const pown_zero = power.pownd(common.D.fromBits(0x0000_0000_0000_0000), 65);
    try std.testing.expect(!pown_zero.ok);
    try std.testing.expectEqual(canonical, common.D.toBits(pown_zero.value));

    const root_odd_degree = power.rootnd(common.D.fromBits(0xbff0_0000_0000_0000), 65);
    try std.testing.expect(!root_odd_degree.ok);
    try std.testing.expectEqual(canonical, common.D.toBits(root_odd_degree.value));
    const root_negative_degree = power.rootnd(common.D.fromBits(0x8000_0000_0000_0000), -65);
    try std.testing.expect(!root_negative_degree.ok);
    try std.testing.expectEqual(canonical, common.D.toBits(root_negative_degree.value));

    const compound_domain = power.compoundnd(common.D.fromBits(0xc000_0000_0000_0000), -65);
    try std.testing.expect(!compound_domain.ok);
    try std.testing.expectEqual(canonical, common.D.toBits(compound_domain.value));
    const compound_snan = power.compoundnd(common.D.fromBits(0x7ff0_0000_0000_0123), 65);
    try std.testing.expect(!compound_snan.ok);
    try std.testing.expectEqual(canonical, common.D.toBits(compound_snan.value));
}

test "binary64 power fast paths degree minus two keeps one exact root" {
    // 2^(-1/2) sits just above the midpoint between
    // 0x3fe6_a09e_667f_3bcc and 0x3fe6_a09e_667f_3bcd, so the once
    // rounded root is the upper neighbour.
    try expectRootD(0x4000_0000_0000_0000, -2, 0x3fe6_a09e_667f_3bcd);
    // sqrt(1 + 2^-52) rounds to one and 1/1 is one, but the exact
    // root is above 1-2^-53 and rounds to it.
    try expectRootD(0x3ff0_0000_0000_0001, -2, 0x3fef_ffff_ffff_ffff);
    // sqrt(2^600 + 2^548) rounds to 2^300 and 1/2^300 is 2^-300, but
    // the exact root is just above 2^-300 - 2^-353 and rounds to it.
    try expectRootD(0x6570_0000_0000_0001, -2, 0x2d2f_ffff_ffff_ffff);
}

test "binary64 exact power optimization positive powers tie to even significands" {
    // (3/2)^34 is 3^34 / 2^34 and 3^34 needs 54 significant bits, so the
    // exact value is precisely the midpoint of 0x412d_9fe7_7988_1944 and
    // 0x412d_9fe7_7988_1945.  The lower neighbour holds the even
    // significand, so ties to even keeps that one.
    try expectPownD(0x3ff8_0000_0000_0000, 34, 0x412d_9fe7_7988_1944);
    // (7/4)^19 is 7^19 / 2^38 and 7^19 also needs 54 significant bits, so
    // this is a second exact midpoint, here between an odd lower neighbour
    // and an even upper one.
    try expectPownD(0x3ffc_0000_0000_0000, 19, 0x40e4_3f9e_0d2d_93ec);
}

test "binary64 exact power optimization truncation keeps guard and sticky bits" {
    // (5/4)^41 is 5^41 / 2^82 with 96 significant bits, and the discarded
    // remainder sits above the midpoint, so the upper neighbour wins.
    try expectPownD(0x3ff4_0000_0000_0000, 41, 0x40c2_5dfa_371a_19e7);
    // 1.6 already fills all 53 significand bits, so its degree-20 power
    // needs about a thousand of them and again lands above a midpoint.
    try expectPownD(0x3ff9_9999_9999_999a, 20, 0x40c7_9ca1_0c92_422b);
    // (1 + 2^-52)^63 is not exactly 1 + 63*2^-52; its retained fraction
    // lies below half an ulp, so the result truncates to 63.
    try expectPownD(0x3ff0_0000_0000_0001, 63, 0x3ff0_0000_0000_003f);
    // The same base at degree 19 keeps a retained fraction below half an
    // ulp as well.
    try expectPownD(0x3ff0_0000_0000_0001, 19, 0x3ff0_0000_0000_0013);
}

test "binary64 exact power optimization subnormal powers round to whole ulps" {
    // 1.5 * 2^-358 cubed is exactly 3.375 minimum subnormal ulps.
    try expectPownD(0x2998_0000_0000_0000, 3, 0x0000_0000_0000_0003);
    // 1.75 * 2^-359 cubed is 0.669921875 of a minimum subnormal ulp, which
    // is past the half-ulp midpoint and so keeps one ulp.
    try expectPownD(0x298c_0000_0000_0000, 3, 0x0000_0000_0000_0001);
}

test "binary64 exact power optimization negative odd degrees keep magnitude rounding" {
    // The negative base reuses the degree-19 midpoint; the magnitude rounds
    // exactly as the positive base does and only the sign follows parity.
    try expectPownD(0xbffc_0000_0000_0000, 19, 0xc0e4_3f9e_0d2d_93ec);
    try expectPownD(0xbff8_0000_0000_0000, 35, 0xc136_37ed_9b26_12f3);
    // The subnormal rows keep their single-ulp significands under negation.
    try expectPownD(0xa998_0000_0000_0000, 3, 0x8000_0000_0000_0003);
    try expectPownD(0xa98c_0000_0000_0000, 3, 0x8000_0000_0000_0001);
}

test "binary64 exact power optimization compound powers expand once before rounding" {
    // (1 + 2^-54)^2 is 1 + 2^-53 + 2^-108, strictly above the midpoint
    // 1 + 2^-53, so it rounds up.  Rounding 1 + x first collapses it to
    // one, and squaring that gives one.
    try expectCompoundD(0x3c90_0000_0000_0000, 2, 0x3ff0_0000_0000_0001);
    // At degree 21 the same base retains a fraction below half an ulp
    // above 1 + 5.25 ulps, giving five ulps.
    try expectCompoundD(0x3c90_0000_0000_0000, 21, 0x3ff0_0000_0000_0005);
    try expectCompoundD(0x3ca0_0000_0000_0000, 19, 0x3ff0_0000_0000_000a);
    // (2 + 2^-52)^19 is 2^19 + 9.5 ulps plus a positive second-order
    // term, so it clears the midpoint to the even upper neighbour.
    try expectCompoundD(0x3ff0_0000_0000_0001, 19, 0x4120_0000_0000_000a);
    // Here 1 + x is exactly 2^-53, so the degree-20 value is the exact
    // power 2^-1060 and the degree-21 value is the exact power 2^-1113,
    // which is below half the minimum subnormal.
    try expectCompoundD(0xbfef_ffff_ffff_ffff, 20, 0x0000_0000_0000_4000);
    try expectCompoundD(0xbfef_ffff_ffff_ffff, 21, 0x0000_0000_0000_0000);
}

test "binary64 exact power optimization rootn stays exact across wide exponent offsets" {
    // The cube root of the minimum subnormal is the exact power 2^-358.
    try expectRootD(0x0000_0000_0000_0001, 3, 0x2990_0000_0000_0000);
    // The minimum normal is 2^-1022, so its square root and its inverse
    // square root are the exact powers 2^-511 and 2^511.
    try expectRootD(0x0010_0000_0000_0000, 2, 0x2000_0000_0000_0000);
    try expectRootD(0x0010_0000_0000_0000, -2, 0x5fe0_0000_0000_0000);
    // The cube roots of the minimum normal and of the largest finite value
    // share one significand, because the 2^-52 shortfall of the largest
    // finite value is far below half an ulp of its cube root.
    try expectRootD(0x0010_0000_0000_0000, 3, 0x2aa4_28a2_f98d_728b);
    try expectRootD(0x7fef_ffff_ffff_ffff, 3, 0x5544_28a2_f98d_728b);
    // The negative degrees of those same two bases are both inexact.
    try expectRootD(0x0010_0000_0000_0000, -3, 0x5539_65fe_a53d_6e3d);
    try expectRootD(0x7fef_ffff_ffff_ffff, -3, 0x2a99_65fe_a53d_6e3d);
    // Degree -63 inverts the exponent scale, putting the largest finite
    // value near 2^-16 and the largest subnormal near 2^16, and neither
    // root is exact.
    try expectRootD(0x7fef_ffff_ffff_ffff, -63, 0x3eea_d5b3_a4a1_6c79);
    try expectRootD(0x000f_ffff_ffff_ffff, -63, 0x40f2_aa1a_5aad_04f5);
    // The exact cube root of 2^600 + 2^548 lies between a quarter and a
    // half of an ulp above 2^200, so it rounds back down to 2^200.
    try expectRootD(0x6570_0000_0000_0001, 3, 0x4c70_0000_0000_0000);
}
