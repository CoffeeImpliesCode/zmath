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
    const pown_d_snan = power.pownd(common.D.fromBits(0x7ff0_0000_0000_0123), 4);
    try std.testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0123), common.D.toBits(pown_d_snan.value));

    const root_h_zero = power.rootnh(common.H.fromBits(0x7e55), 0);
    try std.testing.expect(root_h_zero.ok);
    try std.testing.expectEqual(@as(u16, 0x7e00), common.H.toBits(root_h_zero.value));
    const root_f_zero = power.rootnf(common.F.fromBits(0x7fc0_0123), 0);
    try std.testing.expectEqual(@as(u32, 0x7fc0_0000), common.F.toBits(root_f_zero.value));
    const root_d_zero = power.rootnd(common.D.fromBits(0x7ff0_0000_0000_0123), 0);
    try std.testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), common.D.toBits(root_d_zero.value));
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
    try expectCompoundH(0x7c00, 3, 0x7c00);
    try expectCompoundH(0x7c00, -3, 0x0000);

    try expectCompoundF(0x3380_0000, 2, 0x3f80_0001);
    try expectCompoundF(0x3f80_0000, 64, 0x5f80_0000);
    try expectCompoundF(0xbf80_0000, 3, 0x0000_0000);

    try expectCompoundD(0x3ca0_0000_0000_0000, 2, 0x3ff0_0000_0000_0001);
    try expectCompoundD(0x3ff0_0000_0000_0000, -64, 0x3bf0_0000_0000_0000);
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
