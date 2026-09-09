const std = @import("std");
const pi = @import("scalar_pi_inverse.zig");
const c23 = @import("scalar_c23_types.zig");

fn expectH(expected_bits: u16, actual: f16) !void {
    try std.testing.expectEqual(expected_bits, @as(u16, @bitCast(actual)));
}

fn expectFBits(expected_bits: u32, actual: f32) !void {
    try std.testing.expectEqual(expected_bits, @as(u32, @bitCast(actual)));
}

fn expectDBits(expected_bits: u64, actual: f64) !void {
    try std.testing.expectEqual(expected_bits, @as(u64, @bitCast(actual)));
}

fn expectAtan2H(expected_bits: u16, actual: c23.Atan2Result(f16)) !void {
    try std.testing.expect(actual.resolved);
    try expectH(expected_bits, actual.value);
}

fn expectAtan2F(expected_bits: u32, actual: c23.Atan2Result(f32)) !void {
    try std.testing.expect(actual.resolved);
    try expectFBits(expected_bits, actual.value);
}

fn expectAtan2D(expected_bits: u64, actual: c23.Atan2Result(f64)) !void {
    try std.testing.expect(actual.resolved);
    try expectDBits(expected_bits, actual.value);
}

fn h(bits: u16) f16 {
    return @bitCast(bits);
}

fn f(bits: u32) f32 {
    return @bitCast(bits);
}

fn d(bits: u64) f64 {
    return @bitCast(bits);
}

test "inverse pi exports exact H F D signatures" {
    comptime {
        std.debug.assert(@TypeOf(pi.asinpih) == fn (f16) f16);
        std.debug.assert(@TypeOf(pi.asinpif) == fn (f32) f32);
        std.debug.assert(@TypeOf(pi.asinpid) == fn (f64) f64);
        std.debug.assert(@TypeOf(pi.acospih) == fn (f16) f16);
        std.debug.assert(@TypeOf(pi.acospif) == fn (f32) f32);
        std.debug.assert(@TypeOf(pi.acospid) == fn (f64) f64);
        std.debug.assert(@TypeOf(pi.atanpih) == fn (f16) f16);
        std.debug.assert(@TypeOf(pi.atanpif) == fn (f32) f32);
        std.debug.assert(@TypeOf(pi.atanpid) == fn (f64) f64);
        std.debug.assert(@TypeOf(pi.atan2pih) == fn (f16, f16) c23.Atan2Result(f16));
        std.debug.assert(@TypeOf(pi.atan2pif) == fn (f32, f32) c23.Atan2Result(f32));
        std.debug.assert(@TypeOf(pi.atan2pid) == fn (f64, f64) c23.Atan2Result(f64));
    }
}

test "inverse pi finite kernels match pinned upstream hard cases" {
    try expectH(0x321e, pi.asinpih(h(0x3885)));
    try expectH(0x37ef, pi.acospih(h(0x22ad)));
    try expectH(0x3809, pi.acospih(h(0xa2ad)));
    try expectH(0x3537, pi.atanpih(h(0x3e92)));

    try expectFBits(0x3e93_3802, pi.atanpif(f(0x3fa2_67dd)));
    try expectFBits(0x3e70_d331, pi.atanpif(f(0x3f69_3531)));

    try expectDBits(0x3d77_fe08_f528_4726, pi.asinpid(d(0x3d92_d7f1_6888_8139)));
    try expectDBits(0x3fe2_2530_cb92_c8f1, pi.acospid(d(0xbfca_c26c_c49b_0264)));
    try expectDBits(0x3d90_5d3a_a746_e163, pi.atanpid(d(0x3da9_b46c_fb0f_10ca)));
}

test "asinpi preserves signed zero and exact domain endpoints" {
    try expectH(0x0000, pi.asinpih(h(0x0000)));
    try expectH(0x8000, pi.asinpih(h(0x8000)));
    try expectH(0x3800, pi.asinpih(h(0x3c00)));
    try expectH(0xb800, pi.asinpih(h(0xbc00)));

    try expectFBits(0x0000_0000, pi.asinpif(f(0x0000_0000)));
    try expectFBits(0x8000_0000, pi.asinpif(f(0x8000_0000)));
    try expectFBits(0x3f00_0000, pi.asinpif(f(0x3f80_0000)));
    try expectFBits(0xbf00_0000, pi.asinpif(f(0xbf80_0000)));

    try expectDBits(0x0000_0000_0000_0000, pi.asinpid(d(0x0000_0000_0000_0000)));
    try expectDBits(0x8000_0000_0000_0000, pi.asinpid(d(0x8000_0000_0000_0000)));
    try expectDBits(0x3fe0_0000_0000_0000, pi.asinpid(d(0x3ff0_0000_0000_0000)));
    try expectDBits(0xbfe0_0000_0000_0000, pi.asinpid(d(0xbff0_0000_0000_0000)));
}

test "asinpi rejects invalid finite values and infinities" {
    try expectH(0x7e00, pi.asinpih(h(0x4000)));
    try expectH(0x7e00, pi.asinpih(h(0xc000)));
    try expectH(0x7e00, pi.asinpih(h(0x7c00)));
    try expectH(0x7e00, pi.asinpih(h(0xfc00)));

    try expectFBits(0x7fc0_0000, pi.asinpif(f(0x4000_0000)));
    try expectFBits(0x7fc0_0000, pi.asinpif(f(0xc000_0000)));
    try expectFBits(0x7fc0_0000, pi.asinpif(f(0x7f80_0000)));
    try expectFBits(0x7fc0_0000, pi.asinpif(f(0xff80_0000)));

    try expectDBits(0x7ff8_0000_0000_0000, pi.asinpid(d(0x4000_0000_0000_0000)));
    try expectDBits(0x7ff8_0000_0000_0000, pi.asinpid(d(0xc000_0000_0000_0000)));
    try expectDBits(0x7ff8_0000_0000_0000, pi.asinpid(d(0x7ff0_0000_0000_0000)));
    try expectDBits(0x7ff8_0000_0000_0000, pi.asinpid(d(0xfff0_0000_0000_0000)));
}

test "asinpi quiets NaNs in parameter order" {
    try expectH(0x7f55, pi.asinpih(h(0x7d55)));
    try expectH(0xff55, pi.asinpih(h(0xfd55)));
    try expectFBits(0x7fc1_2345, pi.asinpif(f(0x7f81_2345)));
    try expectFBits(0xffc1_2345, pi.asinpif(f(0xff81_2345)));
    try expectDBits(0x7ff8_0000_0000_1234, pi.asinpid(d(0x7ff0_0000_0000_1234)));
    try expectDBits(0xfff8_0000_0000_1234, pi.asinpid(d(0xfff0_0000_0000_1234)));
}

test "acospi preserves exact endpoints and maps both zeros to one half" {
    try expectH(0x3800, pi.acospih(h(0x0000)));
    try expectH(0x3800, pi.acospih(h(0x8000)));
    try expectH(0x0000, pi.acospih(h(0x3c00)));
    try expectH(0x3c00, pi.acospih(h(0xbc00)));

    try expectFBits(0x3f00_0000, pi.acospif(f(0x0000_0000)));
    try expectFBits(0x3f00_0000, pi.acospif(f(0x8000_0000)));
    try expectFBits(0x0000_0000, pi.acospif(f(0x3f80_0000)));
    try expectFBits(0x3f80_0000, pi.acospif(f(0xbf80_0000)));

    try expectDBits(0x3fe0_0000_0000_0000, pi.acospid(d(0x0000_0000_0000_0000)));
    try expectDBits(0x3fe0_0000_0000_0000, pi.acospid(d(0x8000_0000_0000_0000)));
    try expectDBits(0x0000_0000_0000_0000, pi.acospid(d(0x3ff0_0000_0000_0000)));
    try expectDBits(0x3ff0_0000_0000_0000, pi.acospid(d(0xbff0_0000_0000_0000)));
}

test "acospi rejects invalid values and preserves first NaN" {
    try expectH(0x7e00, pi.acospih(h(0x4000)));
    try expectH(0x7e00, pi.acospih(h(0x7c00)));
    try expectFBits(0x7fc0_0000, pi.acospif(f(0xc000_0000)));
    try expectFBits(0x7fc0_0000, pi.acospif(f(0xff80_0000)));
    try expectDBits(0x7ff8_0000_0000_0000, pi.acospid(d(0x4000_0000_0000_0000)));
    try expectDBits(0x7ff8_0000_0000_0000, pi.acospid(d(0xfff0_0000_0000_0000)));

    try expectH(0xff55, pi.acospih(h(0xfd55)));
    try expectFBits(0xffc1_2345, pi.acospif(f(0xff81_2345)));
    try expectDBits(0x7ff8_0000_0000_1234, pi.acospid(d(0x7ff0_0000_0000_1234)));
}

test "atanpi preserves zero signs and exact branch and infinity values" {
    try expectH(0x0000, pi.atanpih(h(0x0000)));
    try expectH(0x8000, pi.atanpih(h(0x8000)));
    try expectH(0x3400, pi.atanpih(h(0x3c00)));
    try expectH(0xb400, pi.atanpih(h(0xbc00)));
    try expectH(0x3800, pi.atanpih(h(0x7c00)));
    try expectH(0xb800, pi.atanpih(h(0xfc00)));

    try expectFBits(0x0000_0000, pi.atanpif(f(0x0000_0000)));
    try expectFBits(0x8000_0000, pi.atanpif(f(0x8000_0000)));
    try expectFBits(0x3e80_0000, pi.atanpif(f(0x3f80_0000)));
    try expectFBits(0xbe80_0000, pi.atanpif(f(0xbf80_0000)));
    try expectFBits(0x3f00_0000, pi.atanpif(f(0x7f80_0000)));
    try expectFBits(0xbf00_0000, pi.atanpif(f(0xff80_0000)));

    try expectDBits(0x0000_0000_0000_0000, pi.atanpid(d(0x0000_0000_0000_0000)));
    try expectDBits(0x8000_0000_0000_0000, pi.atanpid(d(0x8000_0000_0000_0000)));
    try expectDBits(0x3fd0_0000_0000_0000, pi.atanpid(d(0x3ff0_0000_0000_0000)));
    try expectDBits(0xbfd0_0000_0000_0000, pi.atanpid(d(0xbff0_0000_0000_0000)));
    try expectDBits(0x3fe0_0000_0000_0000, pi.atanpid(d(0x7ff0_0000_0000_0000)));
    try expectDBits(0xbfe0_0000_0000_0000, pi.atanpid(d(0xfff0_0000_0000_0000)));
}

test "atanpi quiets every NaN format without changing payload" {
    try expectH(0x7f55, pi.atanpih(h(0x7d55)));
    try expectFBits(0xffc1_2345, pi.atanpif(f(0xff81_2345)));
    try expectDBits(0x7ff8_0000_0000_1234, pi.atanpid(d(0x7ff0_0000_0000_1234)));
}

test "atan2pi preserves signed-zero axes in every quadrant" {
    try expectAtan2H(0x0000, pi.atan2pih(h(0x0000), h(0x0000)));
    try expectAtan2H(0x8000, pi.atan2pih(h(0x8000), h(0x0000)));
    try expectAtan2H(0x3c00, pi.atan2pih(h(0x0000), h(0x8000)));
    try expectAtan2H(0xbc00, pi.atan2pih(h(0x8000), h(0x8000)));
    try expectAtan2H(0x3800, pi.atan2pih(h(0x3c00), h(0x0000)));
    try expectAtan2H(0xb800, pi.atan2pih(h(0xbc00), h(0x8000)));

    try expectAtan2F(0x0000_0000, pi.atan2pif(f(0x0000_0000), f(0x0000_0000)));
    try expectAtan2F(0x8000_0000, pi.atan2pif(f(0x8000_0000), f(0x0000_0000)));
    try expectAtan2F(0x3f80_0000, pi.atan2pif(f(0x0000_0000), f(0x8000_0000)));
    try expectAtan2F(0xbf80_0000, pi.atan2pif(f(0x8000_0000), f(0x8000_0000)));
    try expectAtan2F(0x3f00_0000, pi.atan2pif(f(0x3f80_0000), f(0x0000_0000)));
    try expectAtan2F(0xbf00_0000, pi.atan2pif(f(0xbf80_0000), f(0x8000_0000)));

    try expectAtan2D(
        0x0000_0000_0000_0000,
        pi.atan2pid(d(0x0000_0000_0000_0000), d(0x0000_0000_0000_0000)),
    );
    try expectAtan2D(
        0x8000_0000_0000_0000,
        pi.atan2pid(d(0x8000_0000_0000_0000), d(0x0000_0000_0000_0000)),
    );
    try expectAtan2D(
        0x3ff0_0000_0000_0000,
        pi.atan2pid(d(0x0000_0000_0000_0000), d(0x8000_0000_0000_0000)),
    );
    try expectAtan2D(
        0xbff0_0000_0000_0000,
        pi.atan2pid(d(0x8000_0000_0000_0000), d(0x8000_0000_0000_0000)),
    );
    try expectAtan2D(
        0x3fe0_0000_0000_0000,
        pi.atan2pid(d(0x3ff0_0000_0000_0000), d(0x0000_0000_0000_0000)),
    );
    try expectAtan2D(
        0xbfe0_0000_0000_0000,
        pi.atan2pid(d(0xbff0_0000_0000_0000), d(0x8000_0000_0000_0000)),
    );
}

test "atan2pi maps infinity quadrants and finite axes exactly" {
    try expectAtan2H(0x3400, pi.atan2pih(h(0x7c00), h(0x7c00)));
    try expectAtan2H(0x3a00, pi.atan2pih(h(0x7c00), h(0xfc00)));
    try expectAtan2H(0xb400, pi.atan2pih(h(0xfc00), h(0x7c00)));
    try expectAtan2H(0xba00, pi.atan2pih(h(0xfc00), h(0xfc00)));
    try expectAtan2H(0x3800, pi.atan2pih(h(0x7c00), h(0x3c00)));
    try expectAtan2H(0xb800, pi.atan2pih(h(0xfc00), h(0x3c00)));
    try expectAtan2H(0x0000, pi.atan2pih(h(0x3c00), h(0x7c00)));
    try expectAtan2H(0x8000, pi.atan2pih(h(0xbc00), h(0x7c00)));
    try expectAtan2H(0x3c00, pi.atan2pih(h(0x3c00), h(0xfc00)));
    try expectAtan2H(0xbc00, pi.atan2pih(h(0xbc00), h(0xfc00)));

    try expectAtan2F(0x3e80_0000, pi.atan2pif(f(0x7f80_0000), f(0x7f80_0000)));
    try expectAtan2F(0x3f40_0000, pi.atan2pif(f(0x7f80_0000), f(0xff80_0000)));
    try expectAtan2F(0xbe80_0000, pi.atan2pif(f(0xff80_0000), f(0x7f80_0000)));
    try expectAtan2F(0xbf40_0000, pi.atan2pif(f(0xff80_0000), f(0xff80_0000)));
    try expectAtan2F(0x3f00_0000, pi.atan2pif(f(0x7f80_0000), f(0x3f80_0000)));
    try expectAtan2F(0xbf00_0000, pi.atan2pif(f(0xff80_0000), f(0x3f80_0000)));
    try expectAtan2F(0x0000_0000, pi.atan2pif(f(0x3f80_0000), f(0x7f80_0000)));
    try expectAtan2F(0x8000_0000, pi.atan2pif(f(0xbf80_0000), f(0x7f80_0000)));
    try expectAtan2F(0x3f80_0000, pi.atan2pif(f(0x3f80_0000), f(0xff80_0000)));
    try expectAtan2F(0xbf80_0000, pi.atan2pif(f(0xbf80_0000), f(0xff80_0000)));

    try expectAtan2D(
        0x3fd0_0000_0000_0000,
        pi.atan2pid(d(0x7ff0_0000_0000_0000), d(0x7ff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0x3fe8_0000_0000_0000,
        pi.atan2pid(d(0x7ff0_0000_0000_0000), d(0xfff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0xbfd0_0000_0000_0000,
        pi.atan2pid(d(0xfff0_0000_0000_0000), d(0x7ff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0xbfe8_0000_0000_0000,
        pi.atan2pid(d(0xfff0_0000_0000_0000), d(0xfff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0x3fe0_0000_0000_0000,
        pi.atan2pid(d(0x7ff0_0000_0000_0000), d(0x3ff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0xbfe0_0000_0000_0000,
        pi.atan2pid(d(0xfff0_0000_0000_0000), d(0x3ff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0x0000_0000_0000_0000,
        pi.atan2pid(d(0x3ff0_0000_0000_0000), d(0x7ff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0x8000_0000_0000_0000,
        pi.atan2pid(d(0xbff0_0000_0000_0000), d(0x7ff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0x3ff0_0000_0000_0000,
        pi.atan2pid(d(0x3ff0_0000_0000_0000), d(0xfff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0xbff0_0000_0000_0000,
        pi.atan2pid(d(0xbff0_0000_0000_0000), d(0xfff0_0000_0000_0000)),
    );
}

test "atan2pi preserves first NaN in y then x order" {
    try expectAtan2H(0x7f01, pi.atan2pih(h(0x7d01), h(0xfc55)));
    try expectAtan2H(0xff55, pi.atan2pih(h(0x3c00), h(0xfd55)));
    try expectAtan2F(0x7fc0_00a5, pi.atan2pif(f(0x7f80_00a5), f(0xff80_0055)));
    try expectAtan2F(0xffc0_0055, pi.atan2pif(f(0x3f80_0000), f(0xff80_0055)));
    try expectAtan2D(0x7ff8_0000_0000_00a5, pi.atan2pid(
        d(0x7ff0_0000_0000_00a5),
        d(0xfff0_0000_0000_0055),
    ));
    try expectAtan2D(0xfff8_0000_0000_0055, pi.atan2pid(
        d(0x3ff0_0000_0000_0000),
        d(0xfff0_0000_0000_0055),
    ));
}

test "atan2pi returns exact principal branch boundaries" {
    try expectAtan2H(0x3400, pi.atan2pih(h(0x3c00), h(0x3c00)));
    try expectAtan2H(0xb400, pi.atan2pih(h(0xbc00), h(0x3c00)));
    try expectAtan2H(0x3a00, pi.atan2pih(h(0x3c00), h(0xbc00)));
    try expectAtan2H(0xba00, pi.atan2pih(h(0xbc00), h(0xbc00)));

    try expectAtan2F(0x3e80_0000, pi.atan2pif(f(0x3f80_0000), f(0x3f80_0000)));
    try expectAtan2F(0xbe80_0000, pi.atan2pif(f(0xbf80_0000), f(0x3f80_0000)));
    try expectAtan2F(0x3f40_0000, pi.atan2pif(f(0x3f80_0000), f(0xbf80_0000)));
    try expectAtan2F(0xbf40_0000, pi.atan2pif(f(0xbf80_0000), f(0xbf80_0000)));

    try expectAtan2D(
        0x3fd0_0000_0000_0000,
        pi.atan2pid(d(0x3ff0_0000_0000_0000), d(0x3ff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0xbfd0_0000_0000_0000,
        pi.atan2pid(d(0xbff0_0000_0000_0000), d(0x3ff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0x3fe8_0000_0000_0000,
        pi.atan2pid(d(0x3ff0_0000_0000_0000), d(0xbff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0xbfe8_0000_0000_0000,
        pi.atan2pid(d(0xbff0_0000_0000_0000), d(0xbff0_0000_0000_0000)),
    );
}

test "atan2pi finite kernels report resolved values" {
    try expectAtan2H(0x30b9, pi.atan2pih(h(0x3c00), h(0x4000)));
    try expectAtan2F(0x3e17_2029, pi.atan2pif(f(0x3f80_0000), f(0x4000_0000)));
    try expectAtan2D(
        0x3fc2_e405_1d9d_f308,
        pi.atan2pid(d(0x3ff0_0000_0000_0000), d(0x4000_0000_0000_0000)),
    );
}

test "atan2pi binary32 selects the quotient operands for unequal magnitudes" {
    // The kernel divides the larger magnitude by the smaller one, so these
    // binary32 pairs cover both comparisons with unequal magnitudes: adjacent
    // encodings, an exact ratio of two, and a wide ratio in every quadrant.
    try expectAtan2F(0x3e80_0001, pi.atan2pif(f(0x3f80_0001), f(0x3f80_0000)));
    try expectAtan2F(0x3e80_0000, pi.atan2pif(f(0x3f80_0000), f(0x3f7f_ffff)));
    try expectAtan2F(0x3eb4_6fec, pi.atan2pif(f(0x4000_0000), f(0x3f80_0000)));
    try expectAtan2F(0x3ea0_2ba4, pi.atan2pif(f(0x3fc0_0000), f(0x3f80_0000)));
    try expectAtan2F(0x3f25_c80a, pi.atan2pif(f(0x4000_0000), f(0xbf80_0000)));
    try expectAtan2F(0xbeb4_6fec, pi.atan2pif(f(0xc000_0000), f(0x3f80_0000)));
    try expectAtan2F(0xbf25_c80a, pi.atan2pif(f(0xc000_0000), f(0xbf80_0000)));
    try expectAtan2F(0x3eb4_6fec, pi.atan2pif(f(0x3f00_0000), f(0x3e80_0000)));
}

test "atan2pi binary32 rounds exponent-gap axis neighbors" {
    try expectAtan2F(0x3f00_0000, pi.atan2pif(f(0x3f80_0000), f(0x3100_0000)));
    try expectAtan2F(0x3f00_0000, pi.atan2pif(f(0x3f80_0000), f(0x317f_ffff)));
    try expectAtan2F(0x3f00_0000, pi.atan2pif(f(0x3f80_0000), f(0x3180_0000)));
    try expectAtan2F(0x3f00_0000, pi.atan2pif(f(0x3f80_0000), f(0x0000_0001)));
    try expectAtan2F(0xbf00_0000, pi.atan2pif(f(0xbf80_0000), f(0x3100_0000)));
    try expectAtan2F(0x3f80_0000, pi.atan2pif(f(0x3100_0000), f(0xbf80_0000)));
    try expectAtan2F(0xbf80_0000, pi.atan2pif(f(0xb100_0000), f(0xbf80_0000)));
    try expectAtan2F(0x3f00_0000, pi.atan2pif(f(0x3f80_0000), f(0x3300_0000)));
    try expectAtan2F(0x3f00_0000, pi.atan2pif(f(0x3f80_0000), f(0x3200_0000)));
    try expectAtan2F(0x3f80_0000, pi.atan2pif(f(0x3200_0000), f(0xbf80_0000)));
    try expectAtan2F(0x3eff_ffff, pi.atan2pif(f(0x3f80_0000), f(0x3380_0000)));
    try expectAtan2F(0x3f80_0000, pi.atan2pif(f(0x3380_0000), f(0xbf80_0000)));
    try expectAtan2F(0x3f7f_ffff, pi.atan2pif(f(0x3400_0000), f(0xbf80_0000)));
}

test "atan2pi wide gaps use certified principal bounds" {
    try expectAtan2D(
        0x3fe0_0000_0000_0000,
        pi.atan2pid(d(0x3ff0_0000_0000_0000), d(0x3c80_0000_0000_0000)),
    );
    try expectAtan2D(
        0x3fe0_0000_0000_0000,
        pi.atan2pid(d(0x3ff0_0000_0000_0000), d(0xbc80_0000_0000_0000)),
    );
    try expectAtan2D(
        0x3ff0_0000_0000_0000,
        pi.atan2pid(d(0x3c80_0000_0000_0000), d(0xbff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0x39c4_5f30_6dc9_c883,
        pi.atan2pid(d(0x39e0_0000_0000_0000), d(0x3ff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0xb9c4_5f30_6dc9_c883,
        pi.atan2pid(d(0xb9e0_0000_0000_0000), d(0x3ff0_0000_0000_0000)),
    );
    try expectAtan2D(
        0x0000_0000_0000_0001,
        pi.atan2pid(d(0x0000_0000_0000_0002), d(0x3ff0_0000_0000_0000)),
    );
}

test "atan2pi normalizes subnormal operands before fast interval" {
    try expectAtan2D(
        0x3ff0_0000_0000_0000,
        pi.atan2pid(d(0x0000_0000_0000_0001), d(0x8020_0000_0000_0000)),
    );
    try expectAtan2D(
        0xbff0_0000_0000_0000,
        pi.atan2pid(d(0x8000_0000_0000_0001), d(0x8020_0000_0000_0000)),
    );
}

test "atan2pi binary64 forced ambiguity reports unresolved" {
    const positive = pi.testing.forcedAmbiguity(1.0, 1.0);
    const negative = pi.testing.forcedAmbiguity(-1.0, -1.0);
    const positive_bits: u64 = @bitCast(positive.value);
    const negative_bits: u64 = @bitCast(negative.value);

    try std.testing.expect(!positive.resolved);
    try std.testing.expect(!negative.resolved);
    try std.testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), positive_bits);
    try std.testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), negative_bits);
}
