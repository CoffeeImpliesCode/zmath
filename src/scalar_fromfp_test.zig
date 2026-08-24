const std = @import("std");
const common = @import("scalar_c23_common.zig");
const fromfp = @import("scalar_fromfp.zig");

fn expectHBits(expected: u16, actual: f16) !void {
    try std.testing.expectEqual(expected, common.H.toBits(actual));
}

fn expectFBits(expected: u32, actual: f32) !void {
    try std.testing.expectEqual(expected, common.F.toBits(actual));
}

fn expectDBits(expected: u64, actual: f64) !void {
    try std.testing.expectEqual(expected, common.D.toBits(actual));
}

test "fromfph rounds both signs in every direction" {
    const positive = common.H.fromBits(0x4100);
    const negative = common.H.fromBits(0xc100);

    try expectHBits(0x4000, fromfp.fromfph(positive, .nearest_even, 16));
    try expectHBits(0xc000, fromfp.fromfph(negative, .nearest_even, 16));
    try expectHBits(0x4200, fromfp.fromfph(positive, .nearest_from_zero, 16));
    try expectHBits(0xc200, fromfp.fromfph(negative, .nearest_from_zero, 16));
    try expectHBits(0x4000, fromfp.fromfph(positive, .toward_zero, 16));
    try expectHBits(0xc000, fromfp.fromfph(negative, .toward_zero, 16));
    try expectHBits(0x4200, fromfp.fromfph(positive, .toward_positive, 16));
    try expectHBits(0xc000, fromfp.fromfph(negative, .toward_positive, 16));
    try expectHBits(0x4000, fromfp.fromfph(positive, .toward_negative, 16));
    try expectHBits(0xc200, fromfp.fromfph(negative, .toward_negative, 16));
}

test "fromfpf rounds both signs in every direction" {
    const positive = common.F.fromBits(0x4020_0000);
    const negative = common.F.fromBits(0xc020_0000);

    try expectFBits(0x4000_0000, fromfp.fromfpf(positive, .nearest_even, 32));
    try expectFBits(0xc000_0000, fromfp.fromfpf(negative, .nearest_even, 32));
    try expectFBits(0x4040_0000, fromfp.fromfpf(positive, .nearest_from_zero, 32));
    try expectFBits(0xc040_0000, fromfp.fromfpf(negative, .nearest_from_zero, 32));
    try expectFBits(0x4000_0000, fromfp.fromfpf(positive, .toward_zero, 32));
    try expectFBits(0xc000_0000, fromfp.fromfpf(negative, .toward_zero, 32));
    try expectFBits(0x4040_0000, fromfp.fromfpf(positive, .toward_positive, 32));
    try expectFBits(0xc000_0000, fromfp.fromfpf(negative, .toward_positive, 32));
    try expectFBits(0x4000_0000, fromfp.fromfpf(positive, .toward_negative, 32));
    try expectFBits(0xc040_0000, fromfp.fromfpf(negative, .toward_negative, 32));
}

test "fromfpd rounds both signs in every direction" {
    const positive = common.D.fromBits(0x4004_0000_0000_0000);
    const negative = common.D.fromBits(0xc004_0000_0000_0000);

    try expectDBits(
        0x4000_0000_0000_0000,
        fromfp.fromfpd(positive, .nearest_even, 64),
    );
    try expectDBits(
        0xc000_0000_0000_0000,
        fromfp.fromfpd(negative, .nearest_even, 64),
    );
    try expectDBits(
        0x4008_0000_0000_0000,
        fromfp.fromfpd(positive, .nearest_from_zero, 64),
    );
    try expectDBits(
        0xc008_0000_0000_0000,
        fromfp.fromfpd(negative, .nearest_from_zero, 64),
    );
    try expectDBits(
        0x4000_0000_0000_0000,
        fromfp.fromfpd(positive, .toward_zero, 64),
    );
    try expectDBits(
        0xc000_0000_0000_0000,
        fromfp.fromfpd(negative, .toward_zero, 64),
    );
    try expectDBits(
        0x4008_0000_0000_0000,
        fromfp.fromfpd(positive, .toward_positive, 64),
    );
    try expectDBits(
        0xc000_0000_0000_0000,
        fromfp.fromfpd(negative, .toward_positive, 64),
    );
    try expectDBits(
        0x4000_0000_0000_0000,
        fromfp.fromfpd(positive, .toward_negative, 64),
    );
    try expectDBits(
        0xc008_0000_0000_0000,
        fromfp.fromfpd(negative, .toward_negative, 64),
    );
}

test "fromfp preserves signed zero results" {
    try expectHBits(
        0x8000,
        fromfp.fromfph(common.H.fromBits(0xb800), .nearest_even, 4),
    );
    try expectFBits(
        0x8000_0000,
        fromfp.fromfpf(common.F.fromBits(0xbe80_0000), .toward_zero, 4),
    );
    try expectDBits(
        0x8000_0000_0000_0000,
        fromfp.fromfpd(
            common.D.fromBits(0xbfd0_0000_0000_0000),
            .toward_positive,
            4,
        ),
    );
}

test "fromfp signed width endpoints are exact in every format" {
    try expectHBits(0xe800, fromfp.fromfph(common.H.fromBits(0xe800), .toward_zero, 12));
    try expectHBits(0x67ff, fromfp.fromfph(common.H.fromBits(0x67ff), .toward_zero, 12));
    try expectHBits(0x7e00, fromfp.fromfph(common.H.fromBits(0x6800), .toward_zero, 12));
    try expectHBits(0x7e00, fromfp.fromfph(common.H.fromBits(0xe801), .toward_zero, 12));

    try expectFBits(
        0xcb80_0000,
        fromfp.fromfpf(common.F.fromBits(0xcb80_0000), .toward_zero, 25),
    );
    try expectFBits(
        0x4b7f_ffff,
        fromfp.fromfpf(common.F.fromBits(0x4b7f_ffff), .toward_zero, 25),
    );
    try expectFBits(
        0x7fc0_0000,
        fromfp.fromfpf(common.F.fromBits(0x4b80_0000), .toward_zero, 25),
    );
    try expectFBits(
        0x7fc0_0000,
        fromfp.fromfpf(common.F.fromBits(0xcb80_0001), .toward_zero, 25),
    );

    try expectDBits(
        0xc340_0000_0000_0000,
        fromfp.fromfpd(common.D.fromBits(0xc340_0000_0000_0000), .toward_zero, 54),
    );
    try expectDBits(
        0x433f_ffff_ffff_ffff,
        fromfp.fromfpd(common.D.fromBits(0x433f_ffff_ffff_ffff), .toward_zero, 54),
    );
    try expectDBits(
        0x7ff8_0000_0000_0000,
        fromfp.fromfpd(common.D.fromBits(0x4340_0000_0000_0000), .toward_zero, 54),
    );
    try expectDBits(
        0x7ff8_0000_0000_0000,
        fromfp.fromfpd(common.D.fromBits(0xc340_0000_0000_0001), .toward_zero, 54),
    );

    try expectHBits(0xbc00, fromfp.fromfph(common.H.fromBits(0xbc00), .toward_zero, 1));
    try expectHBits(0x0000, fromfp.fromfph(common.H.fromBits(0x0000), .toward_zero, 1));
    try expectHBits(0x7e00, fromfp.fromfph(common.H.fromBits(0x3c00), .toward_zero, 1));
}

test "ufromfp enforces unsigned width endpoints and rejects negative results" {
    try expectHBits(0x67ff, fromfp.ufromfph(common.H.fromBits(0x67ff), .toward_zero, 11));
    try expectHBits(0x7e00, fromfp.ufromfph(common.H.fromBits(0x6800), .toward_zero, 11));
    try expectHBits(0x7e00, fromfp.ufromfph(common.H.fromBits(0xbc00), .toward_zero, 11));

    try expectFBits(
        0x4b7f_ffff,
        fromfp.ufromfpf(common.F.fromBits(0x4b7f_ffff), .toward_zero, 24),
    );
    try expectFBits(
        0x7fc0_0000,
        fromfp.ufromfpf(common.F.fromBits(0x4b80_0000), .toward_zero, 24),
    );
    try expectFBits(
        0x7fc0_0000,
        fromfp.ufromfpf(common.F.fromBits(0xbf80_0000), .toward_zero, 24),
    );

    try expectDBits(
        0x433f_ffff_ffff_ffff,
        fromfp.ufromfpd(common.D.fromBits(0x433f_ffff_ffff_ffff), .toward_zero, 53),
    );
    try expectDBits(
        0x7ff8_0000_0000_0000,
        fromfp.ufromfpd(common.D.fromBits(0x4340_0000_0000_0000), .toward_zero, 53),
    );
    try expectDBits(
        0x7ff8_0000_0000_0000,
        fromfp.ufromfpd(common.D.fromBits(0xbff0_0000_0000_0000), .toward_zero, 53),
    );

    try expectHBits(0x0000, fromfp.ufromfph(common.H.fromBits(0x0000), .toward_zero, 1));
    try expectHBits(0x3c00, fromfp.ufromfph(common.H.fromBits(0x3c00), .toward_zero, 1));
    try expectHBits(0x7e00, fromfp.ufromfph(common.H.fromBits(0x4000), .toward_zero, 1));
}

test "fromfp validates the rounded value against width" {
    const positive_tie = common.H.fromBits(0x4780);
    const negative_tie = common.H.fromBits(0xc840);

    try expectHBits(0x4700, fromfp.fromfph(positive_tie, .toward_zero, 4));
    try expectHBits(0x7e00, fromfp.fromfph(positive_tie, .nearest_even, 4));
    try expectHBits(0xc800, fromfp.fromfph(negative_tie, .nearest_even, 4));
    try expectHBits(0x7e00, fromfp.fromfph(negative_tie, .nearest_from_zero, 4));
}

test "zero width and nonfinite inputs return positive canonical NaNs" {
    try expectHBits(0x7e00, fromfp.fromfph(common.H.fromBits(0x3c00), .toward_zero, 0));
    try expectFBits(
        0x7fc0_0000,
        fromfp.fromfpf(common.F.fromBits(0x3f80_0000), .toward_zero, 0),
    );
    try expectDBits(
        0x7ff8_0000_0000_0000,
        fromfp.fromfpd(common.D.fromBits(0x3ff0_0000_0000_0000), .toward_zero, 0),
    );

    try expectHBits(0x7e00, fromfp.fromfph(common.H.fromBits(0xfe55), .toward_zero, 16));
    try expectHBits(0x7e00, fromfp.fromfph(common.H.fromBits(0x7c00), .toward_zero, 16));
    try expectHBits(0x7e00, fromfp.fromfph(common.H.fromBits(0xfc00), .toward_zero, 16));

    try expectFBits(
        0x7fc0_0000,
        fromfp.fromfpf(common.F.fromBits(0xffc0_0123), .toward_zero, 32),
    );
    try expectFBits(
        0x7fc0_0000,
        fromfp.fromfpf(common.F.fromBits(0x7f80_0000), .toward_zero, 32),
    );
    try expectFBits(
        0x7fc0_0000,
        fromfp.fromfpf(common.F.fromBits(0xff80_0000), .toward_zero, 32),
    );

    try expectDBits(
        0x7ff8_0000_0000_0000,
        fromfp.fromfpd(common.D.fromBits(0xfff8_0000_0000_0123), .toward_zero, 64),
    );
    try expectDBits(
        0x7ff8_0000_0000_0000,
        fromfp.fromfpd(common.D.fromBits(0x7ff0_0000_0000_0000), .toward_zero, 64),
    );
    try expectDBits(
        0x7ff8_0000_0000_0000,
        fromfp.fromfpd(common.D.fromBits(0xfff0_0000_0000_0000), .toward_zero, 64),
    );
}

test "fromfpx and ufromfpx have the same exact value semantics" {
    try expectHBits(
        0xc200,
        fromfp.fromfpxh(common.H.fromBits(0xc100), .nearest_from_zero, 16),
    );
    try expectFBits(
        0x4040_0000,
        fromfp.fromfpxf(common.F.fromBits(0x4020_0000), .toward_positive, 32),
    );
    try expectDBits(
        0xc008_0000_0000_0000,
        fromfp.fromfpxd(common.D.fromBits(0xc004_0000_0000_0000), .toward_negative, 64),
    );

    try expectHBits(
        0x4200,
        fromfp.ufromfpxh(common.H.fromBits(0x4100), .nearest_from_zero, 16),
    );
    try expectFBits(
        0x4000_0000,
        fromfp.ufromfpxf(common.F.fromBits(0x4020_0000), .toward_zero, 32),
    );
    try expectDBits(
        0x4008_0000_0000_0000,
        fromfp.ufromfpxd(common.D.fromBits(0x4004_0000_0000_0000), .toward_positive, 64),
    );

    try expectHBits(
        0x7e00,
        fromfp.fromfpxh(common.H.fromBits(0x3c00), .nearest_even, 0),
    );
    try expectFBits(
        0x7fc0_0000,
        fromfp.ufromfpxf(common.F.fromBits(0xbf80_0000), .nearest_even, 32),
    );
    try expectDBits(
        0x7ff8_0000_0000_0000,
        fromfp.fromfpxd(common.D.fromBits(0x7ff0_0000_0000_0000), .nearest_even, 64),
    );
}

test "large widths accept every finite floating-point integer" {
    try expectHBits(
        0xfbff,
        fromfp.fromfph(common.H.fromBits(0xfbff), .toward_zero, std.math.maxInt(u32)),
    );
    try expectFBits(
        0x7f7f_ffff,
        fromfp.ufromfpf(
            common.F.fromBits(0x7f7f_ffff),
            .toward_zero,
            std.math.maxInt(u32),
        ),
    );
    try expectDBits(
        0xffef_ffff_ffff_ffff,
        fromfp.fromfpxd(
            common.D.fromBits(0xffef_ffff_ffff_ffff),
            .toward_zero,
            std.math.maxInt(u32),
        ),
    );
}
