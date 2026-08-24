const std = @import("std");
const common = @import("scalar_c23_common.zig");
const core = @import("ieee_c23_core.zig");

fn expectHBits(expected: u16, value: f16) !void {
    try std.testing.expectEqual(expected, common.H.toBits(value));
}

fn expectFBits(expected: u32, value: f32) !void {
    try std.testing.expectEqual(expected, common.F.toBits(value));
}

fn expectDBits(expected: u64, value: f64) !void {
    try std.testing.expectEqual(expected, common.D.toBits(value));
}

test "canonicalizeh preserves values and quiets both NaN classes" {
    const finite = core.canonicalizeh(common.H.fromBits(0x3c00));
    try std.testing.expect(finite.ok);
    try expectHBits(0x3c00, finite.value);

    const quiet = core.canonicalizeh(common.H.fromBits(0xfe55));
    try std.testing.expect(quiet.ok);
    try expectHBits(0xfe55, quiet.value);

    const signaling = core.canonicalizeh(common.H.fromBits(0xfd55));
    try std.testing.expect(signaling.ok);
    try expectHBits(0xff55, signaling.value);
}

test "canonicalizef preserves values and quiets both NaN classes" {
    const finite = core.canonicalizef(common.F.fromBits(0x3f80_0000));
    try std.testing.expect(finite.ok);
    try expectFBits(0x3f80_0000, finite.value);

    const quiet = core.canonicalizef(common.F.fromBits(0xffc0_0123));
    try std.testing.expect(quiet.ok);
    try expectFBits(0xffc0_0123, quiet.value);

    const signaling = core.canonicalizef(common.F.fromBits(0xff80_0123));
    try std.testing.expect(signaling.ok);
    try expectFBits(0xffc0_0123, signaling.value);
}

test "canonicalized preserves values and quiets both NaN classes" {
    const finite = core.canonicalized(common.D.fromBits(0x3ff0_0000_0000_0000));
    try std.testing.expect(finite.ok);
    try expectDBits(0x3ff0_0000_0000_0000, finite.value);

    const quiet = core.canonicalized(common.D.fromBits(0xfff8_0000_0000_0123));
    try std.testing.expect(quiet.ok);
    try expectDBits(0xfff8_0000_0000_0123, quiet.value);

    const signaling = core.canonicalized(common.D.fromBits(0xfff0_0000_0000_0123));
    try std.testing.expect(signaling.ok);
    try expectDBits(0xfff8_0000_0000_0123, signaling.value);
}

test "nanh parses strict hexadecimal diagnostic payloads" {
    try expectHBits(0x7e00, core.nanh("0"));
    try expectHBits(0x7e01, core.nanh("0001"));
    try expectHBits(0x7eaf, core.nanh("aF"));
    try expectHBits(0x7fff, core.nanh("1ff"));

    try expectHBits(0x7e00, core.nanh(""));
    try expectHBits(0x7e00, core.nanh("-1"));
    try expectHBits(0x7e00, core.nanh("0x1"));
    try expectHBits(0x7e00, core.nanh(" 1"));
    try expectHBits(0x7e00, core.nanh("200"));
    try expectHBits(0x7e00, core.nanh("\x80"));
}

test "nanf parses strict hexadecimal diagnostic payloads" {
    try expectFBits(0x7fc0_0000, core.nanf("0"));
    try expectFBits(0x7fc0_0001, core.nanf("0001"));
    try expectFBits(0x7fca_bcde, core.nanf("aBcDe"));
    try expectFBits(0x7fff_ffff, core.nanf("3fffff"));

    try expectFBits(0x7fc0_0000, core.nanf(""));
    try expectFBits(0x7fc0_0000, core.nanf("-1"));
    try expectFBits(0x7fc0_0000, core.nanf("0x1"));
    try expectFBits(0x7fc0_0000, core.nanf("1_0"));
    try expectFBits(0x7fc0_0000, core.nanf("400000"));
    try expectFBits(0x7fc0_0000, core.nanf("\x80"));
}

test "nand parses strict hexadecimal diagnostic payloads" {
    try expectDBits(0x7ff8_0000_0000_0000, core.nand("0"));
    try expectDBits(0x7ff8_0000_0000_0001, core.nand("0001"));
    try expectDBits(0x7ff8_0001_23ab_cdef, core.nand("123AbCdEf"));
    try expectDBits(0x7fff_ffff_ffff_ffff, core.nand("7ffffffffffff"));

    try expectDBits(0x7ff8_0000_0000_0000, core.nand(""));
    try expectDBits(0x7ff8_0000_0000_0000, core.nand("-1"));
    try expectDBits(0x7ff8_0000_0000_0000, core.nand("0x1"));
    try expectDBits(0x7ff8_0000_0000_0000, core.nand("1g"));
    try expectDBits(0x7ff8_0000_0000_0000, core.nand("8000000000000"));
    try expectDBits(0x7ff8_0000_0000_0000, core.nand("\x80"));
}

test "nextTowardh compares through f64 and preserves zero target signs" {
    try expectHBits(0x3c00, core.nextTowardh(common.H.fromBits(0x3c00), 1.0));
    try expectHBits(0x3c01, core.nextTowardh(common.H.fromBits(0x3c00), 1.1));
    try expectHBits(0x3bff, core.nextTowardh(common.H.fromBits(0x3c00), 0.9));
    try expectHBits(0x8000, core.nextTowardh(common.H.fromBits(0x0000), -0.0));
    try expectHBits(0x0000, core.nextTowardh(common.H.fromBits(0x8000), 0.0));
}

test "nextTowardf compares through f64 and preserves zero target signs" {
    try expectFBits(0x3f80_0000, core.nextTowardf(common.F.fromBits(0x3f80_0000), 1.0));
    try expectFBits(0x3f80_0001, core.nextTowardf(common.F.fromBits(0x3f80_0000), 1.1));
    try expectFBits(0x3f7f_ffff, core.nextTowardf(common.F.fromBits(0x3f80_0000), 0.9));
    try expectFBits(0x8000_0000, core.nextTowardf(common.F.fromBits(0x0000_0000), -0.0));
    try expectFBits(0x0000_0000, core.nextTowardf(common.F.fromBits(0x8000_0000), 0.0));
}

test "nextTowardd compares equality and adjacent f64 values" {
    try expectDBits(
        0x3ff0_0000_0000_0000,
        core.nextTowardd(common.D.fromBits(0x3ff0_0000_0000_0000), 1.0),
    );
    try expectDBits(
        0x3ff0_0000_0000_0001,
        core.nextTowardd(common.D.fromBits(0x3ff0_0000_0000_0000), 1.1),
    );
    try expectDBits(
        0x3fef_ffff_ffff_ffff,
        core.nextTowardd(common.D.fromBits(0x3ff0_0000_0000_0000), 0.9),
    );
}

test "nextToward steps away from both infinities" {
    try expectHBits(
        0x7bff,
        core.nextTowardh(common.H.fromBits(0x7c00), -std.math.inf(f64)),
    );
    try expectHBits(
        0xfbff,
        core.nextTowardh(common.H.fromBits(0xfc00), std.math.inf(f64)),
    );

    try expectFBits(
        0x7f7f_ffff,
        core.nextTowardf(common.F.fromBits(0x7f80_0000), -std.math.inf(f64)),
    );
    try expectFBits(
        0xff7f_ffff,
        core.nextTowardf(common.F.fromBits(0xff80_0000), std.math.inf(f64)),
    );

    try expectDBits(
        0x7fef_ffff_ffff_ffff,
        core.nextTowardd(common.D.fromBits(0x7ff0_0000_0000_0000), -std.math.inf(f64)),
    );
    try expectDBits(
        0xffef_ffff_ffff_ffff,
        core.nextTowardd(common.D.fromBits(0xfff0_0000_0000_0000), std.math.inf(f64)),
    );
}

test "nextToward keeps native first NaN in parameter order" {
    try expectHBits(
        0xff55,
        core.nextTowardh(
            common.H.fromBits(0xfd55),
            common.D.fromBits(0x7ffc_a800_0000_0055),
        ),
    );
    try expectHBits(
        0x7e2a,
        core.nextTowardh(
            common.H.fromBits(0x7e2a),
            common.D.fromBits(0xfff3_5400_0000_012a),
        ),
    );

    try expectFBits(
        0xffc0_0123,
        core.nextTowardf(
            common.F.fromBits(0xff80_0123),
            common.D.fromBits(0x7ff5_579b_c015_4321),
        ),
    );
    try expectFBits(
        0x7fc0_0456,
        core.nextTowardf(
            common.F.fromBits(0x7fc0_0456),
            common.D.fromBits(0xfffa_b4b4_a02a_bcde),
        ),
    );

    try expectDBits(
        0xfff8_0000_0000_0123,
        core.nextTowardd(
            common.D.fromBits(0xfff0_0000_0000_0123),
            common.D.fromBits(0x7ff8_0000_0000_0456),
        ),
    );
    try expectDBits(
        0x7ff8_0000_0000_0456,
        core.nextTowardd(
            common.D.fromBits(0x7ff8_0000_0000_0456),
            common.D.fromBits(0xfff0_0000_0000_0123),
        ),
    );
}

test "nextTowardh keeps high f64 NaN payload bits" {
    // The low slices differ so these assertions reject low-order truncation.
    try expectHBits(
        0x7f2a,
        core.nextTowardh(
            common.H.fromBits(0x3c00),
            common.D.fromBits(0x7ffc_a800_0000_0055),
        ),
    );
    try expectHBits(
        0xfed5,
        core.nextTowardh(
            common.H.fromBits(0x3c00),
            common.D.fromBits(0xfff3_5400_0000_012a),
        ),
    );
}

test "nextTowardf keeps high f64 NaN payload bits" {
    // The low slices differ so these assertions reject low-order truncation.
    try expectFBits(
        0x7fea_bcde,
        core.nextTowardf(
            common.F.fromBits(0x3f80_0000),
            common.D.fromBits(0x7ff5_579b_c015_4321),
        ),
    );
    try expectFBits(
        0xffd5_a5a5,
        core.nextTowardf(
            common.F.fromBits(0x3f80_0000),
            common.D.fromBits(0xfffa_b4b4_a02a_bcde),
        ),
    );
}

test "nextTowardd preserves target NaN bits and quiets signaling inputs" {
    try expectDBits(
        0x7ff8_0000_0000_0123,
        core.nextTowardd(
            common.D.fromBits(0x3ff0_0000_0000_0000),
            common.D.fromBits(0x7ff8_0000_0000_0123),
        ),
    );
    try expectDBits(
        0xfff8_0000_0000_0456,
        core.nextTowardd(
            common.D.fromBits(0x3ff0_0000_0000_0000),
            common.D.fromBits(0xfff0_0000_0000_0456),
        ),
    );
}

test "nextUp and nextDown cover every half-precision boundary" {
    try expectHBits(0x0001, core.nextUph(common.H.fromBits(0x0000)));
    try expectHBits(0x0001, core.nextUph(common.H.fromBits(0x8000)));
    try expectHBits(0x8001, core.nextDownh(common.H.fromBits(0x0000)));
    try expectHBits(0x8001, core.nextDownh(common.H.fromBits(0x8000)));

    try expectHBits(0x7c00, core.nextUph(common.H.fromBits(0x7c00)));
    try expectHBits(0x7bff, core.nextDownh(common.H.fromBits(0x7c00)));
    try expectHBits(0xfbff, core.nextUph(common.H.fromBits(0xfc00)));
    try expectHBits(0xfc00, core.nextDownh(common.H.fromBits(0xfc00)));
    try expectHBits(0x7c00, core.nextUph(common.H.fromBits(0x7bff)));
    try expectHBits(0xfc00, core.nextDownh(common.H.fromBits(0xfbff)));

    try expectHBits(0x7e55, core.nextUph(common.H.fromBits(0x7e55)));
    try expectHBits(0x7e55, core.nextDownh(common.H.fromBits(0x7e55)));
    try expectHBits(0xfe55, core.nextUph(common.H.fromBits(0xfe55)));
    try expectHBits(0xfe55, core.nextDownh(common.H.fromBits(0xfe55)));
    try expectHBits(0x7f55, core.nextUph(common.H.fromBits(0x7d55)));
    try expectHBits(0x7f55, core.nextDownh(common.H.fromBits(0x7d55)));
    try expectHBits(0xff55, core.nextUph(common.H.fromBits(0xfd55)));
    try expectHBits(0xff55, core.nextDownh(common.H.fromBits(0xfd55)));
}

test "nextUp and nextDown cover every single-precision boundary" {
    try expectFBits(0x0000_0001, core.nextUpf(common.F.fromBits(0x0000_0000)));
    try expectFBits(0x0000_0001, core.nextUpf(common.F.fromBits(0x8000_0000)));
    try expectFBits(0x8000_0001, core.nextDownf(common.F.fromBits(0x0000_0000)));
    try expectFBits(0x8000_0001, core.nextDownf(common.F.fromBits(0x8000_0000)));

    try expectFBits(0x7f80_0000, core.nextUpf(common.F.fromBits(0x7f80_0000)));
    try expectFBits(0x7f7f_ffff, core.nextDownf(common.F.fromBits(0x7f80_0000)));
    try expectFBits(0xff7f_ffff, core.nextUpf(common.F.fromBits(0xff80_0000)));
    try expectFBits(0xff80_0000, core.nextDownf(common.F.fromBits(0xff80_0000)));
    try expectFBits(0x7f80_0000, core.nextUpf(common.F.fromBits(0x7f7f_ffff)));
    try expectFBits(0xff80_0000, core.nextDownf(common.F.fromBits(0xff7f_ffff)));

    try expectFBits(0x7fc0_0123, core.nextUpf(common.F.fromBits(0x7fc0_0123)));
    try expectFBits(0x7fc0_0123, core.nextDownf(common.F.fromBits(0x7fc0_0123)));
    try expectFBits(0xffc0_0123, core.nextUpf(common.F.fromBits(0xffc0_0123)));
    try expectFBits(0xffc0_0123, core.nextDownf(common.F.fromBits(0xffc0_0123)));
    try expectFBits(0x7fc0_0123, core.nextUpf(common.F.fromBits(0x7f80_0123)));
    try expectFBits(0x7fc0_0123, core.nextDownf(common.F.fromBits(0x7f80_0123)));
    try expectFBits(0xffc0_0123, core.nextUpf(common.F.fromBits(0xff80_0123)));
    try expectFBits(0xffc0_0123, core.nextDownf(common.F.fromBits(0xff80_0123)));
}

test "nextUp and nextDown cover every double-precision boundary" {
    try expectDBits(0x0000_0000_0000_0001, core.nextUpd(common.D.fromBits(0x0000_0000_0000_0000)));
    try expectDBits(0x0000_0000_0000_0001, core.nextUpd(common.D.fromBits(0x8000_0000_0000_0000)));
    try expectDBits(0x8000_0000_0000_0001, core.nextDownd(common.D.fromBits(0x0000_0000_0000_0000)));
    try expectDBits(0x8000_0000_0000_0001, core.nextDownd(common.D.fromBits(0x8000_0000_0000_0000)));

    try expectDBits(0x7ff0_0000_0000_0000, core.nextUpd(common.D.fromBits(0x7ff0_0000_0000_0000)));
    try expectDBits(0x7fef_ffff_ffff_ffff, core.nextDownd(common.D.fromBits(0x7ff0_0000_0000_0000)));
    try expectDBits(0xffef_ffff_ffff_ffff, core.nextUpd(common.D.fromBits(0xfff0_0000_0000_0000)));
    try expectDBits(0xfff0_0000_0000_0000, core.nextDownd(common.D.fromBits(0xfff0_0000_0000_0000)));
    try expectDBits(0x7ff0_0000_0000_0000, core.nextUpd(common.D.fromBits(0x7fef_ffff_ffff_ffff)));
    try expectDBits(0xfff0_0000_0000_0000, core.nextDownd(common.D.fromBits(0xffef_ffff_ffff_ffff)));

    try expectDBits(0x7ff8_0000_0000_0123, core.nextUpd(common.D.fromBits(0x7ff8_0000_0000_0123)));
    try expectDBits(0x7ff8_0000_0000_0123, core.nextDownd(common.D.fromBits(0x7ff8_0000_0000_0123)));
    try expectDBits(0xfff8_0000_0000_0123, core.nextUpd(common.D.fromBits(0xfff8_0000_0000_0123)));
    try expectDBits(0xfff8_0000_0000_0123, core.nextDownd(common.D.fromBits(0xfff8_0000_0000_0123)));
    try expectDBits(0x7ff8_0000_0000_0123, core.nextUpd(common.D.fromBits(0x7ff0_0000_0000_0123)));
    try expectDBits(0x7ff8_0000_0000_0123, core.nextDownd(common.D.fromBits(0x7ff0_0000_0000_0123)));
    try expectDBits(0xfff8_0000_0000_0123, core.nextUpd(common.D.fromBits(0xfff0_0000_0000_0123)));
    try expectDBits(0xfff8_0000_0000_0123, core.nextDownd(common.D.fromBits(0xfff0_0000_0000_0123)));
}
