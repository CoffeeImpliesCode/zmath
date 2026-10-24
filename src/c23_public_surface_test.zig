const std = @import("std");

// Root passes its namespace to this module so the public-surface test has no import cycle.
const unary_float_families = [_][]const u8{
    "abs",
    "ceil",
    "floor",
    "nearbyint",
    "rint",
    "round",
    "roundeven",
    "trunc",
    "exp",
    "exp2",
    "exp10",
    "expm1",
    "exp2m1",
    "exp10m1",
    "log",
    "log2",
    "log10",
    "log1p",
    "logp1",
    "log2p1",
    "log10p1",
    "sin",
    "cos",
    "tan",
    "asin",
    "acos",
    "atan",
    "sinpi",
    "cospi",
    "tanpi",
    "asinpi",
    "acospi",
    "atanpi",
    "sinh",
    "cosh",
    "tanh",
    "asinh",
    "acosh",
    "atanh",
    "sqrt",
    "cbrt",
    "rsqrt",
    "erf",
    "erfc",
    "tgamma",
    "logb",
};

const binary_float_families = [_][]const u8{
    "pow",
    "powr",
    "hypot",
    "fmod",
    "remainder",
    "fdim",
    "fmax",
    "fmin",
    "fmaximum",
    "fminimum",
    "fmaximum_mag",
    "fminimum_mag",
    "fmaximum_num",
    "fminimum_num",
    "fmaximum_mag_num",
    "fminimum_mag_num",
};

const atan2_result_families = [_][]const u8{
    "atan2",
    "atan2pi",
};

const power_result_families = [_][]const u8{
    "pown",
    "compoundn",
    "rootn",
};

const integer_result_families = [_][]const u8{
    "lrint",
    "llrint",
    "lround",
    "llround",
};

const fromfp_families = [_][]const u8{
    "fromfp",
    "ufromfp",
    "fromfpx",
    "ufromfpx",
};

const scale_i32_families = [_][]const u8{
    "ldexp",
    "scalbn",
};

const ieee_families = [_][]const u8{
    "canonicalize",
    "copySign",
    "nan",
    "nextAfter",
    "nextToward",
    "nextUp",
    "nextDown",
    "totalOrder",
    "totalOrderMag",
    "getPayload",
    "setPayload",
    "setPayloadSig",
};

pub fn run(comptime zmath: type) !void {
    comptime {
        // The exhaustive C23 reflection matrix intentionally exceeds the default quota.
        @setEvalBranchQuota(2_000_000);
        assertResultContainers(zmath.scalar);
        assertResultContainers(zmath);
        assertIeeeSurface(zmath.scalar, zmath.ieee);
        assertRootSurface(zmath);
        assertUnsupportedSurface(zmath);
    }

    try expectPublicStatuses(zmath);
    try expectRootTgammaSpecials(zmath);
}

fn declarationName(comptime base: []const u8, comptime suffix: []const u8) []const u8 {
    return std.fmt.comptimePrint("{s}{s}", .{ base, suffix });
}

fn narrowingName(
    comptime destination: []const u8,
    comptime operation: []const u8,
    comptime source: []const u8,
) []const u8 {
    return std.fmt.comptimePrint("{s}{s}{s}", .{ destination, operation, source });
}

fn requireDecl(comptime Namespace: type, comptime name: []const u8) void {
    if (!@hasDecl(Namespace, name)) {
        @compileError(std.fmt.comptimePrint("missing public C23 symbol `{s}`", .{name}));
    }
}

fn rejectDecl(comptime Namespace: type, comptime name: []const u8) void {
    if (@hasDecl(Namespace, name)) {
        @compileError(std.fmt.comptimePrint("unsupported public C23 symbol `{s}`", .{name}));
    }
}

fn requireFunctionType(
    comptime Namespace: type,
    comptime Expected: type,
    comptime name: []const u8,
) void {
    requireDecl(Namespace, name);
    const Actual = @TypeOf(@field(Namespace, name));
    if (Actual != Expected) {
        @compileError(std.fmt.comptimePrint(
            "public C23 symbol `{s}` has type `{s}`, expected `{s}`",
            .{ name, @typeName(Actual), @typeName(Expected) },
        ));
    }
}

fn requireGenericCallType(
    comptime Namespace: type,
    comptime Expected: type,
    comptime name: []const u8,
    comptime args: anytype,
) void {
    requireDecl(Namespace, name);
    const Function = @TypeOf(@field(Namespace, name));
    const function_info = @typeInfo(Function).@"fn";
    if (!function_info.is_generic) {
        @compileError(std.fmt.comptimePrint(
            "root C23 symbol `{s}` must reject unsupported types at compile time",
            .{name},
        ));
    }

    const Actual = @TypeOf(@call(.auto, @field(Namespace, name), args));
    if (Actual != Expected) {
        @compileError(std.fmt.comptimePrint(
            "root C23 call `{s}` has type `{s}`, expected `{s}`",
            .{ name, @typeName(Actual), @typeName(Expected) },
        ));
    }
}

fn assertScalarUnaryFloatFamily(comptime Scalar: type, comptime name: []const u8) void {
    requireFunctionType(Scalar, fn (f16) f16, declarationName(name, "h"));
    requireFunctionType(Scalar, fn (f32) f32, declarationName(name, "f"));
    requireFunctionType(Scalar, fn (f64) f64, declarationName(name, "d"));
}

fn assertScalarBinaryFloatFamily(comptime Scalar: type, comptime name: []const u8) void {
    requireFunctionType(Scalar, fn (f16, f16) f16, declarationName(name, "h"));
    requireFunctionType(Scalar, fn (f32, f32) f32, declarationName(name, "f"));
    requireFunctionType(Scalar, fn (f64, f64) f64, declarationName(name, "d"));
}

fn assertScalarSurface(comptime Scalar: type) void {
    inline for (unary_float_families) |name| {
        assertScalarUnaryFloatFamily(Scalar, name);
    }
    inline for (binary_float_families) |name| {
        assertScalarBinaryFloatFamily(Scalar, name);
    }

    requireFunctionType(Scalar, fn (f16, f16, f16) f16, "fmah");
    requireFunctionType(Scalar, fn (f32, f32, f32) f32, "fmaf");
    requireFunctionType(Scalar, fn (f64, f64, f64) f64, "fmad");

    inline for (atan2_result_families) |name| {
        requireFunctionType(
            Scalar,
            fn (f16, f16) Scalar.Atan2Result(f16),
            declarationName(name, "h"),
        );
        requireFunctionType(
            Scalar,
            fn (f32, f32) Scalar.Atan2Result(f32),
            declarationName(name, "f"),
        );
        requireFunctionType(
            Scalar,
            fn (f64, f64) Scalar.Atan2Result(f64),
            declarationName(name, "d"),
        );
    }

    inline for (power_result_families) |name| {
        requireFunctionType(
            Scalar,
            fn (f16, i64) Scalar.PowerResult(f16),
            declarationName(name, "h"),
        );
        requireFunctionType(
            Scalar,
            fn (f32, i64) Scalar.PowerResult(f32),
            declarationName(name, "f"),
        );
        requireFunctionType(
            Scalar,
            fn (f64, i64) Scalar.PowerResult(f64),
            declarationName(name, "d"),
        );
    }

    inline for (integer_result_families) |name| {
        requireFunctionType(
            Scalar,
            fn (f16) Scalar.IntegerResult,
            declarationName(name, "h"),
        );
        requireFunctionType(
            Scalar,
            fn (f32) Scalar.IntegerResult,
            declarationName(name, "f"),
        );
        requireFunctionType(
            Scalar,
            fn (f64) Scalar.IntegerResult,
            declarationName(name, "d"),
        );
    }

    inline for (fromfp_families) |name| {
        requireFunctionType(
            Scalar,
            fn (f16, Scalar.FromFpRound, u32) f16,
            declarationName(name, "h"),
        );
        requireFunctionType(
            Scalar,
            fn (f32, Scalar.FromFpRound, u32) f32,
            declarationName(name, "f"),
        );
        requireFunctionType(
            Scalar,
            fn (f64, Scalar.FromFpRound, u32) f64,
            declarationName(name, "d"),
        );
    }

    requireFunctionType(Scalar, fn (f16) Scalar.LgammaResult(f16), "lgammah");
    requireFunctionType(Scalar, fn (f32) Scalar.LgammaResult(f32), "lgammaf");
    requireFunctionType(Scalar, fn (f64) Scalar.LgammaResult(f64), "lgammad");

    requireFunctionType(Scalar, fn (f16) Scalar.FrexpResult(f16), "frexph");
    requireFunctionType(Scalar, fn (f32) Scalar.FrexpResult(f32), "frexpf");
    requireFunctionType(Scalar, fn (f64) Scalar.FrexpResult(f64), "frexpd");

    requireFunctionType(Scalar, fn (f16) Scalar.ModfResult(f16), "modfh");
    requireFunctionType(Scalar, fn (f32) Scalar.ModfResult(f32), "modff");
    requireFunctionType(Scalar, fn (f64) Scalar.ModfResult(f64), "modfd");

    requireFunctionType(Scalar, fn (f16, f16) Scalar.RemquoResult(f16), "remquoh");
    requireFunctionType(Scalar, fn (f32, f32) Scalar.RemquoResult(f32), "remquof");
    requireFunctionType(Scalar, fn (f64, f64) Scalar.RemquoResult(f64), "remquod");

    requireFunctionType(Scalar, fn (f16) i32, "ilogbh");
    requireFunctionType(Scalar, fn (f32) i32, "ilogbf");
    requireFunctionType(Scalar, fn (f64) i32, "ilogbd");
    requireFunctionType(Scalar, fn (f16) i64, "llogbh");
    requireFunctionType(Scalar, fn (f32) i64, "llogbf");
    requireFunctionType(Scalar, fn (f64) i64, "llogbd");

    inline for (scale_i32_families) |name| {
        requireFunctionType(Scalar, fn (f16, i32) f16, declarationName(name, "h"));
        requireFunctionType(Scalar, fn (f32, i32) f32, declarationName(name, "f"));
        requireFunctionType(Scalar, fn (f64, i32) f64, declarationName(name, "d"));
    }
    requireFunctionType(Scalar, fn (f16, i64) f16, "scalblnh");
    requireFunctionType(Scalar, fn (f32, i64) f32, "scalblnf");
    requireFunctionType(Scalar, fn (f64, i64) f64, "scalblnd");

    assertNarrowingPair(Scalar, f32, f16, "f16", "f32");
    assertNarrowingPair(Scalar, f64, f16, "f16", "f64");
    assertNarrowingPair(Scalar, f64, f32, "f32", "f64");
}

fn assertNarrowingPair(
    comptime Namespace: type,
    comptime Source: type,
    comptime Destination: type,
    comptime destination: []const u8,
    comptime source: []const u8,
) void {
    requireFunctionType(
        Namespace,
        fn (Source, Source) Destination,
        narrowingName(destination, "add", source),
    );
    requireFunctionType(
        Namespace,
        fn (Source, Source) Destination,
        narrowingName(destination, "sub", source),
    );
    requireFunctionType(
        Namespace,
        fn (Source, Source) Destination,
        narrowingName(destination, "mul", source),
    );
    requireFunctionType(
        Namespace,
        fn (Source, Source) Destination,
        narrowingName(destination, "div", source),
    );
    requireFunctionType(
        Namespace,
        fn (Source, Source, Source) Destination,
        narrowingName(destination, "fma", source),
    );
    requireFunctionType(
        Namespace,
        fn (Source) Destination,
        narrowingName(destination, "sqrt", source),
    );
}

fn assertIeeeSurface(comptime Scalar: type, comptime Ieee: type) void {
    requireFunctionType(
        Ieee,
        fn (f16) Scalar.CanonicalizeResult(f16),
        "canonicalizeh",
    );
    requireFunctionType(
        Ieee,
        fn (f32) Scalar.CanonicalizeResult(f32),
        "canonicalizef",
    );
    requireFunctionType(
        Ieee,
        fn (f64) Scalar.CanonicalizeResult(f64),
        "canonicalized",
    );

    assertScalarBinaryFloatFamily(Ieee, "copySign");
    requireFunctionType(Ieee, fn ([]const u8) f16, "nanh");
    requireFunctionType(Ieee, fn ([]const u8) f32, "nanf");
    requireFunctionType(Ieee, fn ([]const u8) f64, "nand");
    assertScalarBinaryFloatFamily(Ieee, "nextAfter");
    requireFunctionType(Ieee, fn (f16, f64) f16, "nextTowardh");
    requireFunctionType(Ieee, fn (f32, f64) f32, "nextTowardf");
    requireFunctionType(Ieee, fn (f64, f64) f64, "nextTowardd");
    assertScalarUnaryFloatFamily(Ieee, "nextUp");
    assertScalarUnaryFloatFamily(Ieee, "nextDown");

    requireFunctionType(Ieee, fn (f16, f16) bool, "totalOrderh");
    requireFunctionType(Ieee, fn (f32, f32) bool, "totalOrderf");
    requireFunctionType(Ieee, fn (f64, f64) bool, "totalOrderd");
    requireFunctionType(Ieee, fn (f16, f16) bool, "totalOrderMagh");
    requireFunctionType(Ieee, fn (f32, f32) bool, "totalOrderMagf");
    requireFunctionType(Ieee, fn (f64, f64) bool, "totalOrderMagd");
    assertScalarUnaryFloatFamily(Ieee, "getPayload");

    requireFunctionType(Ieee, fn (f16) Scalar.PayloadResult(f16), "setPayloadh");
    requireFunctionType(Ieee, fn (f32) Scalar.PayloadResult(f32), "setPayloadf");
    requireFunctionType(Ieee, fn (f64) Scalar.PayloadResult(f64), "setPayloadd");
    requireFunctionType(Ieee, fn (f16) Scalar.PayloadResult(f16), "setPayloadSigh");
    requireFunctionType(Ieee, fn (f32) Scalar.PayloadResult(f32), "setPayloadSigf");
    requireFunctionType(Ieee, fn (f64) Scalar.PayloadResult(f64), "setPayloadSigd");
}

fn assertRootUnaryFloatFamily(comptime zmath: type, comptime name: []const u8) void {
    inline for (.{ f16, f32, f64 }) |T| {
        requireGenericCallType(zmath, T, name, .{@as(T, 1.0)});
    }
}

fn assertRootBinaryFloatFamily(comptime zmath: type, comptime name: []const u8) void {
    inline for (.{ f16, f32, f64 }) |T| {
        requireGenericCallType(zmath, T, name, .{ @as(T, 1.0), @as(T, 1.0) });
    }
}

fn assertRootSurface(comptime zmath: type) void {
    const Scalar = zmath.scalar;

    inline for (unary_float_families) |name| {
        assertRootUnaryFloatFamily(zmath, name);
    }
    inline for (binary_float_families) |name| {
        assertRootBinaryFloatFamily(zmath, name);
    }

    inline for (.{ f16, f32, f64 }) |T| {
        requireGenericCallType(
            zmath,
            T,
            "fma",
            .{ @as(T, 1.0), @as(T, 1.0), @as(T, 1.0) },
        );

        inline for (atan2_result_families) |name| {
            requireGenericCallType(
                zmath,
                Scalar.Atan2Result(T),
                name,
                .{ @as(T, 0.0), @as(T, 1.0) },
            );
        }
        inline for (power_result_families) |name| {
            requireGenericCallType(
                zmath,
                Scalar.PowerResult(T),
                name,
                .{ @as(T, 1.0), @as(i64, 1) },
            );
        }
        inline for (integer_result_families) |name| {
            requireGenericCallType(
                zmath,
                Scalar.IntegerResult,
                name,
                .{@as(T, 1.0)},
            );
        }
        inline for (fromfp_families) |name| {
            requireGenericCallType(
                zmath,
                T,
                name,
                .{ @as(T, 1.0), @as(Scalar.FromFpRound, .nearest_even), @as(u32, 8) },
            );
        }

        requireGenericCallType(zmath, Scalar.LgammaResult(T), "lgamma", .{@as(T, 1.0)});
        requireGenericCallType(zmath, Scalar.FrexpResult(T), "frexp", .{@as(T, 1.0)});
        requireGenericCallType(zmath, Scalar.ModfResult(T), "modf", .{@as(T, 1.0)});
        requireGenericCallType(
            zmath,
            Scalar.RemquoResult(T),
            "remquo",
            .{ @as(T, 1.0), @as(T, 1.0) },
        );
        requireGenericCallType(zmath, i32, "ilogb", .{@as(T, 1.0)});
        requireGenericCallType(zmath, i64, "llogb", .{@as(T, 1.0)});

        inline for (scale_i32_families) |name| {
            requireGenericCallType(
                zmath,
                T,
                name,
                .{ @as(T, 1.0), @as(i32, 1) },
            );
        }
        requireGenericCallType(
            zmath,
            T,
            "scalbln",
            .{ @as(T, 1.0), @as(i64, 1) },
        );
    }

    assertRootIeeeSurface(zmath);
    assertRootNarrowingSurface(zmath);
}

fn assertRootIeeeSurface(comptime zmath: type) void {
    const Scalar = zmath.scalar;

    inline for (.{ f16, f32, f64 }) |T| {
        requireGenericCallType(
            zmath,
            Scalar.CanonicalizeResult(T),
            "canonicalize",
            .{@as(T, 1.0)},
        );
        requireGenericCallType(
            zmath,
            T,
            "copySign",
            .{ @as(T, 1.0), @as(T, -1.0) },
        );
        requireGenericCallType(zmath, T, "nan", .{ T, "1" });
        requireGenericCallType(
            zmath,
            T,
            "nextAfter",
            .{ @as(T, 1.0), @as(T, 2.0) },
        );
        requireGenericCallType(
            zmath,
            T,
            "nextToward",
            .{ @as(T, 1.0), @as(f64, 2.0) },
        );
        requireGenericCallType(zmath, T, "nextUp", .{@as(T, 1.0)});
        requireGenericCallType(zmath, T, "nextDown", .{@as(T, 1.0)});
        requireGenericCallType(
            zmath,
            bool,
            "totalOrder",
            .{ @as(T, 1.0), @as(T, 2.0) },
        );
        requireGenericCallType(
            zmath,
            bool,
            "totalOrderMag",
            .{ @as(T, 1.0), @as(T, 2.0) },
        );
        requireGenericCallType(zmath, T, "getPayload", .{@as(T, 1.0)});
        requireGenericCallType(
            zmath,
            Scalar.PayloadResult(T),
            "setPayload",
            .{@as(T, 1.0)},
        );
        requireGenericCallType(
            zmath,
            Scalar.PayloadResult(T),
            "setPayloadSig",
            .{@as(T, 1.0)},
        );
    }
}

fn assertRootNarrowingSurface(comptime zmath: type) void {
    requireFunctionType(zmath, fn (f64, f64) f32, "fadd");
    requireFunctionType(zmath, fn (f64, f64) f32, "fsub");
    requireFunctionType(zmath, fn (f64, f64) f32, "fmul");
    requireFunctionType(zmath, fn (f64, f64) f32, "fdiv");
    requireFunctionType(zmath, fn (f64, f64, f64) f32, "ffma");
    requireFunctionType(zmath, fn (f64) f32, "fsqrt");
}

fn assertField(
    comptime FieldType: type,
    comptime index: usize,
    comptime name: []const u8,
    comptime fields: anytype,
) void {
    if (!std.mem.eql(u8, fields.field_names[index], name)) {
        @compileError(std.fmt.comptimePrint(
            "C23 result field {d} is `{s}`, expected `{s}`",
            .{ index, fields.field_names[index], name },
        ));
    }
    if (fields.field_types[index] != FieldType) {
        @compileError(std.fmt.comptimePrint(
            "C23 result field `{s}` has type `{s}`, expected `{s}`",
            .{ name, @typeName(fields.field_types[index]), @typeName(FieldType) },
        ));
    }
}

fn requireTwoFields(comptime Result: type) std.lang.Type.Struct {
    const fields = @typeInfo(Result).@"struct";
    if (fields.field_names.len != 2) {
        @compileError(std.fmt.comptimePrint(
            "C23 result `{s}` must contain exactly two fields",
            .{@typeName(Result)},
        ));
    }
    return fields;
}

fn assertResultContainers(comptime Scalar: type) void {
    inline for (.{
        "FromFpRound",
        "IntegerResult",
        "Atan2Result",
        "PowerResult",
        "FrexpResult",
        "ModfResult",
        "RemquoResult",
        "LgammaResult",
        "CanonicalizeResult",
        "PayloadResult",
    }) |name| {
        requireDecl(Scalar, name);
    }

    const directions = @typeInfo(Scalar.FromFpRound).@"enum".field_names;
    if (directions.len != 5) @compileError("FromFpRound must contain five directions");
    if (!std.mem.eql(u8, directions[0], "nearest_even")) @compileError("bad FromFpRound");
    if (!std.mem.eql(u8, directions[1], "nearest_from_zero")) @compileError("bad FromFpRound");
    if (!std.mem.eql(u8, directions[2], "toward_zero")) @compileError("bad FromFpRound");
    if (!std.mem.eql(u8, directions[3], "toward_positive")) @compileError("bad FromFpRound");
    if (!std.mem.eql(u8, directions[4], "toward_negative")) @compileError("bad FromFpRound");

    const integer_fields = requireTwoFields(Scalar.IntegerResult);
    assertField(i64, 0, "value", integer_fields);
    assertField(bool, 1, "ok", integer_fields);

    inline for (.{ f16, f32, f64 }) |T| {
        const Atan2 = Scalar.Atan2Result(T);
        const atan2_fields = requireTwoFields(Atan2);
        assertField(T, 0, "value", atan2_fields);
        assertField(bool, 1, "resolved", atan2_fields);

        const Power = Scalar.PowerResult(T);
        const power_fields = requireTwoFields(Power);
        assertField(T, 0, "value", power_fields);
        assertField(bool, 1, "ok", power_fields);

        const Frexp = Scalar.FrexpResult(T);
        const frexp_fields = requireTwoFields(Frexp);
        assertField(T, 0, "fraction", frexp_fields);
        assertField(i32, 1, "exponent", frexp_fields);

        const Modf = Scalar.ModfResult(T);
        const modf_fields = requireTwoFields(Modf);
        assertField(T, 0, "fractional", modf_fields);
        assertField(T, 1, "integral", modf_fields);

        const Remquo = Scalar.RemquoResult(T);
        const remquo_fields = requireTwoFields(Remquo);
        assertField(T, 0, "value", remquo_fields);
        assertField(i32, 1, "quotient", remquo_fields);

        const Lgamma = Scalar.LgammaResult(T);
        const lgamma_fields = requireTwoFields(Lgamma);
        assertField(T, 0, "value", lgamma_fields);
        assertField(i2, 1, "sign", lgamma_fields);

        const Canonicalize = Scalar.CanonicalizeResult(T);
        const canonicalize_fields = requireTwoFields(Canonicalize);
        assertField(T, 0, "value", canonicalize_fields);
        assertField(bool, 1, "ok", canonicalize_fields);

        const Payload = Scalar.PayloadResult(T);
        const payload_fields = requireTwoFields(Payload);
        assertField(T, 0, "value", payload_fields);
        assertField(bool, 1, "ok", payload_fields);
    }
}

fn rejectUnsupportedSuffixes(comptime Namespace: type, comptime families: anytype) void {
    inline for (families) |name| {
        inline for (.{ "H", "F", "D", "q", "f128", "bf16", "bfloat16" }) |suffix| {
            rejectDecl(Namespace, declarationName(name, suffix));
        }
    }
}

fn rejectNarrowingPair(
    comptime Namespace: type,
    comptime destination: []const u8,
    comptime source: []const u8,
) void {
    inline for (.{ "add", "sub", "mul", "div", "fma", "sqrt" }) |operation| {
        rejectDecl(Namespace, narrowingName(destination, operation, source));
    }
}

fn assertUnsupportedSurface(comptime zmath: type) void {
    const Scalar = zmath.scalar;
    const Ieee = zmath.ieee;

    rejectUnsupportedSuffixes(Scalar, unary_float_families);
    rejectUnsupportedSuffixes(Scalar, binary_float_families);
    rejectUnsupportedSuffixes(Scalar, atan2_result_families);
    rejectUnsupportedSuffixes(Scalar, power_result_families);
    rejectUnsupportedSuffixes(Scalar, integer_result_families);
    rejectUnsupportedSuffixes(Scalar, fromfp_families);
    rejectUnsupportedSuffixes(Scalar, scale_i32_families);
    rejectUnsupportedSuffixes(Scalar, .{ "fma", "lgamma", "frexp", "modf", "remquo" });
    rejectUnsupportedSuffixes(Scalar, .{ "ilogb", "llogb", "scalbln" });
    rejectUnsupportedSuffixes(Ieee, ieee_families);

    rejectDecl(zmath, "fabs");
    rejectDecl(Scalar, "fabs");
    rejectDecl(Scalar, "fabsh");
    rejectDecl(Scalar, "fabsf");
    rejectDecl(Scalar, "fabsd");

    inline for (.{ "fadd", "fsub", "fmul", "fdiv", "ffma", "fsqrt" }) |name| {
        rejectDecl(Scalar, name);
    }

    rejectNarrowingPair(zmath, "f16", "f32");
    rejectNarrowingPair(zmath, "f16", "f64");
    rejectNarrowingPair(zmath, "f32", "f64");

    rejectNarrowingPair(Scalar, "f16", "f16");
    rejectNarrowingPair(Scalar, "f32", "f16");
    rejectNarrowingPair(Scalar, "f32", "f32");
    rejectNarrowingPair(Scalar, "f64", "f16");
    rejectNarrowingPair(Scalar, "f64", "f32");
    rejectNarrowingPair(Scalar, "f64", "f64");
}

fn FloatBits(comptime T: type) type {
    if (T == f16) return u16;
    if (T == f32) return u32;
    if (T == f64) return u64;
    @compileError("unsupported public C23 floating type");
}

fn canonicalNaNBits(comptime T: type) FloatBits(T) {
    if (T == f16) return 0x7e00;
    if (T == f32) return 0x7fc0_0000;
    if (T == f64) return 0x7ff8_0000_0000_0000;
    unreachable;
}

fn expectPublicStatuses(comptime zmath: type) !void {
    inline for (.{ f16, f32, f64 }) |T| {
        try expectAtan2Statuses(zmath, T);
        try expectIntegerStatuses(zmath, T);
        try expectPowerRangeStatuses(zmath, T);
        try expectIeeeStatuses(zmath, T);
    }
}

fn expectAtan2Statuses(comptime zmath: type, comptime T: type) !void {
    const Scalar = zmath.scalar;

    inline for (atan2_result_families) |name| {
        const result: Scalar.Atan2Result(T) = @call(
            .auto,
            @field(zmath, name),
            .{ @as(T, 0.0), @as(T, 1.0) },
        );
        try std.testing.expect(result.resolved);
        _ = result.value;
    }
}

fn expectIntegerStatuses(comptime zmath: type, comptime T: type) !void {
    const Scalar = zmath.scalar;

    inline for (integer_result_families) |name| {
        const result: Scalar.IntegerResult = @call(
            .auto,
            @field(zmath, name),
            .{@as(T, 1.0)},
        );
        try std.testing.expect(result.ok);
        try std.testing.expectEqual(@as(i64, 1), result.value);
    }
}

fn expectPowerRangeStatuses(comptime zmath: type, comptime T: type) !void {
    const Scalar = zmath.scalar;

    inline for (power_result_families) |name| {
        inline for (.{ @as(i64, -65), @as(i64, 65) }) |exponent| {
            const result: Scalar.PowerResult(T) = @call(
                .auto,
                @field(zmath, name),
                .{ @as(T, 1.0), exponent },
            );
            try std.testing.expect(!result.ok);
            const result_bits: FloatBits(T) = @bitCast(result.value);
            try std.testing.expectEqual(canonicalNaNBits(T), result_bits);
        }
    }
}

fn expectIeeeStatuses(comptime zmath: type, comptime T: type) !void {
    const Scalar = zmath.scalar;

    const canonicalized: Scalar.CanonicalizeResult(T) = zmath.canonicalize(@as(T, -0.0));
    try std.testing.expect(canonicalized.ok);
    _ = canonicalized.value;

    const payload: Scalar.PayloadResult(T) = zmath.setPayload(@as(T, 1.0));
    try std.testing.expect(payload.ok);
    _ = payload.value;

    const signaling_payload: Scalar.PayloadResult(T) = zmath.setPayloadSig(@as(T, 1.0));
    try std.testing.expect(signaling_payload.ok);
    _ = signaling_payload.value;
}

fn expectTgammaBits(
    comptime zmath: type,
    comptime T: type,
    input_bits: FloatBits(T),
    expected_bits: FloatBits(T),
) !void {
    const input: T = @bitCast(input_bits);
    const actual_bits: FloatBits(T) = @bitCast(zmath.tgamma(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectRootTgammaSpecials(comptime zmath: type) !void {
    const h_cases = [_][2]u16{
        .{ 0x7e55, 0x7e55 },
        .{ 0xfe55, 0xfe55 },
        .{ 0x7d55, 0x7f55 },
        .{ 0xfd55, 0xff55 },
        .{ 0x7c00, 0x7c00 },
        .{ 0xfc00, 0x7e00 },
        .{ 0x0000, 0x7c00 },
        .{ 0x8000, 0xfc00 },
        .{ 0xbc00, 0x7e00 },
        .{ 0xc000, 0x7e00 },
    };
    const f_cases = [_][2]u32{
        .{ 0x7fc1_2345, 0x7fc1_2345 },
        .{ 0xffc1_2345, 0xffc1_2345 },
        .{ 0x7fa1_2345, 0x7fe1_2345 },
        .{ 0xffa1_2345, 0xffe1_2345 },
        .{ 0x7f80_0000, 0x7f80_0000 },
        .{ 0xff80_0000, 0x7fc0_0000 },
        .{ 0x0000_0000, 0x7f80_0000 },
        .{ 0x8000_0000, 0xff80_0000 },
        .{ 0xbf80_0000, 0x7fc0_0000 },
        .{ 0xc000_0000, 0x7fc0_0000 },
    };
    const d_cases = [_][2]u64{
        .{ 0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc },
        .{ 0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc },
        .{ 0x7ff1_2345_6789_abcd, 0x7ff9_2345_6789_abcd },
        .{ 0xfff1_2345_6789_abcd, 0xfff9_2345_6789_abcd },
        .{ 0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000 },
        .{ 0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000 },
        .{ 0x0000_0000_0000_0000, 0x7ff0_0000_0000_0000 },
        .{ 0x8000_0000_0000_0000, 0xfff0_0000_0000_0000 },
        .{ 0xbff0_0000_0000_0000, 0x7ff8_0000_0000_0000 },
        .{ 0xc000_0000_0000_0000, 0x7ff8_0000_0000_0000 },
    };

    for (h_cases) |case| try expectTgammaBits(zmath, f16, case[0], case[1]);
    for (f_cases) |case| try expectTgammaBits(zmath, f32, case[0], case[1]);
    for (d_cases) |case| try expectTgammaBits(zmath, f64, case[0], case[1]);
}
