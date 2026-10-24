const std = @import("std");
const c23 = @import("scalar_c23_types.zig");

fn assertField(
    comptime _: type,
    comptime FieldType: type,
    comptime index: usize,
    comptime name: []const u8,
    comptime fields: anytype,
) void {
    std.debug.assert(std.mem.eql(u8, fields.field_names[index], name));
    std.debug.assert(fields.field_types[index] == FieldType);
}

fn assertResultShapes(comptime T: type) void {
    const Frexp = c23.FrexpResult(T);
    const frexp_fields = @typeInfo(Frexp).@"struct";
    std.debug.assert(frexp_fields.field_names.len == 2);
    assertField(Frexp, T, 0, "fraction", frexp_fields);
    assertField(Frexp, i32, 1, "exponent", frexp_fields);
    const frexp: Frexp = .{ .fraction = @as(T, 0.75), .exponent = -3 };
    std.debug.assert(frexp.fraction == @as(T, 0.75));
    std.debug.assert(frexp.exponent == -3);

    const Modf = c23.ModfResult(T);
    const modf_fields = @typeInfo(Modf).@"struct";
    std.debug.assert(modf_fields.field_names.len == 2);
    assertField(Modf, T, 0, "fractional", modf_fields);
    assertField(Modf, T, 1, "integral", modf_fields);
    const modf: Modf = .{ .fractional = @as(T, -0.5), .integral = @as(T, 3.0) };
    std.debug.assert(modf.fractional == @as(T, -0.5));
    std.debug.assert(modf.integral == @as(T, 3.0));

    const Remquo = c23.RemquoResult(T);
    const remquo_fields = @typeInfo(Remquo).@"struct";
    std.debug.assert(remquo_fields.field_names.len == 2);
    assertField(Remquo, T, 0, "value", remquo_fields);
    assertField(Remquo, i32, 1, "quotient", remquo_fields);
    const remquo: Remquo = .{ .value = @as(T, -0.25), .quotient = -7 };
    std.debug.assert(remquo.value == @as(T, -0.25));
    std.debug.assert(remquo.quotient == -7);

    const Lgamma = c23.LgammaResult(T);
    const lgamma = @typeInfo(Lgamma).@"struct";
    std.debug.assert(lgamma.field_names.len == 2);
    assertField(Lgamma, T, 0, "value", lgamma);
    assertField(Lgamma, i2, 1, "sign", lgamma);
    const lgamma_value: Lgamma = .{ .value = @as(T, 1.5), .sign = -1 };
    std.debug.assert(lgamma_value.value == @as(T, 1.5));
    std.debug.assert(lgamma_value.sign == -1);

    const Atan2 = c23.Atan2Result(T);
    const atan2_fields = @typeInfo(Atan2).@"struct";
    std.debug.assert(atan2_fields.field_names.len == 2);
    assertField(Atan2, T, 0, "value", atan2_fields);
    assertField(Atan2, bool, 1, "resolved", atan2_fields);
    const atan2_resolved: Atan2 = .{ .value = @as(T, 0.5), .resolved = true };
    const atan2_unresolved: Atan2 = .{ .value = @as(T, -0.5), .resolved = false };
    std.debug.assert(atan2_resolved.value == @as(T, 0.5));
    std.debug.assert(atan2_resolved.resolved);
    std.debug.assert(atan2_unresolved.value == @as(T, -0.5));
    std.debug.assert(!atan2_unresolved.resolved);

    const Power = c23.PowerResult(T);
    const power_fields = @typeInfo(Power).@"struct";
    std.debug.assert(power_fields.field_names.len == 2);
    assertField(Power, T, 0, "value", power_fields);
    assertField(Power, bool, 1, "ok", power_fields);
    const power_ok: Power = .{ .value = @as(T, 2), .ok = true };
    const power_not_ok: Power = .{ .value = @as(T, -2), .ok = false };
    std.debug.assert(power_ok.value == @as(T, 2));
    std.debug.assert(power_ok.ok);
    std.debug.assert(power_not_ok.value == @as(T, -2));
    std.debug.assert(!power_not_ok.ok);

    const Payload = c23.PayloadResult(T);
    const payload_fields = @typeInfo(Payload).@"struct";
    std.debug.assert(payload_fields.field_names.len == 2);
    assertField(Payload, T, 0, "value", payload_fields);
    assertField(Payload, bool, 1, "ok", payload_fields);
    const payload: Payload = .{ .value = @as(T, 0.25), .ok = true };
    std.debug.assert(payload.value == @as(T, 0.25));
    std.debug.assert(payload.ok);

    const Canonicalize = c23.CanonicalizeResult(T);
    const canonicalize_fields = @typeInfo(Canonicalize).@"struct";
    std.debug.assert(canonicalize_fields.field_names.len == 2);
    assertField(Canonicalize, T, 0, "value", canonicalize_fields);
    assertField(Canonicalize, bool, 1, "ok", canonicalize_fields);
    const canonicalize: Canonicalize = .{ .value = @as(T, -0.75), .ok = false };
    std.debug.assert(canonicalize.value == @as(T, -0.75));
    std.debug.assert(!canonicalize.ok);
}

test "FromFpRound has five distinct directions" {
    comptime {
        const fields = @typeInfo(c23.FromFpRound).@"enum".field_names;
        const values = [_]c23.FromFpRound{
            .nearest_even,
            .nearest_from_zero,
            .toward_zero,
            .toward_positive,
            .toward_negative,
        };

        std.debug.assert(fields.len == 5);
        std.debug.assert(values.len == fields.len);
        std.debug.assert(std.mem.eql(u8, fields[0], "nearest_even"));
        std.debug.assert(std.mem.eql(u8, fields[1], "nearest_from_zero"));
        std.debug.assert(std.mem.eql(u8, fields[2], "toward_zero"));
        std.debug.assert(std.mem.eql(u8, fields[3], "toward_positive"));
        std.debug.assert(std.mem.eql(u8, fields[4], "toward_negative"));
        std.debug.assert(values[0] != values[1]);
        std.debug.assert(values[0] != values[2]);
        std.debug.assert(values[0] != values[3]);
        std.debug.assert(values[0] != values[4]);
        std.debug.assert(values[1] != values[2]);
        std.debug.assert(values[1] != values[3]);
        std.debug.assert(values[1] != values[4]);
        std.debug.assert(values[2] != values[3]);
        std.debug.assert(values[2] != values[4]);
        std.debug.assert(values[3] != values[4]);
    }
}

test "IntegerResult has exact fields and constructs both states" {
    comptime {
        const fields = @typeInfo(c23.IntegerResult).@"struct";
        std.debug.assert(fields.field_names.len == 2);
        assertField(c23.IntegerResult, i64, 0, "value", fields);
        assertField(c23.IntegerResult, bool, 1, "ok", fields);

        const valid: c23.IntegerResult = .{ .value = 42, .ok = true };
        const invalid: c23.IntegerResult = .{ .value = -7, .ok = false };
        std.debug.assert(valid.value == 42);
        std.debug.assert(valid.ok);
        std.debug.assert(invalid.value == -7);
        std.debug.assert(!invalid.ok);
    }
}

test "C23 result factories preserve f16 f32 and f64 fields" {
    comptime {
        assertResultShapes(f16);
        assertResultShapes(f32);
        assertResultShapes(f64);
    }
}
