//! Generalized per-lane mapping core for vector (`@Vector`) operations.
//!
//! Every helper here applies the *identical* accepted scalar kernel to each
//! lane, so lane `i` of the result is bit-for-bit equal to the scalar result for
//! lane `i` of the input. That is what preserves the schema-3 `max_ulp = 0`
//! contract and the immutable corpus without new acceptance evidence: the vector
//! path is a per-lane application of the same scalar code, never a re-derived
//! polynomial.
//!
//! Unlike the private helpers in `vector.zig`, the result element type `R` is
//! explicit. The same helpers therefore serve float-to-float, float-to-struct,
//! float-to-integer, float-to-boolean and narrowing (`f64` in, `f32` out)
//! results. The lane type `T` is the accepted scalar kernel's input type, so a
//! vector entry point can check that its lanes match the kernel it selected.
//!
//! Result structs in `scalar_c23_types.zig` are bare `struct` values with auto
//! layout. They are never `extern` or `packed`, so `@bitCast` and `@offsetOf`
//! do not apply. Each field is assigned explicitly, never reinterpreted.
//!
//! Width is inferred from the argument, so one helper covers every
//! `@Vector(N, T)`. A non-vector argument, a mismatched lane type, or mismatched
//! argument types or widths are `@compileError` failures, never silent casts.

// ---------------------------------------------------------------------------
// Shared width, lane and result types
// ---------------------------------------------------------------------------

/// Returns the lane count of the vector type `X`.
fn widthOf(comptime X: type) usize {
    return switch (@typeInfo(X)) {
        .vector => |v| v.len,
        else => @compileError("zmath vector operations require a @Vector argument"),
    };
}

/// Fails unless the vector type `X` has lane type exactly `T`.
fn requireLaneType(comptime X: type, comptime T: type) void {
    switch (@typeInfo(X)) {
        .vector => |v| if (v.child != T) {
            @compileError("zmath vector lane type does not match the selected scalar kernel");
        },
        else => @compileError("zmath vector operations require a @Vector argument"),
    }
}

/// The `@Vector(N, R)` result type for the width of the vector type `X`.
fn ResultVector(comptime R: type, comptime X: anytype) type {
    return @Vector(widthOf(X), R);
}

// ---------------------------------------------------------------------------
// Core per-lane application
// ---------------------------------------------------------------------------

/// Applies the unary scalar kernel `f` to every lane of `x`.
///
/// The result is bit-identical to calling `f` on each lane, because each lane
/// is computed by the same scalar function on the same value.
pub fn mapUnaryR(
    comptime T: type,
    comptime R: type,
    comptime f: fn (T) R,
    x: anytype,
) ResultVector(R, @TypeOf(x)) {
    const n = comptime widthOf(@TypeOf(x));
    requireLaneType(@TypeOf(x), T);
    var out: ResultVector(R, @TypeOf(x)) = undefined;
    inline for (0..n) |i| out[i] = f(x[i]);
    return out;
}

/// Applies the binary scalar kernel `f` to corresponding lanes of `a` and `b`.
pub fn mapBinaryR(
    comptime T: type,
    comptime R: type,
    comptime f: fn (T, T) R,
    a: anytype,
    b: anytype,
) ResultVector(R, @TypeOf(a)) {
    const n = comptime widthOf(@TypeOf(a));
    requireLaneType(@TypeOf(a), T);
    if (@TypeOf(b) != @TypeOf(a)) {
        @compileError("zmath vector binary operations require both arguments to have the same type");
    }
    var out: ResultVector(R, @TypeOf(a)) = undefined;
    inline for (0..n) |i| out[i] = f(a[i], b[i]);
    return out;
}

/// Applies the ternary scalar kernel `f` to lanes `i` of `a`, `b` and `c`.
pub fn mapTernaryR(
    comptime T: type,
    comptime R: type,
    comptime f: fn (T, T, T) R,
    a: anytype,
    b: anytype,
    c: anytype,
) ResultVector(R, @TypeOf(a)) {
    const n = comptime widthOf(@TypeOf(a));
    requireLaneType(@TypeOf(a), T);
    if (@TypeOf(b) != @TypeOf(a)) {
        @compileError("zmath vector ternary operations require all arguments to have the same type");
    }
    if (@TypeOf(c) != @TypeOf(a)) {
        @compileError("zmath vector ternary operations require all arguments to have the same type");
    }
    var out: ResultVector(R, @TypeOf(a)) = undefined;
    inline for (0..n) |i| out[i] = f(a[i], b[i], c[i]);
    return out;
}

/// Applies `f` to every lane of `x` with the extra runtime scalar `p`.
///
/// `p` is a plain runtime scalar, passed unchanged to each lane.
pub fn mapUnaryExtra(
    comptime T: type,
    comptime R: type,
    comptime P: type,
    comptime f: fn (T, P) R,
    p: P,
    x: anytype,
) ResultVector(R, @TypeOf(x)) {
    const n = comptime widthOf(@TypeOf(x));
    requireLaneType(@TypeOf(x), T);
    var out: ResultVector(R, @TypeOf(x)) = undefined;
    inline for (0..n) |i| out[i] = f(x[i], p);
    return out;
}

/// Applies `f` to every lane of `x` with the extra runtime scalars `p` and `q`.
///
/// `p` and `q` are plain runtime scalars, passed unchanged to each lane.
pub fn mapUnaryExtra2(
    comptime T: type,
    comptime R: type,
    comptime P: type,
    comptime Q: type,
    comptime f: fn (T, P, Q) R,
    p: P,
    q: Q,
    x: anytype,
) ResultVector(R, @TypeOf(x)) {
    const n = comptime widthOf(@TypeOf(x));
    requireLaneType(@TypeOf(x), T);
    var out: ResultVector(R, @TypeOf(x)) = undefined;
    inline for (0..n) |i| out[i] = f(x[i], p, q);
    return out;
}

/// Applies the kernel `f` to corresponding lanes of `x` and `e`.
///
/// `e` is a per-lane integer vector, unlike the plain scalars of
/// `mapUnaryExtra`. Its width must equal the width of `x` and its lane type
/// must be exactly `I`.
pub fn mapIntLanes(
    comptime T: type,
    comptime R: type,
    comptime I: type,
    comptime f: fn (T, I) R,
    x: anytype,
    e: anytype,
) ResultVector(R, @TypeOf(x)) {
    const n = comptime widthOf(@TypeOf(x));
    requireLaneType(@TypeOf(x), T);
    requireLaneType(@TypeOf(e), I);
    if (comptime widthOf(@TypeOf(e)) != n) {
        @compileError("zmath vector per-lane operations require both arguments to have the same width");
    }
    var out: ResultVector(R, @TypeOf(x)) = undefined;
    inline for (0..n) |i| out[i] = f(x[i], e[i]);
    return out;
}

comptime {
    _ = mapUnaryR;
    _ = mapBinaryR;
    _ = mapTernaryR;
    _ = mapUnaryExtra;
    _ = mapUnaryExtra2;
    _ = mapIntLanes;
}
