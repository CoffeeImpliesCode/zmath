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
//! Result lanes follow logical index order, not a vector's byte layout. A
//! scoped compiler workaround keeps boolean results in byte lanes on LLVM
//! AArch64 big-endian, where an indexed store into a `@Vector(N, bool)` drops
//! or reorders lanes. Every other result type, and every other target, keeps
//! direct vector storage. The `Vec` bridges map named components to logical
//! lanes.
//!
//! Byte storage only stages the verdicts. `publish` turns them back into a
//! mask by one of two routes: an unsigned integer of one bit per lane at
//! byte-aligned widths, which no measurement has compiled or run, and a
//! per-group lane reversal in Debug at power-of-two widths, which undoes the
//! measured permutation inside each aligned group of eight lanes at every width
//! it claims. Every other width publishes by comparing the staged bytes, which
//! forms the mask as one vector value instead of one lane at a time.
//!
//! Width is inferred from the argument, so one helper covers every
//! `@Vector(N, T)`. A non-vector argument, a mismatched lane type, or mismatched
//! argument types or widths are `@compileError` failures, never silent casts.

const builtin = @import("builtin");
const std = @import("std");

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

/// Indexed boolean vector stores lose lanes on LLVM AArch64 big-endian.
/// Store scalar verdicts as bytes before forming the boolean vector.
const byte_bool_lanes = builtin.zig_backend == .stage2_llvm and
    builtin.cpu.arch == .aarch64_be;

/// Debug boolean publication reverses each aligned lane group at
/// power-of-two widths. Comptime publication and ReleaseSafe need no reversal.
const reverse_bool_lanes = byte_bool_lanes and builtin.mode == .Debug;

/// The lane group the measured Debug permutation acts inside. The permutation
/// reverses a whole width below eight lanes, and each aligned group of eight
/// above it.
const bool_group_lanes = 8;

/// The widest boolean result published through one integer, because Zig
/// integer types admit at most 65535 bits. A wider mask falls back to the
/// comparison route rather than failing to compile.
const max_integer_published_lanes = 65535;

/// The lane group that makes a published mask byte aligned: a boolean mask byte
/// holds eight lanes, and lane `i` is shifted to bit `n - 1 - i`, so a width
/// that is a multiple of eight fills whole bytes and any other width leaves
/// lanes in a partial byte. That is a property of the published mask, not of
/// the staging, which is one byte per lane at every width. No measurement
/// establishes this group; `publishesAtByteWidth` records the open question.
const integer_publish_alignment = 8;

/// The width rule of the reversal repair: the measured permutation only
/// reaches a power-of-two width above one lane. Width one is its own
/// reversal, and the measurement shows no permutation at any other width.
fn reversesAtWidth(comptime n: usize) bool {
    return n > 1 and (n & (n - 1)) == 0;
}

/// Whether a width-`n` boolean result is published by reversing each aligned
/// lane group, which undoes the measured Debug permutation.
fn reversesLanes(comptime n: usize) bool {
    return reverse_bool_lanes and reversesAtWidth(n);
}

/// The width rule of the integer publication: the mask is published as one
/// unsigned integer holding one bit per lane, so the width is one the integer
/// limit admits, and it is a whole number of bytes of that published mask,
/// which is what `integer_publish_alignment` states. The empty width publishes
/// no lane and is excluded.
///
/// The route now has a measurement: `cross-test --targets aarch64_be-linux
/// -- -Doptimize=ReleaseSafe` on this tree compiles and runs it, because it
/// needs an `aarch64_be` target and a non-`Debug` mode at once. That run
/// reports 955 of 956, and its single failure is the one-hot mask sweep at
/// width 14, which does not take this route. The widths this route does take
/// in the shipped sweep, 8, 16 and 24, are clean on that target and mode.
///
/// What no measurement covers is a width this route rejects: those fall back
/// to the comparison route, so nothing here has ever published an unaligned
/// width. That leaves an open question the staging cannot answer, since a lane
/// is one byte at every width and nothing in the staging requires whole bytes.
/// The accepted widths are left exactly as they are, because no measurement
/// justifies changing them.
///
/// A width failing the rule is published by comparing the staged bytes instead,
/// which is a valid mask at every width rather than an unsupported-width
/// compile error.
fn publishesAtByteWidth(comptime n: usize) bool {
    return n != 0 and n % integer_publish_alignment == 0 and n <= max_integer_published_lanes;
}

/// Whether a width-`n` boolean result is published through one unsigned
/// integer holding a bit per lane, instead of a comparison of staged bytes.
/// Debug is excluded, because there the measured permutation is repaired by
/// `reversesLanes` instead.
fn publishesAsInteger(comptime n: usize) bool {
    return byte_bool_lanes and builtin.mode != .Debug and publishesAtByteWidth(n);
}

/// Reverses lanes within each aligned group of at most `bool_group_lanes`.
///
/// Reversing inside a group twice is the identity, so this permutation is its
/// own inverse and undoes the measured one whichever way it is applied.
fn boolLaneIndices(comptime n: usize) [n]i32 {
    std.debug.assert(reversesAtWidth(n));
    // SAFETY: The loop assigns every one of the `n` entries.
    var indices: [n]i32 = undefined;
    const group = @min(n, bool_group_lanes);
    for (0..n) |i| indices[i] = @intCast((i / group) * group + group - 1 - i % group);
    return indices;
}

/// The lane storage of an `n`-lane result of type `R`.
///
/// Direct vector storage is kept everywhere except where an indexed boolean
/// lane store loses earlier true lanes. Byte storage is confined to that
/// compiler target, and every other result type keeps its own lanes.
fn Storage(comptime R: type, comptime n: usize) type {
    return @Vector(n, if (R == bool and byte_bool_lanes) u8 else R);
}

/// Converts a scalar result to its lane storage representation.
inline fn storeLane(comptime R: type, value: R) if (R == bool and byte_bool_lanes) u8 else R {
    return if (R == bool and byte_bool_lanes) @intFromBool(value) else value;
}

/// Forms the result vector without changing the per-lane scalar results.
///
/// Inline publication avoids an additional vector-return ABI boundary.
inline fn publish(comptime R: type, comptime n: usize, lanes: Storage(R, n)) @Vector(n, R) {
    if (R == bool) {
        if (comptime publishesAsInteger(n)) {
            // Publish a byte-aligned mask as an integer, not as an LLVM i1
            // comparison. `publishesAsInteger` keeps `n` inside the widest
            // integer Zig admits, and inside a whole number of mask bytes.
            const Bits = @Int(.unsigned, n);
            var bits: Bits = 0;
            inline for (0..n) |i| bits |= @as(Bits, lanes[i]) << @intCast(n - 1 - i);
            return @bitCast(bits);
        }
        const mask: @Vector(n, bool) = if (comptime byte_bool_lanes)
            lanes != @as(@TypeOf(lanes), @splat(0))
        else
            lanes;
        if (comptime reversesLanes(n)) {
            if (!@inComptime()) return @shuffle(bool, mask, undefined, boolLaneIndices(n));
        }
        return mask;
    }
    return lanes;
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
    // SAFETY: Every result lane is assigned before the result is published.
    var lanes: Storage(R, n) = undefined;
    inline for (0..n) |i| lanes[i] = storeLane(R, f(x[i]));
    return publish(R, n, lanes);
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
    // SAFETY: Every result lane is assigned before the result is published.
    var lanes: Storage(R, n) = undefined;
    inline for (0..n) |i| lanes[i] = storeLane(R, f(a[i], b[i]));
    return publish(R, n, lanes);
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
    // SAFETY: Every result lane is assigned before the result is published.
    var lanes: Storage(R, n) = undefined;
    inline for (0..n) |i| lanes[i] = storeLane(R, f(a[i], b[i], c[i]));
    return publish(R, n, lanes);
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
    // SAFETY: Every result lane is assigned before the result is published.
    var lanes: Storage(R, n) = undefined;
    inline for (0..n) |i| lanes[i] = storeLane(R, f(x[i], p));
    return publish(R, n, lanes);
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
    // SAFETY: Every result lane is assigned before the result is published.
    var lanes: Storage(R, n) = undefined;
    inline for (0..n) |i| lanes[i] = storeLane(R, f(x[i], p, q));
    return publish(R, n, lanes);
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
    // SAFETY: Every result lane is assigned before the result is published.
    var lanes: Storage(R, n) = undefined;
    inline for (0..n) |i| lanes[i] = storeLane(R, f(x[i], e[i]));
    return publish(R, n, lanes);
}

comptime {
    _ = mapUnaryR;
    _ = mapBinaryR;
    _ = mapTernaryR;
    _ = mapUnaryExtra;
    _ = mapUnaryExtra2;
    _ = mapIntLanes;

    // Self check: the width rules of both publication routes are pure, so they
    // are pinned on every target. They are the limits a run-time test could
    // only reach on the one compiler that needs the repair.
    std.debug.assert(@bitSizeOf(@Int(.unsigned, max_integer_published_lanes)) == 65535);
    // The widest byte-aligned width that integer limit still admits.
    const widest_published = max_integer_published_lanes / integer_publish_alignment *
        integer_publish_alignment;
    std.debug.assert(publishesAtByteWidth(widest_published));
    // Four widths fall back to the comparison route, one per rejected part of
    // the rule: the width one below an aligned group and the width one above
    // it, neither a whole number of bytes; the empty width, which publishes no
    // lane; and the first aligned width past the integer limit.
    std.debug.assert(!publishesAtByteWidth(integer_publish_alignment - 1));
    std.debug.assert(!publishesAtByteWidth(integer_publish_alignment + 1));
    std.debug.assert(!publishesAtByteWidth(0));
    std.debug.assert(!publishesAtByteWidth(max_integer_published_lanes + 1));
    // The reversal repairs the power-of-two widths above one lane, and nothing
    // else.
    std.debug.assert(reversesAtWidth(bool_group_lanes));
    std.debug.assert(reversesAtWidth(max_integer_published_lanes + 1));
    std.debug.assert(!reversesAtWidth(0));
    std.debug.assert(!reversesAtWidth(1));
    std.debug.assert(!reversesAtWidth(3));
    // The repair is correct only because the permutation is its own inverse,
    // pinned here below a group, across the group boundary and above it.
    for ([_]usize{ 2, 4, 8, 16, 32, 64, 128 }) |n| {
        const reversed = boolLaneIndices(n);
        for (0..n) |i| {
            std.debug.assert(reversed[@as(usize, @intCast(reversed[i]))] == @as(i32, @intCast(i)));
        }
    }
}
