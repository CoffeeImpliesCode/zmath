//! Named-component vector value types: `Vec2f`, `Vec3f`, `Vec4f`, `Vec2d`,
//! `Vec3d`, and `Vec4d`.
//!
//! This is the value-type layer that `API.md` gives to `zmath.vector`. It is
//! not the `@Vector` elementwise layer that already lives in `vector.zig`. The
//! two share a module path and nothing else. Every entry point of that layer
//! ends in `v` (`sinfv`, `sinhv`, `hypotv`), takes a `@Vector`, and returns a
//! `@Vector`. Every operation here is a method on a struct with named
//! components, takes and returns a `Vec*`, and carries no `v` suffix.
//!
//! The two surfaces also compute different things. A `@Vector` add is a
//! lane-wise hardware add. `dot` sums in a documented left-to-right order that
//! the compiler may neither reassociate nor contract into an FMA, and `length`
//! composes the accepted `hypot` kernel instead of a naive square root.
//! `toVector` and `fromVector` are the only bridges between the layers, and
//! they move components without making the arithmetic interchangeable.
//!
//! `w` is the fourth component. On `Vec4f` and `Vec4d` it is the homogeneous
//! coordinate, and it matches `Quatf.w` and `Quatd.w`, where `w` is the
//! quaternion's scalar component (`API.md:910`). Nothing here reads `w`
//! differently from the others: `dot`, `length`, and `normalize` include it
//! like any other component.
//!
//! These are fixed-size value types. No operation allocates.
//!
//! ## Operations
//!
//! `init` builds a value from named components. `fromArray` and `toArray`
//! convert to and from `[N]T` in `x, y, z, w` order. `zero` is the additive
//! identity and `one` is the `scaleVec` multiplicative identity; `one` is not
//! a unit vector, because `one().length()` is `sqrt(N)`.
//!
//! `add`, `sub`, `neg`, `scale`, `scaleVec`, and `divScalar` are the
//! arithmetic. `dot` is the scalar product. `cross` is on the three-component
//! types only, because there is no cross product in two or four dimensions
//! and none is invented here. `lengthSquared`, `length`, `distanceSquared`,
//! and `distance` are the magnitudes. `normalize` returns a unit vector, or
//! its input unchanged when the input is zero. `lerp` interpolates for a `t`
//! in `[0, 1]`, which is a precondition and is not enforced. `min`, `max`,
//! `clamp`, and `abs` work componentwise. `isFinite` is true when every
//! component is finite. `eq` is the bit equality below. `toVector` and
//! `fromVector` bridge to `@Vector`.
//!
//! ## Special values
//!
//! `spec/core-math-2026-08-22.md:769-801` pins the encodings, and this module
//! meets it. A NaN operand propagates through the arithmetic with its sign and
//! diagnostic payload intact and its quiet bit set, so a signaling NaN becomes
//! quiet. An invalid operation on non-NaN operands gives the POSITIVE
//! canonical quiet NaN of the format. The second rule is not free: the host
//! unit produces the NEGATIVE indefinite for an invalid operation, 0x7fc00000
//! with the sign bit set, so every arithmetic kernel here checks its result
//! and substitutes the encoding the specification names. That check is a
//! single comparison, and it is what makes the result the same on every
//! target. zmath never reads, validates, sets, or restores `fenv`; the caller
//! must keep round-to-nearest, ties-to-even active.
//!
//! `min`, `max`, and `clamp` are selection rather than arithmetic, so they do
//! not propagate NaN. They are the accepted C23 `fmin` and `fmax` kernels
//! applied componentwise: a NaN component loses to the other operand, and
//! only a pair of NaNs yields a NaN. `clamp` is exactly
//! `min(max(v, lower), upper)`, so it inherits that rule and returns `upper`
//! when `lower` is above `upper`. `abs` is the accepted C23 `abs` kernel,
//! which is a sign-bit mask: it maps `-0` to `+0` and leaves a NaN's payload
//! and quiet bit exactly as they were.
//!
//! `distance` is the norm of the rounded difference, so it is not correctly
//! rounded. Take the difference first when the exact difference matters.
//!
//! `normalize` is defined for a finite nonzero input. A zero vector has no
//! direction and comes back unchanged. An infinite component gives an
//! infinite length, and the division is then the invalid operation
//! `infinity / infinity`, so that component comes back as the positive
//! canonical quiet NaN. A caller who needs a defined result for a non-finite
//! input has to check `isFinite` first.
//!
//! ## Equality
//!
//! `eq` compares bit patterns, with one exception: any NaN equals any other
//! NaN, whatever its sign or payload. A NaN is never equal to a number. `+0`
//! and `-0` are different components and are not equal, and every other pair
//! must match bit for bit. That makes `eq` a total equivalence relation, which
//! is what a caller needs to check a result against a re-derived one. Plain
//! `==` fails the job twice over: it is not reflexive on NaN, and it cannot
//! tell the two zeros apart.

const std = @import("std");
const scalar = @import("scalar.zig");

// ---------------------------------------------------------------------------
// Shared kernels
//
// Every public method below is a one-line delegation to a kernel here. The
// kernels are generic over the value type, so the six types cannot drift apart.
// Each walks the components by name with `inline for`, which unrolls at
// comptime and therefore fixes the order of every reduction.
// ---------------------------------------------------------------------------

/// The component names, in the order the types declare them.
const component_names = [4][]const u8{ "x", "y", "z", "w" };

/// The number of components in `V`.
fn componentCount(comptime V: type) usize {
    return @typeInfo(V).@"struct".fields.len;
}

/// The component type of `V`.
fn Component(comptime V: type) type {
    return @typeInfo(V).@"struct".fields[0].type;
}

/// The unsigned integer of the same width as `T`, for bit comparisons.
fn Bits(comptime T: type) type {
    return switch (@bitSizeOf(T)) {
        32 => u32,
        64 => u64,
        else => @compileError("zmath vector value types are f32 or f64"),
    };
}

/// The accepted `hypot` kernel for `T`.
fn hypotOf(comptime T: type, a: T, b: T) T {
    return switch (T) {
        f32 => scalar.hypotf(a, b),
        f64 => scalar.hypotd(a, b),
        else => @compileError("zmath vector value types are f32 or f64"),
    };
}

/// The accepted C23 `fmin` kernel for `T`.
fn minOf(comptime T: type, a: T, b: T) T {
    return switch (T) {
        f32 => scalar.fminf(a, b),
        f64 => scalar.fmind(a, b),
        else => @compileError("zmath vector value types are f32 or f64"),
    };
}

/// The accepted C23 `fmax` kernel for `T`.
fn maxOf(comptime T: type, a: T, b: T) T {
    return switch (T) {
        f32 => scalar.fmaxf(a, b),
        f64 => scalar.fmaxd(a, b),
        else => @compileError("zmath vector value types are f32 or f64"),
    };
}

/// The accepted C23 `abs` kernel for `T`.
fn absOf(comptime T: type, a: T) T {
    return switch (T) {
        f32 => scalar.absf(a),
        f64 => scalar.absd(a),
        else => @compileError("zmath vector value types are f32 or f64"),
    };
}

/// The quiet bit of `T`, from the masks
/// `spec/core-math-2026-08-22.md:781` gives for F and D.
fn quietMaskOf(comptime T: type) Bits(T) {
    return if (T == f32) 0x0040_0000 else 0x0008_0000_0000_0000;
}

/// The positive canonical quiet NaN of `T`, from
/// `spec/core-math-2026-08-22.md:769-801`.
fn canonicalNaNOf(comptime T: type) T {
    return @bitCast(@as(Bits(T), if (T == f32) 0x7fc0_0000 else 0x7ff8_0000_0000_0000));
}

/// `a` with the quiet bit set, which is a no-op on a quiet NaN. The sign and
/// the diagnostic payload survive, as the specification requires.
fn quietNaNOf(comptime T: type, a: T) T {
    return @bitCast(@as(Bits(T), @bitCast(a)) | quietMaskOf(T));
}

/// The result of an arithmetic operation on `a` and `b`, given the value the
/// host unit produced. A NaN operand keeps its sign and diagnostic payload
/// with the quiet bit set. An invalid operation on non-NaN operands gives the
/// positive canonical quiet NaN, which the host unit does NOT produce on its
/// own: an invalid SSE operation yields the negative indefinite, and the
/// specification pins the positive one.
fn binaryResult(comptime T: type, a: T, b: T, raw: T) T {
    if (!std.math.isNan(raw)) return raw;
    if (std.math.isNan(a)) return quietNaNOf(T, a);
    if (std.math.isNan(b)) return quietNaNOf(T, b);
    return canonicalNaNOf(T);
}

/// The result of an arithmetic operation on `a` alone, given the value the
/// host unit produced.
fn unaryResult(comptime T: type, a: T, raw: T) T {
    if (!std.math.isNan(raw)) return raw;
    if (std.math.isNan(a)) return quietNaNOf(T, a);
    return canonicalNaNOf(T);
}

/// The components of `a` as an array, in declaration order.
fn toArrayOf(comptime V: type, a: V) [componentCount(V)]Component(V) {
    const width = comptime componentCount(V);
    var out: [width]Component(V) = undefined;
    inline for (0..width) |index| {
        out[index] = @field(a, component_names[index]);
    }
    return out;
}

/// The value with every component set to `value`.
fn fillOf(comptime V: type, value: Component(V)) V {
    const width = comptime componentCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        @field(out, component_names[index]) = value;
    }
    return out;
}

fn fromArrayOf(comptime V: type, a: [componentCount(V)]Component(V)) V {
    const width = comptime componentCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        @field(out, component_names[index]) = a[index];
    }
    return out;
}

fn addOf(comptime V: type, a: V, b: V) V {
    @setFloatMode(.strict);
    const T = Component(V);
    const width = comptime componentCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = component_names[index];
        const left = @field(a, name);
        const right = @field(b, name);
        @field(out, name) = binaryResult(T, left, right, left + right);
    }
    return out;
}

fn subOf(comptime V: type, a: V, b: V) V {
    @setFloatMode(.strict);
    const T = Component(V);
    const width = comptime componentCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = component_names[index];
        const left = @field(a, name);
        const right = @field(b, name);
        @field(out, name) = binaryResult(T, left, right, left - right);
    }
    return out;
}

fn negOf(comptime V: type, a: V) V {
    @setFloatMode(.strict);
    const T = Component(V);
    const width = comptime componentCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = component_names[index];
        const component = @field(a, name);
        @field(out, name) = unaryResult(T, component, -component);
    }
    return out;
}

fn scaleOf(comptime V: type, a: V, factor: Component(V)) V {
    @setFloatMode(.strict);
    const T = Component(V);
    const width = comptime componentCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = component_names[index];
        const component = @field(a, name);
        @field(out, name) = binaryResult(T, component, factor, component * factor);
    }
    return out;
}

fn scaleVecOf(comptime V: type, a: V, b: V) V {
    @setFloatMode(.strict);
    const T = Component(V);
    const width = comptime componentCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = component_names[index];
        const left = @field(a, name);
        const right = @field(b, name);
        @field(out, name) = binaryResult(T, left, right, left * right);
    }
    return out;
}

fn divScalarOf(comptime V: type, a: V, divisor: Component(V)) V {
    @setFloatMode(.strict);
    const T = Component(V);
    const width = comptime componentCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = component_names[index];
        const component = @field(a, name);
        @field(out, name) = binaryResult(T, component, divisor, component / divisor);
    }
    return out;
}

/// The scalar product, summed left to right in declaration order. The order is
/// the contract: a caller that groups the terms differently gets a different
/// answer, and `@setFloatMode(.strict)` keeps the compiler from doing the
/// grouping, or contracting a product into the following addition, for them.
fn dotOf(comptime V: type, a: V, b: V) Component(V) {
    @setFloatMode(.strict);
    const T = Component(V);
    const width = comptime componentCount(V);
    var sum: T = 0.0;
    inline for (0..width) |index| {
        const name = component_names[index];
        const left = @field(a, name);
        const right = @field(b, name);
        const product = binaryResult(T, left, right, left * right);
        sum = binaryResult(T, sum, product, sum + product);
    }
    return sum;
}

/// One cross-product component: the difference of two products. Each product
/// and the difference each carry the special-value rule, so an invalid
/// product cannot hide inside a valid difference.
fn crossTermOf(comptime T: type, p: T, q: T, r: T, s: T) T {
    @setFloatMode(.strict);
    const left = binaryResult(T, p, q, p * q);
    const right = binaryResult(T, r, s, r * s);
    return binaryResult(T, left, right, left - right);
}

/// The right-handed cross product. Only the three-component types call this;
/// there is no cross product in two or four dimensions.
fn cross3Of(comptime V: type, a: V, b: V) V {
    @setFloatMode(.strict);
    const T = Component(V);
    return .{
        .x = crossTermOf(T, a.y, b.z, a.z, b.y),
        .y = crossTermOf(T, a.z, b.x, a.x, b.z),
        .z = crossTermOf(T, a.x, b.y, a.y, b.x),
    };
}

/// The Euclidean length, folded through the accepted `hypot` kernel from a
/// positive zero. The fold cannot overflow while the components cannot, and it
/// cannot flush a representable norm to zero, which `sqrt(x*x + y*y)` does in
/// both directions. Nesting the kernel is not correctly rounded for three or
/// four components; it is the accepted kernel applied again, never a
/// re-derived polynomial.
///
/// The first component enters the fold as the value the kernel would have
/// returned for `hypot(+0, x)`, which is that value exactly: a NaN comes back
/// with its sign and payload and its quiet bit set, an infinity comes back
/// positive, and every other value comes back as `|x|`. For finite `x`, binary32
/// uses an exact wider-format square and root; binary64 returns the magnitude
/// bits when the smaller operand is zero. Both routes are exact. The remaining
/// fold is unchanged, and
/// it still decides an infinity against a later signaling or quiet NaN by the
/// kernel's own rules, in the kernel's own order.
fn lengthOf(comptime V: type, a: V) Component(V) {
    const T = Component(V);
    const width = comptime componentCount(V);
    const first = @field(a, component_names[0]);
    var accumulator: T = if (std.math.isNan(first)) quietNaNOf(T, first) else @abs(first);
    inline for (1..width) |index| {
        accumulator = hypotOf(T, accumulator, @field(a, component_names[index]));
    }
    return accumulator;
}

/// The length, or the input unchanged when the length is zero. A zero vector
/// has no direction, and `a / 0` would turn it into a NaN instead of leaving
/// the caller's zero alone, signed zero included.
fn normalizeOf(comptime V: type, a: V) V {
    const length = lengthOf(V, a);
    if (length == 0.0) return a;
    return divScalarOf(V, a, length);
}

fn distanceSquaredOf(comptime V: type, a: V, b: V) Component(V) {
    const difference = subOf(V, a, b);
    return dotOf(V, difference, difference);
}

fn distanceOf(comptime V: type, a: V, b: V) Component(V) {
    return lengthOf(V, subOf(V, a, b));
}

/// Linear interpolation for a `t` in `[0, 1]`. The endpoints are
/// short-circuited because the general form does not reach them: the
/// difference rounds before it is added back, so `a + (b - a) * 1` is not `b`.
fn lerpOf(comptime V: type, a: V, b: V, t: Component(V)) V {
    if (t == 0.0) return a;
    if (t == 1.0) return b;
    return addOf(V, a, scaleOf(V, subOf(V, b, a), t));
}

fn minOfVec(comptime V: type, a: V, b: V) V {
    const T = Component(V);
    const width = comptime componentCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = component_names[index];
        @field(out, name) = minOf(T, @field(a, name), @field(b, name));
    }
    return out;
}

fn maxOfVec(comptime V: type, a: V, b: V) V {
    const T = Component(V);
    const width = comptime componentCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = component_names[index];
        @field(out, name) = maxOf(T, @field(a, name), @field(b, name));
    }
    return out;
}

fn clampOf(comptime V: type, a: V, lower: V, upper: V) V {
    return minOfVec(V, maxOfVec(V, a, lower), upper);
}

fn absOfVec(comptime V: type, a: V) V {
    const T = Component(V);
    const width = comptime componentCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = component_names[index];
        @field(out, name) = absOf(T, @field(a, name));
    }
    return out;
}

fn isFiniteOf(comptime V: type, a: V) bool {
    const width = comptime componentCount(V);
    var finite = true;
    inline for (0..width) |index| {
        finite = finite and std.math.isFinite(@field(a, component_names[index]));
    }
    return finite;
}

/// Bit equality with NaN normalization: two NaNs match, everything else must
/// match bit for bit, so `+0` and `-0` do not.
fn eqComponent(comptime T: type, a: T, b: T) bool {
    // Two NaNs match, whatever their sign or payload. One NaN and one number
    // never match, which is what keeps the relation an equivalence relation on
    // the whole domain rather than on the numbers alone.
    const a_is_nan = std.math.isNan(a);
    if (a_is_nan != std.math.isNan(b)) return false;
    if (a_is_nan) return true;
    return @as(Bits(T), @bitCast(a)) == @as(Bits(T), @bitCast(b));
}

fn eqOf(comptime V: type, a: V, b: V) bool {
    const T = Component(V);
    const width = comptime componentCount(V);
    var same = true;
    inline for (0..width) |index| {
        const name = component_names[index];
        same = same and eqComponent(T, @field(a, name), @field(b, name));
    }
    return same;
}

/// The array and the vector have the same size and alignment, so this is a
/// reinterpretation of the same bits in the same order.
fn toVectorOf(comptime V: type, a: V) @Vector(componentCount(V), Component(V)) {
    return @bitCast(toArrayOf(V, a));
}

fn fromVectorOf(comptime V: type, v: @Vector(componentCount(V), Component(V))) V {
    const array: [componentCount(V)]Component(V) = @bitCast(v);
    return fromArrayOf(V, array);
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// A two-component `f32` vector.
pub const Vec2f = struct {
    /// The first component.
    x: f32,
    /// The second component.
    y: f32,

    /// Builds a value from named components.
    pub fn init(x: f32, y: f32) Vec2f {
        return .{ .x = x, .y = y };
    }

    /// The zero vector, the additive identity.
    pub fn zero() Vec2f {
        return fillOf(Vec2f, 0.0);
    }

    /// The all-ones vector, the `scaleVec` identity and not a unit vector.
    pub fn one() Vec2f {
        return fillOf(Vec2f, 1.0);
    }

    /// Builds a value from `[x, y]`.
    pub fn fromArray(a: [2]f32) Vec2f {
        return fromArrayOf(Vec2f, a);
    }

    /// The components as `[x, y]`.
    pub fn toArray(self: Vec2f) [2]f32 {
        return toArrayOf(Vec2f, self);
    }

    /// Componentwise sum.
    pub fn add(self: Vec2f, other: Vec2f) Vec2f {
        return addOf(Vec2f, self, other);
    }

    /// Componentwise difference.
    pub fn sub(self: Vec2f, other: Vec2f) Vec2f {
        return subOf(Vec2f, self, other);
    }

    /// Componentwise negation, which maps `-0` to `+0`.
    pub fn neg(self: Vec2f) Vec2f {
        return negOf(Vec2f, self);
    }

    /// Every component times `factor`.
    pub fn scale(self: Vec2f, factor: f32) Vec2f {
        return scaleOf(Vec2f, self, factor);
    }

    /// Componentwise product.
    pub fn scaleVec(self: Vec2f, other: Vec2f) Vec2f {
        return scaleVecOf(Vec2f, self, other);
    }

    /// Every component divided by `divisor`.
    pub fn divScalar(self: Vec2f, divisor: f32) Vec2f {
        return divScalarOf(Vec2f, self, divisor);
    }

    /// The scalar product, summed left to right.
    pub fn dot(self: Vec2f, other: Vec2f) f32 {
        return dotOf(Vec2f, self, other);
    }

    /// The squared length, which is the dot product with itself.
    pub fn lengthSquared(self: Vec2f) f32 {
        return dotOf(Vec2f, self, self);
    }

    /// The length, folded through the accepted `hypot` kernel.
    pub fn length(self: Vec2f) f32 {
        return lengthOf(Vec2f, self);
    }

    /// The squared distance to `other`.
    pub fn distanceSquared(self: Vec2f, other: Vec2f) f32 {
        return distanceSquaredOf(Vec2f, self, other);
    }

    /// The distance to `other`, the length of the rounded difference.
    pub fn distance(self: Vec2f, other: Vec2f) f32 {
        return distanceOf(Vec2f, self, other);
    }

    /// The unit vector along `self`, or `self` unchanged when it is zero.
    pub fn normalize(self: Vec2f) Vec2f {
        return normalizeOf(Vec2f, self);
    }

    /// Linear interpolation toward `other` for a `t` in `[0, 1]`.
    pub fn lerp(self: Vec2f, other: Vec2f, t: f32) Vec2f {
        return lerpOf(Vec2f, self, other, t);
    }

    /// Componentwise `fmin`, which prefers `-0` over `+0`.
    pub fn min(self: Vec2f, other: Vec2f) Vec2f {
        return minOfVec(Vec2f, self, other);
    }

    /// Componentwise `fmax`, which prefers `+0` over `-0`.
    pub fn max(self: Vec2f, other: Vec2f) Vec2f {
        return maxOfVec(Vec2f, self, other);
    }

    /// Componentwise `min(max(v, lower), upper)`.
    pub fn clamp(self: Vec2f, lower: Vec2f, upper: Vec2f) Vec2f {
        return clampOf(Vec2f, self, lower, upper);
    }

    /// Componentwise magnitude, which maps `-0` to `+0`.
    pub fn abs(self: Vec2f) Vec2f {
        return absOfVec(Vec2f, self);
    }

    /// True when both components are finite.
    pub fn isFinite(self: Vec2f) bool {
        return isFiniteOf(Vec2f, self);
    }

    /// Bit equality, with every NaN equal to every other NaN.
    pub fn eq(self: Vec2f, other: Vec2f) bool {
        return eqOf(Vec2f, self, other);
    }

    /// The components as `@Vector(2, f32)`, in `x, y` order.
    pub fn toVector(self: Vec2f) @Vector(2, f32) {
        return toVectorOf(Vec2f, self);
    }

    /// The `@Vector(2, f32)` lanes as components.
    pub fn fromVector(v: @Vector(2, f32)) Vec2f {
        return fromVectorOf(Vec2f, v);
    }
};

/// A three-component `f32` vector.
pub const Vec3f = struct {
    /// The first component.
    x: f32,
    /// The second component.
    y: f32,
    /// The third component.
    z: f32,

    /// Builds a value from named components.
    pub fn init(x: f32, y: f32, z: f32) Vec3f {
        return .{ .x = x, .y = y, .z = z };
    }

    /// The zero vector, the additive identity.
    pub fn zero() Vec3f {
        return fillOf(Vec3f, 0.0);
    }

    /// The all-ones vector, the `scaleVec` identity and not a unit vector.
    pub fn one() Vec3f {
        return fillOf(Vec3f, 1.0);
    }

    /// Builds a value from `[x, y, z]`.
    pub fn fromArray(a: [3]f32) Vec3f {
        return fromArrayOf(Vec3f, a);
    }

    /// The components as `[x, y, z]`.
    pub fn toArray(self: Vec3f) [3]f32 {
        return toArrayOf(Vec3f, self);
    }

    /// Componentwise sum.
    pub fn add(self: Vec3f, other: Vec3f) Vec3f {
        return addOf(Vec3f, self, other);
    }

    /// Componentwise difference.
    pub fn sub(self: Vec3f, other: Vec3f) Vec3f {
        return subOf(Vec3f, self, other);
    }

    /// Componentwise negation, which maps `-0` to `+0`.
    pub fn neg(self: Vec3f) Vec3f {
        return negOf(Vec3f, self);
    }

    /// Every component times `factor`.
    pub fn scale(self: Vec3f, factor: f32) Vec3f {
        return scaleOf(Vec3f, self, factor);
    }

    /// Componentwise product.
    pub fn scaleVec(self: Vec3f, other: Vec3f) Vec3f {
        return scaleVecOf(Vec3f, self, other);
    }

    /// Every component divided by `divisor`.
    pub fn divScalar(self: Vec3f, divisor: f32) Vec3f {
        return divScalarOf(Vec3f, self, divisor);
    }

    /// The scalar product, summed left to right.
    pub fn dot(self: Vec3f, other: Vec3f) f32 {
        return dotOf(Vec3f, self, other);
    }

    /// The right-handed cross product.
    pub fn cross(self: Vec3f, other: Vec3f) Vec3f {
        return cross3Of(Vec3f, self, other);
    }

    /// The squared length, which is the dot product with itself.
    pub fn lengthSquared(self: Vec3f) f32 {
        return dotOf(Vec3f, self, self);
    }

    /// The length, folded through the accepted `hypot` kernel.
    pub fn length(self: Vec3f) f32 {
        return lengthOf(Vec3f, self);
    }

    /// The squared distance to `other`.
    pub fn distanceSquared(self: Vec3f, other: Vec3f) f32 {
        return distanceSquaredOf(Vec3f, self, other);
    }

    /// The distance to `other`, the length of the rounded difference.
    pub fn distance(self: Vec3f, other: Vec3f) f32 {
        return distanceOf(Vec3f, self, other);
    }

    /// The unit vector along `self`, or `self` unchanged when it is zero.
    pub fn normalize(self: Vec3f) Vec3f {
        return normalizeOf(Vec3f, self);
    }

    /// Linear interpolation toward `other` for a `t` in `[0, 1]`.
    pub fn lerp(self: Vec3f, other: Vec3f, t: f32) Vec3f {
        return lerpOf(Vec3f, self, other, t);
    }

    /// Componentwise `fmin`, which prefers `-0` over `+0`.
    pub fn min(self: Vec3f, other: Vec3f) Vec3f {
        return minOfVec(Vec3f, self, other);
    }

    /// Componentwise `fmax`, which prefers `+0` over `-0`.
    pub fn max(self: Vec3f, other: Vec3f) Vec3f {
        return maxOfVec(Vec3f, self, other);
    }

    /// Componentwise `min(max(v, lower), upper)`.
    pub fn clamp(self: Vec3f, lower: Vec3f, upper: Vec3f) Vec3f {
        return clampOf(Vec3f, self, lower, upper);
    }

    /// Componentwise magnitude, which maps `-0` to `+0`.
    pub fn abs(self: Vec3f) Vec3f {
        return absOfVec(Vec3f, self);
    }

    /// True when every component is finite.
    pub fn isFinite(self: Vec3f) bool {
        return isFiniteOf(Vec3f, self);
    }

    /// Bit equality, with every NaN equal to every other NaN.
    pub fn eq(self: Vec3f, other: Vec3f) bool {
        return eqOf(Vec3f, self, other);
    }

    /// The components as `@Vector(3, f32)`, in `x, y, z` order.
    pub fn toVector(self: Vec3f) @Vector(3, f32) {
        return toVectorOf(Vec3f, self);
    }

    /// The `@Vector(3, f32)` lanes as components.
    pub fn fromVector(v: @Vector(3, f32)) Vec3f {
        return fromVectorOf(Vec3f, v);
    }
};

/// A four-component `f32` vector. `w` is the homogeneous component.
pub const Vec4f = struct {
    /// The first component.
    x: f32,
    /// The second component.
    y: f32,
    /// The third component.
    z: f32,
    /// The fourth, homogeneous component, matching `Quatf.w`.
    w: f32,

    /// Builds a value from named components.
    pub fn init(x: f32, y: f32, z: f32, w: f32) Vec4f {
        return .{ .x = x, .y = y, .z = z, .w = w };
    }

    /// The zero vector, the additive identity.
    pub fn zero() Vec4f {
        return fillOf(Vec4f, 0.0);
    }

    /// The all-ones vector, the `scaleVec` identity and not a unit vector.
    pub fn one() Vec4f {
        return fillOf(Vec4f, 1.0);
    }

    /// Builds a value from `[x, y, z, w]`.
    pub fn fromArray(a: [4]f32) Vec4f {
        return fromArrayOf(Vec4f, a);
    }

    /// The components as `[x, y, z, w]`.
    pub fn toArray(self: Vec4f) [4]f32 {
        return toArrayOf(Vec4f, self);
    }

    /// Componentwise sum.
    pub fn add(self: Vec4f, other: Vec4f) Vec4f {
        return addOf(Vec4f, self, other);
    }

    /// Componentwise difference.
    pub fn sub(self: Vec4f, other: Vec4f) Vec4f {
        return subOf(Vec4f, self, other);
    }

    /// Componentwise negation, which maps `-0` to `+0`.
    pub fn neg(self: Vec4f) Vec4f {
        return negOf(Vec4f, self);
    }

    /// Every component times `factor`.
    pub fn scale(self: Vec4f, factor: f32) Vec4f {
        return scaleOf(Vec4f, self, factor);
    }

    /// Componentwise product.
    pub fn scaleVec(self: Vec4f, other: Vec4f) Vec4f {
        return scaleVecOf(Vec4f, self, other);
    }

    /// Every component divided by `divisor`.
    pub fn divScalar(self: Vec4f, divisor: f32) Vec4f {
        return divScalarOf(Vec4f, self, divisor);
    }

    /// The scalar product, summed left to right.
    pub fn dot(self: Vec4f, other: Vec4f) f32 {
        return dotOf(Vec4f, self, other);
    }

    /// The squared length, which is the dot product with itself.
    pub fn lengthSquared(self: Vec4f) f32 {
        return dotOf(Vec4f, self, self);
    }

    /// The length, folded through the accepted `hypot` kernel.
    pub fn length(self: Vec4f) f32 {
        return lengthOf(Vec4f, self);
    }

    /// The squared distance to `other`.
    pub fn distanceSquared(self: Vec4f, other: Vec4f) f32 {
        return distanceSquaredOf(Vec4f, self, other);
    }

    /// The distance to `other`, the length of the rounded difference.
    pub fn distance(self: Vec4f, other: Vec4f) f32 {
        return distanceOf(Vec4f, self, other);
    }

    /// The unit vector along `self`, or `self` unchanged when it is zero.
    pub fn normalize(self: Vec4f) Vec4f {
        return normalizeOf(Vec4f, self);
    }

    /// Linear interpolation toward `other` for a `t` in `[0, 1]`.
    pub fn lerp(self: Vec4f, other: Vec4f, t: f32) Vec4f {
        return lerpOf(Vec4f, self, other, t);
    }

    /// Componentwise `fmin`, which prefers `-0` over `+0`.
    pub fn min(self: Vec4f, other: Vec4f) Vec4f {
        return minOfVec(Vec4f, self, other);
    }

    /// Componentwise `fmax`, which prefers `+0` over `-0`.
    pub fn max(self: Vec4f, other: Vec4f) Vec4f {
        return maxOfVec(Vec4f, self, other);
    }

    /// Componentwise `min(max(v, lower), upper)`.
    pub fn clamp(self: Vec4f, lower: Vec4f, upper: Vec4f) Vec4f {
        return clampOf(Vec4f, self, lower, upper);
    }

    /// Componentwise magnitude, which maps `-0` to `+0`.
    pub fn abs(self: Vec4f) Vec4f {
        return absOfVec(Vec4f, self);
    }

    /// True when every component is finite.
    pub fn isFinite(self: Vec4f) bool {
        return isFiniteOf(Vec4f, self);
    }

    /// Bit equality, with every NaN equal to every other NaN.
    pub fn eq(self: Vec4f, other: Vec4f) bool {
        return eqOf(Vec4f, self, other);
    }

    /// The components as `@Vector(4, f32)`, in `x, y, z, w` order.
    pub fn toVector(self: Vec4f) @Vector(4, f32) {
        return toVectorOf(Vec4f, self);
    }

    /// The `@Vector(4, f32)` lanes as components.
    pub fn fromVector(v: @Vector(4, f32)) Vec4f {
        return fromVectorOf(Vec4f, v);
    }
};

/// A two-component `f64` vector.
pub const Vec2d = struct {
    /// The first component.
    x: f64,
    /// The second component.
    y: f64,

    /// Builds a value from named components.
    pub fn init(x: f64, y: f64) Vec2d {
        return .{ .x = x, .y = y };
    }

    /// The zero vector, the additive identity.
    pub fn zero() Vec2d {
        return fillOf(Vec2d, 0.0);
    }

    /// The all-ones vector, the `scaleVec` identity and not a unit vector.
    pub fn one() Vec2d {
        return fillOf(Vec2d, 1.0);
    }

    /// Builds a value from `[x, y]`.
    pub fn fromArray(a: [2]f64) Vec2d {
        return fromArrayOf(Vec2d, a);
    }

    /// The components as `[x, y]`.
    pub fn toArray(self: Vec2d) [2]f64 {
        return toArrayOf(Vec2d, self);
    }

    /// Componentwise sum.
    pub fn add(self: Vec2d, other: Vec2d) Vec2d {
        return addOf(Vec2d, self, other);
    }

    /// Componentwise difference.
    pub fn sub(self: Vec2d, other: Vec2d) Vec2d {
        return subOf(Vec2d, self, other);
    }

    /// Componentwise negation, which maps `-0` to `+0`.
    pub fn neg(self: Vec2d) Vec2d {
        return negOf(Vec2d, self);
    }

    /// Every component times `factor`.
    pub fn scale(self: Vec2d, factor: f64) Vec2d {
        return scaleOf(Vec2d, self, factor);
    }

    /// Componentwise product.
    pub fn scaleVec(self: Vec2d, other: Vec2d) Vec2d {
        return scaleVecOf(Vec2d, self, other);
    }

    /// Every component divided by `divisor`.
    pub fn divScalar(self: Vec2d, divisor: f64) Vec2d {
        return divScalarOf(Vec2d, self, divisor);
    }

    /// The scalar product, summed left to right.
    pub fn dot(self: Vec2d, other: Vec2d) f64 {
        return dotOf(Vec2d, self, other);
    }

    /// The squared length, which is the dot product with itself.
    pub fn lengthSquared(self: Vec2d) f64 {
        return dotOf(Vec2d, self, self);
    }

    /// The length, folded through the accepted `hypot` kernel.
    pub fn length(self: Vec2d) f64 {
        return lengthOf(Vec2d, self);
    }

    /// The squared distance to `other`.
    pub fn distanceSquared(self: Vec2d, other: Vec2d) f64 {
        return distanceSquaredOf(Vec2d, self, other);
    }

    /// The distance to `other`, the length of the rounded difference.
    pub fn distance(self: Vec2d, other: Vec2d) f64 {
        return distanceOf(Vec2d, self, other);
    }

    /// The unit vector along `self`, or `self` unchanged when it is zero.
    pub fn normalize(self: Vec2d) Vec2d {
        return normalizeOf(Vec2d, self);
    }

    /// Linear interpolation toward `other` for a `t` in `[0, 1]`.
    pub fn lerp(self: Vec2d, other: Vec2d, t: f64) Vec2d {
        return lerpOf(Vec2d, self, other, t);
    }

    /// Componentwise `fmin`, which prefers `-0` over `+0`.
    pub fn min(self: Vec2d, other: Vec2d) Vec2d {
        return minOfVec(Vec2d, self, other);
    }

    /// Componentwise `fmax`, which prefers `+0` over `-0`.
    pub fn max(self: Vec2d, other: Vec2d) Vec2d {
        return maxOfVec(Vec2d, self, other);
    }

    /// Componentwise `min(max(v, lower), upper)`.
    pub fn clamp(self: Vec2d, lower: Vec2d, upper: Vec2d) Vec2d {
        return clampOf(Vec2d, self, lower, upper);
    }

    /// Componentwise magnitude, which maps `-0` to `+0`.
    pub fn abs(self: Vec2d) Vec2d {
        return absOfVec(Vec2d, self);
    }

    /// True when both components are finite.
    pub fn isFinite(self: Vec2d) bool {
        return isFiniteOf(Vec2d, self);
    }

    /// Bit equality, with every NaN equal to every other NaN.
    pub fn eq(self: Vec2d, other: Vec2d) bool {
        return eqOf(Vec2d, self, other);
    }

    /// The components as `@Vector(2, f64)`, in `x, y` order.
    pub fn toVector(self: Vec2d) @Vector(2, f64) {
        return toVectorOf(Vec2d, self);
    }

    /// The `@Vector(2, f64)` lanes as components.
    pub fn fromVector(v: @Vector(2, f64)) Vec2d {
        return fromVectorOf(Vec2d, v);
    }
};

/// A three-component `f64` vector.
pub const Vec3d = struct {
    /// The first component.
    x: f64,
    /// The second component.
    y: f64,
    /// The third component.
    z: f64,

    /// Builds a value from named components.
    pub fn init(x: f64, y: f64, z: f64) Vec3d {
        return .{ .x = x, .y = y, .z = z };
    }

    /// The zero vector, the additive identity.
    pub fn zero() Vec3d {
        return fillOf(Vec3d, 0.0);
    }

    /// The all-ones vector, the `scaleVec` identity and not a unit vector.
    pub fn one() Vec3d {
        return fillOf(Vec3d, 1.0);
    }

    /// Builds a value from `[x, y, z]`.
    pub fn fromArray(a: [3]f64) Vec3d {
        return fromArrayOf(Vec3d, a);
    }

    /// The components as `[x, y, z]`.
    pub fn toArray(self: Vec3d) [3]f64 {
        return toArrayOf(Vec3d, self);
    }

    /// Componentwise sum.
    pub fn add(self: Vec3d, other: Vec3d) Vec3d {
        return addOf(Vec3d, self, other);
    }

    /// Componentwise difference.
    pub fn sub(self: Vec3d, other: Vec3d) Vec3d {
        return subOf(Vec3d, self, other);
    }

    /// Componentwise negation, which maps `-0` to `+0`.
    pub fn neg(self: Vec3d) Vec3d {
        return negOf(Vec3d, self);
    }

    /// Every component times `factor`.
    pub fn scale(self: Vec3d, factor: f64) Vec3d {
        return scaleOf(Vec3d, self, factor);
    }

    /// Componentwise product.
    pub fn scaleVec(self: Vec3d, other: Vec3d) Vec3d {
        return scaleVecOf(Vec3d, self, other);
    }

    /// Every component divided by `divisor`.
    pub fn divScalar(self: Vec3d, divisor: f64) Vec3d {
        return divScalarOf(Vec3d, self, divisor);
    }

    /// The scalar product, summed left to right.
    pub fn dot(self: Vec3d, other: Vec3d) f64 {
        return dotOf(Vec3d, self, other);
    }

    /// The right-handed cross product.
    pub fn cross(self: Vec3d, other: Vec3d) Vec3d {
        return cross3Of(Vec3d, self, other);
    }

    /// The squared length, which is the dot product with itself.
    pub fn lengthSquared(self: Vec3d) f64 {
        return dotOf(Vec3d, self, self);
    }

    /// The length, folded through the accepted `hypot` kernel.
    pub fn length(self: Vec3d) f64 {
        return lengthOf(Vec3d, self);
    }

    /// The squared distance to `other`.
    pub fn distanceSquared(self: Vec3d, other: Vec3d) f64 {
        return distanceSquaredOf(Vec3d, self, other);
    }

    /// The distance to `other`, the length of the rounded difference.
    pub fn distance(self: Vec3d, other: Vec3d) f64 {
        return distanceOf(Vec3d, self, other);
    }

    /// The unit vector along `self`, or `self` unchanged when it is zero.
    pub fn normalize(self: Vec3d) Vec3d {
        return normalizeOf(Vec3d, self);
    }

    /// Linear interpolation toward `other` for a `t` in `[0, 1]`.
    pub fn lerp(self: Vec3d, other: Vec3d, t: f64) Vec3d {
        return lerpOf(Vec3d, self, other, t);
    }

    /// Componentwise `fmin`, which prefers `-0` over `+0`.
    pub fn min(self: Vec3d, other: Vec3d) Vec3d {
        return minOfVec(Vec3d, self, other);
    }

    /// Componentwise `fmax`, which prefers `+0` over `-0`.
    pub fn max(self: Vec3d, other: Vec3d) Vec3d {
        return maxOfVec(Vec3d, self, other);
    }

    /// Componentwise `min(max(v, lower), upper)`.
    pub fn clamp(self: Vec3d, lower: Vec3d, upper: Vec3d) Vec3d {
        return clampOf(Vec3d, self, lower, upper);
    }

    /// Componentwise magnitude, which maps `-0` to `+0`.
    pub fn abs(self: Vec3d) Vec3d {
        return absOfVec(Vec3d, self);
    }

    /// True when every component is finite.
    pub fn isFinite(self: Vec3d) bool {
        return isFiniteOf(Vec3d, self);
    }

    /// Bit equality, with every NaN equal to every other NaN.
    pub fn eq(self: Vec3d, other: Vec3d) bool {
        return eqOf(Vec3d, self, other);
    }

    /// The components as `@Vector(3, f64)`, in `x, y, z` order.
    pub fn toVector(self: Vec3d) @Vector(3, f64) {
        return toVectorOf(Vec3d, self);
    }

    /// The `@Vector(3, f64)` lanes as components.
    pub fn fromVector(v: @Vector(3, f64)) Vec3d {
        return fromVectorOf(Vec3d, v);
    }
};

/// A four-component `f64` vector. `w` is the homogeneous component.
pub const Vec4d = struct {
    /// The first component.
    x: f64,
    /// The second component.
    y: f64,
    /// The third component.
    z: f64,
    /// The fourth, homogeneous component, matching `Quatd.w`.
    w: f64,

    /// Builds a value from named components.
    pub fn init(x: f64, y: f64, z: f64, w: f64) Vec4d {
        return .{ .x = x, .y = y, .z = z, .w = w };
    }

    /// The zero vector, the additive identity.
    pub fn zero() Vec4d {
        return fillOf(Vec4d, 0.0);
    }

    /// The all-ones vector, the `scaleVec` identity and not a unit vector.
    pub fn one() Vec4d {
        return fillOf(Vec4d, 1.0);
    }

    /// Builds a value from `[x, y, z, w]`.
    pub fn fromArray(a: [4]f64) Vec4d {
        return fromArrayOf(Vec4d, a);
    }

    /// The components as `[x, y, z, w]`.
    pub fn toArray(self: Vec4d) [4]f64 {
        return toArrayOf(Vec4d, self);
    }

    /// Componentwise sum.
    pub fn add(self: Vec4d, other: Vec4d) Vec4d {
        return addOf(Vec4d, self, other);
    }

    /// Componentwise difference.
    pub fn sub(self: Vec4d, other: Vec4d) Vec4d {
        return subOf(Vec4d, self, other);
    }

    /// Componentwise negation, which maps `-0` to `+0`.
    pub fn neg(self: Vec4d) Vec4d {
        return negOf(Vec4d, self);
    }

    /// Every component times `factor`.
    pub fn scale(self: Vec4d, factor: f64) Vec4d {
        return scaleOf(Vec4d, self, factor);
    }

    /// Componentwise product.
    pub fn scaleVec(self: Vec4d, other: Vec4d) Vec4d {
        return scaleVecOf(Vec4d, self, other);
    }

    /// Every component divided by `divisor`.
    pub fn divScalar(self: Vec4d, divisor: f64) Vec4d {
        return divScalarOf(Vec4d, self, divisor);
    }

    /// The scalar product, summed left to right.
    pub fn dot(self: Vec4d, other: Vec4d) f64 {
        return dotOf(Vec4d, self, other);
    }

    /// The squared length, which is the dot product with itself.
    pub fn lengthSquared(self: Vec4d) f64 {
        return dotOf(Vec4d, self, self);
    }

    /// The length, folded through the accepted `hypot` kernel.
    pub fn length(self: Vec4d) f64 {
        return lengthOf(Vec4d, self);
    }

    /// The squared distance to `other`.
    pub fn distanceSquared(self: Vec4d, other: Vec4d) f64 {
        return distanceSquaredOf(Vec4d, self, other);
    }

    /// The distance to `other`, the length of the rounded difference.
    pub fn distance(self: Vec4d, other: Vec4d) f64 {
        return distanceOf(Vec4d, self, other);
    }

    /// The unit vector along `self`, or `self` unchanged when it is zero.
    pub fn normalize(self: Vec4d) Vec4d {
        return normalizeOf(Vec4d, self);
    }

    /// Linear interpolation toward `other` for a `t` in `[0, 1]`.
    pub fn lerp(self: Vec4d, other: Vec4d, t: f64) Vec4d {
        return lerpOf(Vec4d, self, other, t);
    }

    /// Componentwise `fmin`, which prefers `-0` over `+0`.
    pub fn min(self: Vec4d, other: Vec4d) Vec4d {
        return minOfVec(Vec4d, self, other);
    }

    /// Componentwise `fmax`, which prefers `+0` over `-0`.
    pub fn max(self: Vec4d, other: Vec4d) Vec4d {
        return maxOfVec(Vec4d, self, other);
    }

    /// Componentwise `min(max(v, lower), upper)`.
    pub fn clamp(self: Vec4d, lower: Vec4d, upper: Vec4d) Vec4d {
        return clampOf(Vec4d, self, lower, upper);
    }

    /// Componentwise magnitude, which maps `-0` to `+0`.
    pub fn abs(self: Vec4d) Vec4d {
        return absOfVec(Vec4d, self);
    }

    /// True when every component is finite.
    pub fn isFinite(self: Vec4d) bool {
        return isFiniteOf(Vec4d, self);
    }

    /// Bit equality, with every NaN equal to every other NaN.
    pub fn eq(self: Vec4d, other: Vec4d) bool {
        return eqOf(Vec4d, self, other);
    }

    /// The components as `@Vector(4, f64)`, in `x, y, z, w` order.
    pub fn toVector(self: Vec4d) @Vector(4, f64) {
        return toVectorOf(Vec4d, self);
    }

    /// The `@Vector(4, f64)` lanes as components.
    pub fn fromVector(v: @Vector(4, f64)) Vec4d {
        return fromVectorOf(Vec4d, v);
    }
};
