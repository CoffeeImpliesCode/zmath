//! Hamilton quaternions, and the validated unit quaternions that stand for
//! rotations.
//!
//! `Quatf` and `Quatd` are the raw algebra. `UnitQuatf` and `UnitQuatd` wrap
//! one of them after `init` has proved it is unit, and they are the only types
//! that can rotate a vector. The split exists because "a quaternion" and "a
//! rotation" are different claims: a rotation is only a rotation when the
//! quaternion is unit, and `q * (v, 0) * conjugate(q)` for a non-unit `q`
//! scales `v` by `normSquared(q)` instead of leaving its length alone. Putting
//! the proof in the type is what keeps that mistake off the API.
//!
//! `w` is the scalar component and `x`, `y`, `z` are the vector part, which is
//! the order `API.md:901-905` spells the fields in. Hamilton multiplication is
//! written so the scalar part is the real part of the four-dimensional dot
//! product, so `a.mul(a.conjugate())` is `normSquared(a) * identity()` and
//! nothing else. Every product here is summed in the order `API.md:907` gives,
//! with `@setFloatMode(.strict)`, because a caller that gets a different
//! grouping gets a different answer and the difference is not always harmless.
//!
//! ## Composition
//!
//! `a.then(b)` means apply `a` first and then `b`, and it is `b.mul(a)`. That
//! reversal is the whole reason `then` exists: composition of functions reads
//! left to right, composition of quaternions reads right to left, and a
//! rotation pipeline written the first way and built the second way is a sign
//! error with no test failure to warn about it. `rotate` is the active,
//! right-handed action on a column vector, so `a.then(b).rotate(v)` is
//! `b.rotate(a.rotate(v))` and the two agree to rounding.
//!
//! ## Special values
//!
//! `spec/core-math-2026-08-22.md:769-801` pins the encodings and this module
//! meets it. A NaN anywhere in an operand keeps its sign and diagnostic
//! payload with the quiet bit set, so a signaling NaN becomes quiet. An invalid
//! operation on non-NaN operands gives the POSITIVE canonical quiet NaN, which
//! the host unit does not produce on its own, so every result is checked.
//! zmath never reads, validates, sets, or restores `fenv`; the caller must keep
//! round-to-nearest, ties-to-even active.
//!
//! Every entry point that can be given a value with no meaningful answer
//! returns an optional. `Quatf.inverse` returns `null` for a zero or
//! non-finite norm, `UnitQuatf.init` returns `null` for a value it cannot make
//! unit, `slerp` returns `null` for a non-finite `t`, and `fromAxisAngle`
//! returns `null` for a zero or non-finite axis or a non-finite angle. The
//! algebraic operations do not return optionals: a quaternion is closed under
//! them, and a NaN component is the answer.
//!
//! These are fixed-size value types. No operation allocates.

const std = @import("std");
const scalar = @import("scalar.zig");
const vector_value = @import("vector_value.zig");

const Vec3f = vector_value.Vec3f;
const Vec3d = vector_value.Vec3d;

// ---------------------------------------------------------------------------
// Kernels
//
// Every public method below is a one-line delegation to a kernel here. The
// kernels are generic over the component type, so `f32` and `f64` cannot drift
// apart, and they index by position in `x, y, z, w` order so every reduction
// and product order is fixed at comptime instead of written out twice.
// ---------------------------------------------------------------------------

/// The unsigned integer of the same width as `T`, for bit comparisons.
fn Bits(comptime T: type) type {
    return switch (@bitSizeOf(T)) {
        32 => u32,
        64 => u64,
        else => @compileError("zmath quaternions are f32 or f64"),
    };
}

/// The accepted `hypot` kernel for `T`.
fn hypotOf(comptime T: type, a: T, b: T) T {
    return switch (T) {
        f32 => scalar.hypotf(a, b),
        f64 => scalar.hypotd(a, b),
        else => @compileError("zmath quaternions are f32 or f64"),
    };
}

fn sinOf(comptime T: type, value: T) T {
    return switch (T) {
        f32 => scalar.sinf(value),
        f64 => scalar.sind(value),
        else => @compileError("zmath quaternions are f32 or f64"),
    };
}

fn cosOf(comptime T: type, value: T) T {
    return switch (T) {
        f32 => scalar.cosf(value),
        f64 => scalar.cosd(value),
        else => @compileError("zmath quaternions are f32 or f64"),
    };
}

fn acosOf(comptime T: type, value: T) T {
    return switch (T) {
        f32 => scalar.acosf(value),
        f64 => scalar.acosd(value),
        else => @compileError("zmath quaternions are f32 or f64"),
    };
}

/// The accepted C23 `fmin` and `fmax` kernels for `T`. They are selection and
/// not arithmetic, so a NaN operand loses, which is what lets them clamp a dot
/// product that rounded a hair outside the range it belongs in.
fn minOf(comptime T: type, a: T, b: T) T {
    return switch (T) {
        f32 => scalar.fminf(a, b),
        f64 => scalar.fmind(a, b),
        else => @compileError("zmath quaternions are f32 or f64"),
    };
}

fn maxOf(comptime T: type, a: T, b: T) T {
    return switch (T) {
        f32 => scalar.fmaxf(a, b),
        f64 => scalar.fmaxd(a, b),
        else => @compileError("zmath quaternions are f32 or f64"),
    };
}

/// The quiet bit of `T`, from the masks `spec/core-math-2026-08-22.md:781`
/// gives for F and D.
fn quietMaskOf(comptime T: type) Bits(T) {
    return if (T == f32) 0x0040_0000 else 0x0008_0000_0000_0000;
}

/// The positive canonical quiet NaN of `T`, from
/// `spec/core-math-2026-08-22.md:769-801`.
fn canonicalNaNOf(comptime T: type) T {
    return @bitCast(@as(Bits(T), if (T == f32) 0x7fc0_0000 else 0x7ff8_0000_0000_0000));
}

/// `value` with the quiet bit set, which is a no-op on a quiet NaN. The sign
/// and the diagnostic payload survive, as the specification requires.
fn quietNaNOf(comptime T: type, value: T) T {
    return @bitCast(@as(Bits(T), @bitCast(value)) | quietMaskOf(T));
}

/// The special-value rule for a result that consumed every component of both
/// `a` and `b`. A NaN anywhere in either operand keeps its sign and diagnostic
/// payload with the quiet bit set, so the result never reports a payload that
/// came from nowhere. An invalid operation on non-NaN operands gives the
/// positive canonical quiet NaN, which the host unit does NOT produce on its
/// own: an invalid SSE operation yields the negative indefinite, and the
/// specification pins the positive one. The scan runs `a` then `b`, and within
/// each in declaration order, so when both operands carry a NaN the survivor
/// is the first in that fixed order rather than whatever the hardware found
/// first.
fn binaryResult(comptime T: type, a: [4]T, b: [4]T, raw: T) T {
    if (!std.math.isNan(raw)) return raw;
    inline for (a) |component| {
        if (std.math.isNan(component)) return quietNaNOf(T, component);
    }
    inline for (b) |component| {
        if (std.math.isNan(component)) return quietNaNOf(T, component);
    }
    return canonicalNaNOf(T);
}

/// The special-value rule for a result that consumed every component of `a`.
fn unaryResult(comptime T: type, a: [4]T, raw: T) T {
    if (!std.math.isNan(raw)) return raw;
    inline for (a) |component| {
        if (std.math.isNan(component)) return quietNaNOf(T, component);
    }
    return canonicalNaNOf(T);
}

/// The components as `[x, y, z, w]`, in the order `API.md:901` spells them.
fn toArrayOf(comptime Q: type, a: Q) [4]Component(Q) {
    return .{ a.x, a.y, a.z, a.w };
}

/// The component type of a raw or validated quaternion.
fn Component(comptime Q: type) type {
    return @typeInfo(Q).@"struct".fields[0].type;
}

fn fromArrayOf(comptime Q: type, a: [4]Component(Q)) Q {
    return .{ .x = a[0], .y = a[1], .z = a[2], .w = a[3] };
}

fn addOf(comptime T: type, a: [4]T, b: [4]T) [4]T {
    @setFloatMode(.strict);
    var out: [4]T = undefined;
    inline for (0..4) |index| {
        out[index] = binaryResult(T, a, b, a[index] + b[index]);
    }
    return out;
}

fn subOf(comptime T: type, a: [4]T, b: [4]T) [4]T {
    @setFloatMode(.strict);
    var out: [4]T = undefined;
    inline for (0..4) |index| {
        out[index] = binaryResult(T, a, b, a[index] - b[index]);
    }
    return out;
}

fn negOf(comptime T: type, a: [4]T) [4]T {
    @setFloatMode(.strict);
    var out: [4]T = undefined;
    inline for (0..4) |index| {
        out[index] = unaryResult(T, a, -a[index]);
    }
    return out;
}

fn scaleOf(comptime T: type, a: [4]T, factor: T) [4]T {
    @setFloatMode(.strict);
    const operands = [4]T{ factor, factor, factor, factor };
    var out: [4]T = undefined;
    inline for (0..4) |index| {
        out[index] = binaryResult(T, a, operands, a[index] * factor);
    }
    return out;
}

/// The componentwise product. It is not a quaternion product: `x` multiplies
/// `x` and nothing rotates.
fn scaleQuatOf(comptime T: type, a: [4]T, b: [4]T) [4]T {
    @setFloatMode(.strict);
    var out: [4]T = undefined;
    inline for (0..4) |index| {
        out[index] = binaryResult(T, a, b, a[index] * b[index]);
    }
    return out;
}

/// Every component divided by `divisor`. The caller has already rejected a
/// zero and a non-finite divisor, so a NaN here can only come from a component
/// of `a` and the rule is the unary one.
fn divByOf(comptime T: type, a: [4]T, divisor: T) [4]T {
    @setFloatMode(.strict);
    var out: [4]T = undefined;
    inline for (0..4) |index| {
        out[index] = unaryResult(T, a, a[index] / divisor);
    }
    return out;
}

/// The scalar product, summed left to right in `x, y, z, w` order. The order
/// is the contract, and `@setFloatMode(.strict)` keeps the compiler from
/// reassociating the sum or contracting a product into the following addition.
fn dotOf(comptime T: type, a: [4]T, b: [4]T) T {
    @setFloatMode(.strict);
    var sum: T = 0.0;
    inline for (0..4) |index| {
        sum = binaryResult(T, a, b, sum + a[index] * b[index]);
    }
    return sum;
}

/// The Euclidean length, folded through the accepted `hypot` kernel from a
/// positive zero. The fold cannot overflow while a single component cannot,
/// and it cannot flush a representable norm to zero, which
/// `sqrt(x*x + y*y + z*z + w*w)` does in both directions. Nesting the kernel is
/// not correctly rounded for four components; it is the accepted kernel
/// applied again, never a re-derived polynomial.
fn lengthOf(comptime T: type, a: [4]T) T {
    var accumulator: T = 0.0;
    inline for (0..4) |index| {
        accumulator = hypotOf(T, accumulator, a[index]);
    }
    return accumulator;
}

/// One step of a Hamilton component sum. Every step goes through the
/// special-value rule, so an invalid product cannot hide inside a valid sum.
/// Negation of a product is exact, so negating the product here and negating
/// one operand before multiplying give the same bits.
fn mulStep(comptime T: type, a: [4]T, b: [4]T, sum: T, product: T, negative: bool) T {
    return binaryResult(T, a, b, if (negative) sum - product else sum + product);
}

/// Hamilton multiplication, in the form `API.md:907` names:
///
///     w = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
///     x = a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y
///     y = a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x
///     z = a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w
///
/// Each line is summed strictly left to right and each product is formed
/// before it is added, so a caller that groups the terms differently gets a
/// different answer. The products do not need the special-value rule of their
/// own: a NaN product makes the sum NaN and the next step finds the NaN
/// operand, and an overflow to infinity cannot come back from an addition.
fn mulOf(comptime T: type, a: [4]T, b: [4]T) [4]T {
    @setFloatMode(.strict);

    var w: T = binaryResult(T, a, b, a[3] * b[3]);
    w = mulStep(T, a, b, w, a[0] * b[0], true);
    w = mulStep(T, a, b, w, a[1] * b[1], true);
    w = mulStep(T, a, b, w, a[2] * b[2], true);

    var x: T = binaryResult(T, a, b, a[3] * b[0]);
    x = mulStep(T, a, b, x, a[0] * b[3], false);
    x = mulStep(T, a, b, x, a[1] * b[2], false);
    x = mulStep(T, a, b, x, a[2] * b[1], true);

    var y: T = binaryResult(T, a, b, a[3] * b[1]);
    y = mulStep(T, a, b, y, a[0] * b[2], true);
    y = mulStep(T, a, b, y, a[1] * b[3], false);
    y = mulStep(T, a, b, y, a[2] * b[0], false);

    var z: T = binaryResult(T, a, b, a[3] * b[2]);
    z = mulStep(T, a, b, z, a[0] * b[1], false);
    z = mulStep(T, a, b, z, a[1] * b[0], true);
    z = mulStep(T, a, b, z, a[2] * b[3], false);

    return .{ x, y, z, w };
}

fn isFiniteOf(comptime T: type, a: [4]T) bool {
    var finite = true;
    inline for (0..4) |index| {
        finite = finite and std.math.isFinite(a[index]);
    }
    return finite;
}

/// Bit equality with NaN normalization: two NaNs match, everything else must
/// match bit for bit, so `+0` and `-0` do not match. Plain `==` cannot do this
/// job twice over: it is not reflexive on a NaN, and it cannot tell the two
/// zeros apart.
fn eqOf(comptime T: type, a: [4]T, b: [4]T) bool {
    var same = true;
    inline for (0..4) |index| {
        const left_is_nan = std.math.isNan(a[index]);
        const right_is_nan = std.math.isNan(b[index]);
        if (left_is_nan and right_is_nan) {
            // Any NaN equals any other NaN, whatever its sign or payload.
        } else if (left_is_nan != right_is_nan) {
            // A NaN is never equal to a number, which is what keeps the
            // relation an equivalence relation on the whole domain.
            return false;
        } else {
            same = same and @as(Bits(T), @bitCast(a[index])) == @as(Bits(T), @bitCast(b[index]));
        }
    }
    return same;
}

/// The band inside which a value counts as unit. Three `hypot` folds and one
/// division put a few ulps on the length, and `eps` is one ulp at one, so
/// eight of them clear that arithmetic with room to spare and is still far
/// tighter than any difference a caller could observe in a rotation.
fn unitToleranceOf(comptime T: type) T {
    return 8 * std.math.floatEps(T);
}

/// The value divided by its own length, or `null` when there is no direction to
/// normalize towards. A zero length has no quotient, and a length that
/// overflowed to infinity divides every finite component down to zero, which
/// is not a unit value either, so both are rejected rather than returned. The
/// division cannot overflow, because no component exceeds the length, and a
/// subnormal input needs no special handling: a subnormal component over a
/// subnormal length is an ordinary ratio near one with no exponent range left
/// to leave.
fn normalizedOf(comptime T: type, a: [4]T) ?[4]T {
    const length = lengthOf(T, a);
    if (length == 0.0 or !std.math.isFinite(length)) return null;
    return divByOf(T, a, length);
}

/// The value made unit, or the value unchanged when it cannot be made unit.
/// Every operation of the validated type re-normalizes its result, and for the
/// documented inputs the division always succeeds. The fallback is reachable
/// only through a struct literal, where there is nothing to correct towards
/// and handing back the arithmetic result is more honest than inventing a
/// direction.
fn renormalizedOf(comptime T: type, a: [4]T) [4]T {
    return normalizedOf(T, a) orelse a;
}

/// The `Vec3` of the same precision as `T`.
fn Vec3Of(comptime T: type) type {
    return switch (T) {
        f32 => Vec3f,
        f64 => Vec3d,
        else => @compileError("zmath quaternions are f32 or f64"),
    };
}

/// The unit direction of `axis`, or `null` when the axis has none. A component
/// that overflows the format has no representable direction left once the
/// length is taken, and comes back as zeros, so one zero-length test covers
/// both ways of having no axis.
fn unitAxisOf(comptime T: type, axis: Vec3Of(T)) ?Vec3Of(T) {
    if (!axis.isFinite()) return null;
    const direction = axis.normalize();
    if (direction.length() == 0.0) return null;
    return direction;
}

/// The dot product above which `slerp` stops following the great circle. Below
/// that angle `acos` has already lost its leading digits, and over the arc
/// that is left the normalized lerp is within a few ulps of the spherical one,
/// so the better-conditioned formula is the better answer. The threshold is
/// four square roots of `eps` below one, which for binary32 is 0.9986 and for
/// binary64 is 1 - 6.0e-8.
fn nearParallelOf(comptime T: type) T {
    const threshold: T = comptime 1.0 - 4.0 * @sqrt(std.math.floatEps(T));
    return threshold;
}

/// A dot product clamped into `[-1, 1]`. The dot product of two unit
/// quaternions is a cosine, and rounding can put it a hair outside, where
/// `acos` would answer a NaN. This is selection and not arithmetic, so a NaN
/// operand loses and yields the near end of the range.
fn clampedDotOf(comptime T: type, value: T) T {
    return maxOf(T, minOf(T, value, 1.0), -1.0);
}

/// Linear interpolation toward `to` for a `t` in `[0, 1]`. The endpoints are
/// short-circuited because the general form does not reach them: the
/// difference rounds before it is added back, so `a + (b - a) * 1` is not `b`.
fn lerpOf(comptime T: type, from: [4]T, to: [4]T, t: T) [4]T {
    if (t == 0.0) return from;
    if (t == 1.0) return to;
    return addOf(T, from, scaleOf(T, subOf(T, to, from), t));
}

// ---------------------------------------------------------------------------
// The raw algebra
// ---------------------------------------------------------------------------

/// A Hamilton quaternion of precision `T`. `w` is the scalar component.
///
/// The product is the Hamilton product, the product with a conjugate is a
/// square times the identity, and the unit ball around the identity is where
/// rotations live. This type carries no claim about being a rotation, so it
/// has no `rotate`: `q * (v, 0) * conjugate(q)` for a non-unit `q` scales `v`
/// by `normSquared(q)` and is a scaled rotation, not a rotation. Build a
/// `UnitQuat` when that is what the value means.
fn Quat(comptime T: type) type {
    return struct {
        const Self = @This();

        /// The first component of the vector part.
        x: T,
        /// The second component of the vector part.
        y: T,
        /// The third component of the vector part.
        z: T,
        /// The scalar component.
        w: T,

        /// Builds a value from named components.
        pub fn init(x: T, y: T, z: T, w: T) Self {
            return .{ .x = x, .y = y, .z = z, .w = w };
        }

        /// The value with `w` one and the rest zero, the multiplicative
        /// identity of the Hamilton product.
        pub fn identity() Self {
            return .{ .x = 0.0, .y = 0.0, .z = 0.0, .w = 1.0 };
        }

        /// Every component zero, the additive identity and the only quaternion
        /// with no inverse.
        pub fn zero() Self {
            return .{ .x = 0.0, .y = 0.0, .z = 0.0, .w = 0.0 };
        }

        /// The components as `[x, y, z, w]`.
        pub fn fromArray(a: [4]T) Self {
            return .{ .x = a[0], .y = a[1], .z = a[2], .w = a[3] };
        }

        /// The components as `[x, y, z, w]`.
        pub fn toArray(self: Self) [4]T {
            return toArrayOf(Self, self);
        }

        /// Componentwise sum.
        pub fn add(self: Self, other: Self) Self {
            return fromArrayOf(Self, addOf(T, toArrayOf(Self, self), toArrayOf(Self, other)));
        }

        /// Componentwise difference.
        pub fn sub(self: Self, other: Self) Self {
            return fromArrayOf(Self, subOf(T, toArrayOf(Self, self), toArrayOf(Self, other)));
        }

        /// Componentwise negation. The negated quaternion is `conjugate()`
        /// times `-1` on the scalar part alone, so it is not a conjugation.
        pub fn neg(self: Self) Self {
            return fromArrayOf(Self, negOf(T, toArrayOf(Self, self)));
        }

        /// Every component times `factor`.
        pub fn scale(self: Self, factor: T) Self {
            return fromArrayOf(Self, scaleOf(T, toArrayOf(Self, self), factor));
        }

        /// The componentwise product, which is not a quaternion product: `x`
        /// multiplies `x` and nothing rotates.
        pub fn scaleQuat(self: Self, other: Self) Self {
            return fromArrayOf(Self, scaleQuatOf(T, toArrayOf(Self, self), toArrayOf(Self, other)));
        }

        /// The scalar product, summed left to right in `x, y, z, w` order.
        pub fn dot(self: Self, other: Self) T {
            return dotOf(T, toArrayOf(Self, self), toArrayOf(Self, other));
        }

        /// The squared length, the dot product with itself. It overflows to
        /// infinity long before the length does, so it is a magnitude for
        /// moderate components and a test for a large one, not a scale factor.
        pub fn lengthSquared(self: Self) T {
            return dotOf(T, toArrayOf(Self, self), toArrayOf(Self, self));
        }

        /// The length, folded through the accepted `hypot` kernel, which
        /// cannot overflow while the components cannot and cannot flush a
        /// representable norm to zero.
        pub fn length(self: Self) T {
            return lengthOf(T, toArrayOf(Self, self));
        }

        /// The quaternion with the vector part negated and the scalar part left
        /// alone. Negation is exact, so the conjugate of a unit quaternion is
        /// unit to the bit.
        pub fn conjugate(self: Self) Self {
            const value = toArrayOf(Self, self);
            return fromArrayOf(Self, .{
                unaryResult(T, value, -value[0]),
                unaryResult(T, value, -value[1]),
                unaryResult(T, value, -value[2]),
                value[3],
            });
        }

        /// The squared norm, which is the real part of `self.mul(conjugate())`.
        /// It is that real part computed directly as `dot(self, self)`, whose
        /// terms all carry a positive sign, where the product sums them with
        /// three subtractions, so the two agree in exact arithmetic and to
        /// within a few ulps in floating point. That is why the round trip
        /// `self.mul(self.conjugate()) == normSquared() * identity()` is a
        /// tolerance test and not an equality test.
        pub fn normSquared(self: Self) T {
            return self.lengthSquared();
        }

        /// The multiplicative inverse, `conjugate() / normSquared()`, or `null`
        /// when there is none. A zero norm has no quotient and a norm that
        /// overflowed to infinity would scale the conjugate to zero, so both are
        /// rejected, as is a NaN norm.
        pub fn inverse(self: Self) ?Self {
            const norm = self.normSquared();
            if (norm == 0.0 or !std.math.isFinite(norm)) return null;
            return fromArrayOf(Self, divByOf(T, toArrayOf(Self, self.conjugate()), norm));
        }

        /// The Hamilton product `self * other`. The summation order is
        /// `API.md:907`'s, and it is not reassociated and not contracted into an
        /// FMA, so this is a function of the written formula and not of the
        /// optimizer. The product is not commutative, and a commuting pair is
        /// bit-identical in both orders only when the products in the two
        /// orderings are the same products, which is true of a pair about one
        /// axis and not of a general pair.
        pub fn mul(self: Self, other: Self) Self {
            return fromArrayOf(Self, mulOf(T, toArrayOf(Self, self), toArrayOf(Self, other)));
        }

        /// The quaternion for a rotation of `radians` about `axis`, with the
        /// axis normalized first, so a non-unit axis is accepted. This is an
        /// algebra-building constructor and not a rotation, which is why it can
        /// live on the raw type: the answer is a raw quaternion, and the caller
        /// decides whether it is a rotation by handing it to `UnitQuat.init`.
        /// A zero axis has no direction and gives `null`, as does a non-finite
        /// axis or a non-finite angle.
        pub fn fromAxisAngle(axis: Vec3Of(T), radians: T) ?Self {
            if (!std.math.isFinite(radians)) return null;
            const direction = unitAxisOf(T, axis) orelse return null;
            const half = radians * 0.5;
            const s = sinOf(T, half);
            return Self.init(direction.x * s, direction.y * s, direction.z * s, cosOf(T, half));
        }

        /// True when every component is finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteOf(T, toArrayOf(Self, self));
        }

        /// Bit equality, with every NaN equal to every other NaN. `+0` and `-0`
        /// are different components and are not equal, and every other pair
        /// must match bit for bit. That makes this a total equivalence relation,
        /// which is what a caller needs to check a result against a re-derived
        /// one.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(T, toArrayOf(Self, self), toArrayOf(Self, other));
        }
    };
}

// ---------------------------------------------------------------------------
// The validated rotation type
// ---------------------------------------------------------------------------

/// A quaternion of precision `T` that `init` has proved is unit, and therefore
/// stands for an active, right-handed rotation of a column vector.
///
/// `raw` is private to `zmath.quaternion` as the contract states, and Zig
/// cannot make it so. A struct field that a function of this module can write
/// is a field every consumer can write too, because `UnitQuatf{ .raw = q }`
/// compiles outside this file whatever this file says. What the language gives
/// instead is the constructor, and the constructor is the contract:
/// `UnitQuatf.init(q)` is the only way to build one of these that is
/// documented, it is the only one that validates, and it is the one a caller
/// should use. A value written with the struct literal can be any quaternion
/// at all, and every `rotate` result from such a value is meaningless rather
/// than merely inaccurate, because the rotation it performs scales by the norm
/// squared. `fromAxisAngle` is the second way in, and it routes through
/// `init`, so it cannot produce an unvalidated value either. The field stays
/// named `raw` because that is the name the contract gives it, and it is the
/// read path back to the algebra: read `q.raw` and hand it to a `Quat` when the
/// algebra is what is wanted.
fn Unit(comptime T: type) type {
    const Q = Quat(T);
    const Vec = Vec3Of(T);
    return struct {
        const Self = @This();

        /// The stored quaternion, which `init` made unit.
        raw: Q,

        /// The unit quaternion for `raw`, or `null` when there is none. A
        /// non-finite component, a zero quaternion, and a norm that overflowed
        /// to infinity all give `null`, and so does a value that survives
        /// normalization without landing inside the unit band, which is the
        /// specification's "not unit after normalization" case.
        ///
        /// A quaternion whose four components are all subnormal is one of
        /// those. Its length is subnormal and nonzero, so the division runs and
        /// cannot overflow, but the fold does not scale it: the accepted
        /// `hypot` kernel is correctly rounded only above the subnormal range,
        /// and for two subnormals the exact answer is not representable, so the
        /// fold reports the larger of the two. Four equal components therefore
        /// come back as four ones, whose length is not one, and the band
        /// rejects the value. Nothing on that path produces an infinity. A
        /// quaternion with a single subnormal component and three zeros does
        /// normalize, because folding zeros changes nothing.
        ///
        /// This is the only documented constructor, and it is the only one that
        /// validates; see the type comment.
        pub fn init(raw: Q) ?Self {
            const value = toArrayOf(Q, raw);
            if (!isFiniteOf(T, value)) return null;
            const unit = normalizedOf(T, value) orelse return null;
            if (@abs(lengthOf(T, unit) - 1.0) > unitToleranceOf(T)) return null;
            return .{ .raw = fromArrayOf(Q, unit) };
        }

        /// The scalar product, summed left to right in `x, y, z, w` order. For
        /// two unit quaternions it is the cosine of the angle between them,
        /// which is the quantity that decides whether an interpolation takes
        /// the short way round.
        pub fn dot(self: Self, other: Self) T {
            return dotOf(T, toArrayOf(Q, self.raw), toArrayOf(Q, other.raw));
        }

        /// Bit equality, with every NaN equal to every other NaN. It compares
        /// the stored algebra, so two unit quaternions that are the same
        /// rotation by different routes do not match unless they match to the
        /// bit.
        pub fn eq(self: Self, other: Self) bool {
            return self.raw.eq(other.raw);
        }

        /// The conjugate, which for a unit quaternion is its own inverse. The
        /// three negations are exact, so this result is unit to the bit and is
        /// not re-normalized.
        pub fn conjugate(self: Self) Self {
            return .{ .raw = self.raw.conjugate() };
        }

        /// The inverse, which for a unit quaternion is the conjugate. It is not
        /// a division: `conj(q) / normSquared(q)` and `conj(q)` are the same
        /// number, because the norm squared is one, and the division would only
        /// add a rounding the multiplication does not need. A caller who has a
        /// `UnitQuat` has a proved norm, so the identity is available without
        /// re-deriving it. Build a `Quat` and ask it for `inverse` when the norm
        /// is not known to be one.
        pub fn inverse(self: Self) Self {
            return self.conjugate();
        }

        /// The Hamilton product, re-normalized.
        ///
        /// The product of two unit quaternions is unit in exact arithmetic, so
        /// the division below corrects the few ulps the two products rounded
        /// away, and it is applied rather than assumed because a quaternion
        /// that drifts from the unit ball drifts out of it faster than a
        /// pipeline of rotations can measure. A caller that can only reach
        /// values through `init` never sees the fallback: it is reachable only
        /// through the struct literal, where there is nothing to correct
        /// towards.
        pub fn mul(self: Self, b: Self) Self {
            const product = self.raw.mul(b.raw);
            return .{ .raw = fromArrayOf(Q, renormalizedOf(T, toArrayOf(Q, product))) };
        }

        /// Apply `self` first and then `b`, which is `b.mul(self)`.
        ///
        /// The reversal is the point of the method. Composition of functions
        /// reads left to right and composition of quaternions reads right to
        /// left, so `a.then(b)` is the only spelling in which a rotation
        /// pipeline reads in the order it runs, and
        /// `a.then(b).rotate(v)` is `b.rotate(a.rotate(v))` to rounding.
        pub fn then(self: Self, b: Self) Self {
            return b.mul(self);
        }

        /// The vector `q * (v, 0) * conjugate(q)`, computed through the raw
        /// Hamilton product in that association so the answer is the one the
        /// contract names and not a rearranged matrix form. The result has the
        /// same length as `v`: a unit quaternion is a rotation, and the
        /// conjugate is the inverse, so the two products cancel the scale.
        pub fn rotate(self: Self, v: Vec) Vec {
            const q = self.raw;
            const moved = q.mul(Q.init(v.x, v.y, v.z, 0.0)).mul(q.conjugate());
            return Vec.init(moved.x, moved.y, moved.z);
        }

        /// The unit quaternion for a rotation of `radians` about `axis`, with
        /// the axis normalized first, so a non-unit axis is accepted. A zero
        /// axis has no direction and gives `null`, as does a non-finite axis or
        /// a non-finite angle. The raw quaternion is built by
        /// `Quat.fromAxisAngle` and validated by `init`, so this cannot produce
        /// a value that `init` would have rejected.
        pub fn fromAxisAngle(axis: Vec, radians: T) ?Self {
            return init(Q.fromAxisAngle(axis, radians) orelse return null);
        }

        /// The spherical interpolation toward `b` for a `t` in `[0, 1]`, along
        /// the short arc, or `null` for a non-finite `t`.
        ///
        /// `b` and `-b` name the same rotation, so which arc between them is
        /// the short one is decided by the sign of the dot product and by
        /// nothing else. A negative dot product means the arc toward `b` is the
        /// long one, so the endpoint is negated to `-b` and the dot product
        /// becomes positive with it; negating the endpoint unconditionally
        /// would send every pair the long way round. Below `nearParallel` the
        /// great circle is followed exactly; above it `acos` has lost the
        /// relative accuracy the arc needs and the normalized lerp takes over,
        /// which is the better-conditioned answer over the arc that remains.
        ///
        /// The result is the lift of the arc that starts at `self`, so its dot
        /// product with `self` is positive and its dot product with `b` is
        /// negative when the arc was flipped. Callers that need a
        /// representation-independent comparison have to use the magnitude, or
        /// compare the rotations rather than the quaternions.
        pub fn slerp(self: Self, b: Self, t: T) ?Self {
            if (!std.math.isFinite(t)) return null;
            if (t == 0.0) return self;
            if (t == 1.0) return b;

            const from = toArrayOf(Q, self.raw);
            var to = toArrayOf(Q, b.raw);
            var d = dotOf(T, from, to);
            if (d < 0.0) {
                to = negOf(T, to);
                d = -d;
            }

            if (d > nearParallelOf(T)) {
                return .{ .raw = fromArrayOf(Q, renormalizedOf(T, lerpOf(T, from, to, t))) };
            }

            const theta = acosOf(T, clampedDotOf(T, d));
            const s = sinOf(T, theta);
            if (s == 0.0) {
                // Only a value built with the struct literal can make the two
                // unit directions identical, and then the lerp is the defined
                // answer rather than a division by a zero sine.
                return .{ .raw = fromArrayOf(Q, renormalizedOf(T, lerpOf(T, from, to, t))) };
            }
            const left = sinOf(T, (1.0 - t) * theta) / s;
            const right = sinOf(T, t * theta) / s;
            const mixed = addOf(T, scaleOf(T, from, left), scaleOf(T, to, right));
            return init(fromArrayOf(Q, mixed));
        }

        /// The normalized linear interpolation toward `b` for a `t` in `[0, 1]`.
        ///
        /// The chord is re-normalized, so the result is on the unit sphere, but
        /// it is not on the great circle: the angle advances in `sin` of the
        /// interpolated angle rather than in the angle itself, so `lerp` is for
        /// a blend and `slerp` is for a path. A non-finite `t` has no meaning
        /// here and no optional to report it in, so it propagates as a NaN
        /// quaternion from the component arithmetic; `slerp` is the entry point
        /// that reports it.
        pub fn lerp(self: Self, b: Self, t: T) Self {
            const mixed = lerpOf(T, toArrayOf(Q, self.raw), toArrayOf(Q, b.raw), t);
            return .{ .raw = fromArrayOf(Q, renormalizedOf(T, mixed)) };
        }
    };
}

/// A raw `f32` quaternion, the algebra and not a rotation.
pub const Quatf = Quat(f32);

/// A raw `f64` quaternion, the algebra and not a rotation.
pub const Quatd = Quat(f64);

/// A validated `f32` rotation, an active right-handed rotation of a column
/// vector.
pub const UnitQuatf = Unit(f32);

/// A validated `f64` rotation, an active right-handed rotation of a column
/// vector.
pub const UnitQuatd = Unit(f64);
