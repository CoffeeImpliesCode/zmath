//! Concrete `f32`/`f64` families in two, three, and four dimensions.
//!
//! `API.md:973-991` is the whole contract. This file spells it out.
//!
//! ## Naming
//!
//! Every public type carries its dimension and its precision, in that order:
//! `Point2f`, `Vector3d`, `Line4f`, `Aabb2d`, `Hypersphere4f`. The `f`/`d`
//! suffix is the project's convention (`scalar.hypotf`, `scalar.hypotd`), and
//! the digit is the number of components. A two-component `f32` point is
//! `Point2f`; the same object in binary64 is `Point2d`. Nothing in the public
//! surface is generic.
//!
//! ## Finiteness
//!
//! Every operation that CONSUMES a float, a point, a vector, or a scalar
//! validates it first and returns `error.NonFiniteInput`. A geometry call in
//! this module therefore never propagates an infinity or a NaN, and never
//! returns a NaN where an error belongs. The operations that take no float
//! argument, which are the pure queries `isFinite` and the `eq` of a point or
//! a vector, are infallible by construction.
//!
//! Each constructor adds the one failure mode its own contract names:
//!
//! - `Line`, `Ray`, `Plane`, `Hyperplane`, `Vector.normalize`:
//!   `error{ NonFiniteInput, Degenerate }`. The direction or normal is
//!   normalized, and a zero vector has no direction.
//! - `Circle`, `Sphere`, `Hypersphere`: `error{ NonFiniteInput, NegativeRadius }`.
//! - `Aabb`: `error{ NonFiniteInput, InvertedBounds }`.
//!
//! ## A point and a vector are different types
//!
//! `Point2f` and `Vector2f` are distinct structs. `Point2f.add` takes a
//! `Vector2f`; `Point2f.toVector` converts. `Line2f.init` takes a point and a
//! vector because a point and a direction are what a line is made of, and
//! neither type is accepted where the other is named.
//!
//! ## Predicate policy
//!
//! Every predicate that can DIFFER between policies takes the policy as an
//! argument, and nothing picks one for the caller. `fast` reports the sign of
//! the plain scalar evaluation, in the operation order this module documents,
//! with contraction disabled. `adaptive_exact` reports the exact sign for
//! finite `f32` and `f64` inputs. There is no epsilon anywhere:
//! `adaptive_exact` is exact, not tighter, and it never falls back to the
//! `fast` answer.
//!
//! Predicates that CANNOT differ take no policy, because a policy argument a
//! predicate ignores is a lie. An `Aabb` contains an axis-aligned range
//! comparison, which rounds nothing, so `Aabb.contains`, `Aabb.intersects`,
//! and `Aabb.intersect` take no policy. Every polynomial predicate -
//! orientation, sphere containment, sphere to sphere, sphere to flat, and
//! every affine intersection that reduces to them - takes one.
//!
//! `adaptive_exact` is exact for the values handed to it. The three- and
//! four-dimensional affine intersections first decide coplanarity with an
//! exact predicate on the ORIGINAL coordinates, then project into a plane and
//! run the exact two-dimensional predicate on the PROJECTED coordinates. The
//! coplanarity gate, which is the decision a plain computation gets wrong, is
//! therefore exact on the inputs as given; the in-plane decisions are exact
//! for the frame the projection produced. The intersection PARAMETERS are
//! `f32`/`f64` values, and the range tests on them are ordinary comparisons.
//!
//! ## Exact arithmetic
//!
//! Every finite `f64` is `m * 2^e` with `|m| < 2^53` and
//! `-1074 <= e <= 971`, so `x * 2^1074` is an integer below `2^2098`. An `f32`
//! widens to `f64` exactly, so one representation covers both precisions.
//! Every term of a predicate has the same degree, so every term of its integer
//! value carries the same power of two, and that power is positive. The sign
//! of the integer is the sign of the real value, and nothing rounded.
//!
//! Storage is a fixed `[132]u64` magnitude plus a sign, on the stack. 132
//! limbs is 8448 bits. The worst predicate here is degree four: a four-by-four
//! determinant has 24 terms, each a product of four scaled differences below
//! `2^2099`, so its exact value is below `24 * 2^8396 < 2^8401`. The
//! degree-three and degree-two forms fit in 66 limbs and the exact path
//! allocates the same array for them anyway, so the storage is one size for
//! the whole module. Nothing allocates, and no limb count depends on a value.
//!
//! ## Intersection unions
//!
//! One union per dimension and precision, `Intersection2f` through
//! `Intersection4d`. Each is a tagged union; none of them uses `null` to
//! combine cases, and `none` is an ordinary tag like any other. The tags name
//! the DIMENSION of the intersection's affine span - zero `point`, one `line`,
//! two `plane`, three `hyperplane` - with the bounded and solid forms named
//! separately, plus three tags the contract's nine do not cover:
//!
//! - `box`: an AABB. The intersection of two boxes is a box, and none of the
//!   nine names it.
//! - `circle` in three dimensions and `hypersphere` in four: a flat that cuts
//!   a ball gives back a ball, and that is a curve, not a finite set of
//!   points.
//!
//! `coincident` means the two operands are the same set and that set is not an
//! infinite flat, so it covers identical segments, spheres, and boxes.
//! `contained` means the intersection is exactly one operand and the operands
//! are not the same set. A coincident pair of lines, planes, or hyperplanes
//! returns the tag that names its own span instead, because "the intersection
//! is this line" is the more useful answer than "the two are the same".
//!
//! Every tag is the WHOLE common set, never a part of it, and that is what
//! settles the two cases where a part would have been convenient:
//!
//! - A ray whose origin is inside a ball meets it in the whole run from the
//!   origin to the exit, so the tag is `segment`, not the exit alone. The
//!   origin is a common point, so reporting only the exit would report a set
//!   that is not the intersection. The run is a `point` when the origin is
//!   already on the surface and the ray points away, because then it has
//!   zero length. A ray whose origin is OUTSIDE reports `two_points`, or
//!   `point` when the run is tangent and the two ends coincide.
//! - A segment that merely ENDS on a ball or a flat is not a touch. A
//!   segment that reaches the surface from outside crosses into the solid and
//!   its common set is a proper sub-segment whose far end is that surface
//!   point. A true touch is a DOUBLE root, and only then is the common set
//!   the single tangent point.
//!
//! A ball is the CLOSED solid ball, so its boundary belongs to it and a
//! segment ending exactly on the surface is inside up to that point.
//!
//! ## Deliberate omissions
//!
//! - `Vector.cross` exists in three dimensions only. In four dimensions there
//!   is no cross product, because the orthogonal complement of the plane two
//!   vectors span is a whole plane rather than a line. `Vector4.wedge`
//!   provides the exterior product `u ^ v`, which is the product that
//!   dimension actually has.
//! - `Point2.side` exists in two dimensions only. A line divides the plane
//!   into two sides; a line in three or four dimensions does not.
//! - `Ray` has no `closestPointTo`. The contract names `at` for a ray and
//!   `closestPointTo` for a segment, and that asymmetry is kept rather than
//!   smoothed over.
//!
//! `ContainedSide` names the operand the intersection EQUALS, and it names it
//! by its position in the call, so the receiver of a method is `first` and
//! its argument is `second`. A box that holds its argument has the ARGUMENT
//! as its whole common set, so it reports `contained = .second`, exactly as
//! a segment lying inside the ball or the flat it is tested against does. It
//! is `.first` only when the receiver is the one wholly inside, which for
//! `unit.intersect(inner)` means the reversed call `inner.intersect(unit)`.
const std = @import("std");
const scalar = @import("scalar.zig");

// ---------------------------------------------------------------------------
// Policy, sign, and the shared scalar kernels
// ---------------------------------------------------------------------------

/// Which sign a predicate reports.
///
/// `fast` is the sign of the plain scalar evaluation. `adaptive_exact` is the
/// sign of the exact value of the same polynomial. The two differ only where
/// rounding decides the sign, and the caller chooses which one it wants.
pub const PredicatePolicy = enum { fast, adaptive_exact };

/// The sign of a real value.
pub const Sign = enum { negative, zero, positive };

/// Which operand an intersection equals, when the intersection is one whole
/// operand and the two operands are not the same set.
pub const ContainedSide = enum { first, second };

fn hypotOf(comptime T: type, a: T, b: T) T {
    return switch (T) {
        f32 => scalar.hypotf(a, b),
        f64 => scalar.hypotd(a, b),
        else => @compileError("zmath geometry is f32 or f64"),
    };
}

fn sqrtOf(comptime T: type, a: T) T {
    return switch (T) {
        f32 => scalar.sqrtf(a),
        f64 => scalar.sqrtd(a),
        else => @compileError("zmath geometry is f32 or f64"),
    };
}

fn minOf(comptime T: type, a: T, b: T) T {
    return switch (T) {
        f32 => scalar.fminf(a, b),
        f64 => scalar.fmind(a, b),
        else => @compileError("zmath geometry is f32 or f64"),
    };
}

fn maxOf(comptime T: type, a: T, b: T) T {
    return switch (T) {
        f32 => scalar.fmaxf(a, b),
        f64 => scalar.fmaxd(a, b),
        else => @compileError("zmath geometry is f32 or f64"),
    };
}

fn absOf(comptime T: type, x: T) T {
    return switch (T) {
        f32 => scalar.absf(x),
        f64 => scalar.absd(x),
        else => @compileError("zmath geometry is f32 or f64"),
    };
}

fn isFiniteOf(comptime T: type, x: T) bool {
    return switch (T) {
        f32 => std.math.isFinite(x),
        f64 => std.math.isFinite(x),
        else => @compileError("zmath geometry is f32 or f64"),
    };
}

fn signOfScalar(comptime T: type, x: T) Sign {
    if (x > 0) return .positive;
    if (x < 0) return .negative;
    return .zero;
}

/// `x` widened to `f64`, which is exact for every finite `f32`.
fn widen(comptime T: type, x: T) f64 {
    return @floatCast(x);
}

fn clampRange(comptime T: type, value: T, lower: T, upper: T) T {
    if (value < lower) return lower;
    if (value > upper) return upper;
    return value;
}

fn clampUnitOf(comptime T: type, t: T) T {
    return clampRange(T, t, 0, 1);
}

/// The open end of the parameter range of a line or a ray. It is a range
/// bound only, never an input: nothing here validates it, because the
/// comparisons that use it accept an infinite bound by construction.
fn infinityOf(comptime T: type) T {
    return std.math.inf(T);
}

// ---------------------------------------------------------------------------
// Field-name plumbing
//
// Every value type is a struct of named components, so the shared algorithms
// walk the fields by name. The names and the order are the ones the value
// layer in `vector_value.zig` already uses, so the two layers cannot drift.
// ---------------------------------------------------------------------------

const component_names = [4][]const u8{ "x", "y", "z", "w" };

fn fieldCount(comptime V: type) usize {
    return @typeInfo(V).@"struct".fields.len;
}

fn fieldName(comptime index: usize) []const u8 {
    return component_names[index];
}

fn componentType(comptime V: type) type {
    return @typeInfo(V).@"struct".fields[0].type;
}

fn toArrayOf(comptime V: type, a: V) [fieldCount(V)]componentType(V) {
    const width = comptime fieldCount(V);
    var out: [width]componentType(V) = undefined;
    inline for (0..width) |index| {
        out[index] = @field(a, fieldName(index));
    }
    return out;
}

fn fromArrayOf(comptime V: type, a: [fieldCount(V)]componentType(V)) V {
    const width = comptime fieldCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        @field(out, fieldName(index)) = a[index];
    }
    return out;
}

/// The finiteness gate. It walks the value's OWN fields, so it covers a point,
/// a vector, and every composite - a line, a ray, a segment, a ball, a flat,
/// a box - without the caller naming the field's type at the call site. A
/// composite's fields are not all named after an axis, so the axis-named walk
/// that the point and vector types use cannot be shared with them.
fn requireFiniteField(value: anytype) error{NonFiniteInput}!void {
    const T = @TypeOf(value);
    if (T == f32 or T == f64) {
        if (!std.math.isFinite(value)) return error.NonFiniteInput;
        return;
    }
    const width = comptime @typeInfo(T).@"struct".fields.len;
    inline for (0..width) |index| {
        const name = comptime @typeInfo(T).@"struct".fields[index].name;
        try requireFiniteField(@field(value, name));
    }
}

fn requireFiniteScalar(comptime T: type, a: T) error{NonFiniteInput}!void {
    if (!isFiniteOf(T, a)) return error.NonFiniteInput;
}

fn isFiniteQuiet(comptime V: type, a: V) bool {
    const width = comptime fieldCount(V);
    var finite = true;
    inline for (0..width) |index| {
        finite = finite and isFiniteOf(componentType(V), @field(a, fieldName(index)));
    }
    return finite;
}

/// Bit equality with the NaN normalization the value layer uses: two NaNs
/// match, `+0` and `-0` do not.
fn eqOf(comptime V: type, a: V, b: V) bool {
    const T = componentType(V);
    const width = comptime fieldCount(V);
    var same = true;
    inline for (0..width) |index| {
        const left = @field(a, fieldName(index));
        const right = @field(b, fieldName(index));
        const left_nan = std.math.isNan(left);
        if (left_nan != std.math.isNan(right)) {
            same = false;
        } else if (left_nan) {
            same = same;
        } else {
            same = same and switch (T) {
                f32 => @as(u32, @bitCast(left)) == @as(u32, @bitCast(right)),
                f64 => @as(u64, @bitCast(left)) == @as(u64, @bitCast(right)),
                else => @compileError("zmath geometry is f32 or f64"),
            };
        }
    }
    return same;
}

// ---------------------------------------------------------------------------
// Exact sign arithmetic
// ---------------------------------------------------------------------------

/// Limbs of exact magnitude. The module comment proves 132 is enough for the
/// degree-four worst case.
const exact_limbs = 132;

/// An exact signed integer, used only to decide a sign. A degree-`k`
/// polynomial's value here is the real value times `2^(1074 * k)`, a positive
/// factor, so the sign of the integer is the sign of the real value.
const Big = struct {
    negative: bool,
    used: u8,
    limbs: [exact_limbs]u64,

    const zero_value: Big = .{
        .negative = false,
        .used = 0,
        .limbs = [_]u64{0} ** exact_limbs,
    };

    fn zero() Big {
        return zero_value;
    }

    fn sign(a: Big) Sign {
        if (a.used == 0) return .zero;
        return if (a.negative) .negative else .positive;
    }

    /// The exact integer `value * 2^1074`. The caller has already rejected a
    /// non-finite input, which is the only thing that reaches here.
    fn fromFloat(value: f64) Big {
        const bits: u64 = @bitCast(value);
        const exponent_field: u64 = (bits >> 52) & 0x7ff;
        const fraction = bits & 0x000f_ffff_ffff_ffff;
        if (exponent_field == 0x7ff) return zero_value;
        if (exponent_field == 0 and fraction == 0) return zero_value;
        // A normal value is `significand * 2^(E - 1075)` and a subnormal is
        // `fraction * 2^-1074`; scaling by 2^1074 makes the shift a
        // non-negative limb offset.
        const significand = fraction | (if (exponent_field == 0) @as(u64, 0) else (@as(u64, 1) << 52));
        const offset: u32 = if (exponent_field == 0) 0 else @intCast(exponent_field - 1);
        const word: usize = offset / 64;
        const shift: u6 = @intCast(offset % 64);
        var limbs = [_]u64{0} ** exact_limbs;
        // A 128-bit shift writes the two limbs at once, which avoids a
        // shift-by-64 that the six-bit shift count cannot even name.
        const placed: u128 = @as(u128, significand) << shift;
        limbs[word] = @truncate(placed);
        limbs[word + 1] = @truncate(placed >> 64);
        return .{
            .negative = (bits >> 63) != 0,
            .used = topUsed(&limbs),
            .limbs = limbs,
        };
    }

    fn neg(a: Big) Big {
        if (a.used == 0) return a;
        var out = a;
        out.negative = !a.negative;
        return out;
    }

    fn add(a: Big, b: Big) Big {
        if (a.used == 0) return b;
        if (b.used == 0) return a;
        if (a.negative == b.negative) return magnitudeSum(a, b, a.negative);
        const order = magnitudeCompare(a, b);
        if (order == .eq) return zero_value;
        if (order == .gt) return magnitudeDifference(a, b, a.negative);
        return magnitudeDifference(b, a, b.negative);
    }

    fn sub(a: Big, b: Big) Big {
        return a.add(b.neg());
    }

    /// Schoolbook, truncated at the array. The bounds in the module comment
    /// prove the product fits, so the truncation is a bound check and not a
    /// loss.
    fn mul(a: Big, b: Big) Big {
        if (a.used == 0 or b.used == 0) return zero_value;
        var limbs = [_]u64{0} ** exact_limbs;
        var i: usize = 0;
        while (i < a.used and i < exact_limbs) : (i += 1) {
            const row = @min(b.used, exact_limbs - i);
            var carry: u128 = 0;
            var j: usize = 0;
            while (j < row) : (j += 1) {
                const total: u128 = @as(u128, a.limbs[i]) * @as(u128, b.limbs[j]) +
                    @as(u128, limbs[i + j]) + carry;
                limbs[i + j] = @truncate(total);
                carry = total >> 64;
            }
            if (i + row < exact_limbs) {
                const total: u128 = @as(u128, limbs[i + row]) + carry;
                limbs[i + row] = @truncate(total);
            }
        }
        return .{
            .negative = a.negative != b.negative,
            .used = topUsed(&limbs),
            .limbs = limbs,
        };
    }
};

fn topUsed(limbs: *const [exact_limbs]u64) u8 {
    var index: usize = exact_limbs;
    while (index > 0) {
        index -= 1;
        if (limbs[index] != 0) return @intCast(index + 1);
    }
    return 0;
}

fn magnitudeCompare(a: Big, b: Big) std.math.Order {
    if (a.used != b.used) return if (a.used > b.used) .gt else .lt;
    var index: usize = a.used;
    while (index > 0) {
        index -= 1;
        const left = a.limbs[index];
        const right = b.limbs[index];
        if (left != right) return if (left > right) .gt else .lt;
    }
    return .eq;
}

fn magnitudeSum(a: Big, b: Big, negative: bool) Big {
    var limbs = [_]u64{0} ** exact_limbs;
    const width = @max(a.used, b.used);
    var carry: u64 = 0;
    var index: usize = 0;
    while (index < width) : (index += 1) {
        const left = if (index < a.used) a.limbs[index] else 0;
        const right = if (index < b.used) b.limbs[index] else 0;
        const total: u128 = @as(u128, left) + @as(u128, right) + @as(u128, carry);
        limbs[index] = @truncate(total);
        carry = @intCast(total >> 64);
    }
    if (carry != 0 and width < exact_limbs) limbs[width] = carry;
    return .{ .negative = negative, .used = topUsed(&limbs), .limbs = limbs };
}

/// `a - b` in magnitude. The caller has already ordered the operands, so the
/// borrow always resolves inside the loop. The arithmetic is WRAPPING, not
/// checked: an intermediate `left - right - borrow` goes below zero exactly
/// when the incoming borrow exceeds the two limbs, and a checked unsigned
/// subtraction would turn that ordinary borrow into a panic.
fn magnitudeDifference(a: Big, b: Big, negative: bool) Big {
    var limbs = [_]u64{0} ** exact_limbs;
    var borrow: u64 = 0;
    var index: usize = 0;
    while (index < a.used) : (index += 1) {
        const left = a.limbs[index];
        const right = if (index < b.used) b.limbs[index] else 0;
        const total: u128 = @as(u128, left) -% @as(u128, right) -% @as(u128, borrow);
        limbs[index] = @truncate(total);
        borrow = @intCast((~(total >> 64)) & 1);
    }
    return .{ .negative = negative, .used = topUsed(&limbs), .limbs = limbs };
}

/// The exact sign of `a + b` over `dim` components, degree one.
fn exactSumOfProducts(comptime dim: usize, values: [dim]f64, factors: [dim]f64) Big {
    var total = Big.zero();
    inline for (0..dim) |index| {
        total = total.add(Big.fromFloat(values[index]).mul(Big.fromFloat(factors[index])));
    }
    return total;
}

/// The exact sign of `dot(x, x) - dot(y, y)`, degree two.
fn exactSquaredDifference(comptime dim: usize, x: [dim]f64, y: [dim]f64) Sign {
    const left = exactSumOfProducts(dim, x, x);
    const right = exactSumOfProducts(dim, y, y);
    return left.sub(right).sign();
}

/// The exact sign of `p * q - r * s`, degree two.
fn exactDet2(p: f64, q: f64, r: f64, s: f64) Sign {
    return Big.fromFloat(p).mul(Big.fromFloat(q))
        .sub(Big.fromFloat(r).mul(Big.fromFloat(s)))
        .sign();
}

/// The exact sign of the two-component oriented area `(b - a) x (c - a)`.
fn exactOrient2(comptime T: type, a: [2]T, b: [2]T, c: [2]T) Sign {
    return exactDet2(
        widen(T, b[0]) - widen(T, a[0]),
        widen(T, c[1]) - widen(T, a[1]),
        widen(T, b[1]) - widen(T, a[1]),
        widen(T, c[0]) - widen(T, a[0]),
    );
}

/// The exact sign of `dot(cross(b - a, c - a), d - a)`, degree three. The
/// expression matches the fast path one level down, so the two policies differ
/// only in the arithmetic and never in the algebra.
fn exactOrient3(comptime T: type, a: [3]T, b: [3]T, c: [3]T, d: [3]T) Sign {
    const u = [3]f64{ widen(T, b[0]) - widen(T, a[0]), widen(T, b[1]) - widen(T, a[1]), widen(T, b[2]) - widen(T, a[2]) };
    const v = [3]f64{ widen(T, c[0]) - widen(T, a[0]), widen(T, c[1]) - widen(T, a[1]), widen(T, c[2]) - widen(T, a[2]) };
    const w = [3]f64{ widen(T, d[0]) - widen(T, a[0]), widen(T, d[1]) - widen(T, a[1]), widen(T, d[2]) - widen(T, a[2]) };
    const nx = Big.fromFloat(u[1]).mul(Big.fromFloat(v[2])).sub(Big.fromFloat(u[2]).mul(Big.fromFloat(v[1])));
    const ny = Big.fromFloat(u[2]).mul(Big.fromFloat(v[0])).sub(Big.fromFloat(u[0]).mul(Big.fromFloat(v[2])));
    const nz = Big.fromFloat(u[0]).mul(Big.fromFloat(v[1])).sub(Big.fromFloat(u[1]).mul(Big.fromFloat(v[0])));
    var total = nx.mul(Big.fromFloat(w[0]));
    total = total.add(ny.mul(Big.fromFloat(w[1])));
    total = total.add(nz.mul(Big.fromFloat(w[2])));
    return total.sign();
}

/// The exact three-by-three minor of a four-by-four matrix: the determinant
/// of the three selected rows against the three selected columns, expanded
/// along the third column into three terms of two products each.
fn exactMinor3(rows: [4][4]f64, r0: usize, r1: usize, r2: usize, c0: usize, c1: usize, c2: usize) Big {
    var total = cell(rows, r0, c2).mul(cell(rows, r1, c0).mul(cell(rows, r2, c1)).sub(cell(rows, r1, c1).mul(cell(rows, r2, c0))));
    total = total.sub(cell(rows, r1, c2).mul(cell(rows, r0, c0).mul(cell(rows, r2, c1)).sub(cell(rows, r0, c1).mul(cell(rows, r2, c0)))));
    total = total.add(cell(rows, r2, c2).mul(cell(rows, r0, c0).mul(cell(rows, r1, c1)).sub(cell(rows, r0, c1).mul(cell(rows, r1, c0)))));
    return total;
}

fn cell(rows: [4][4]f64, row: usize, column: usize) Big {
    return Big.fromFloat(rows[row][column]);
}

/// The exact sign of the four-by-four determinant whose rows are the four
/// points, degree four. The expansion runs along the first row, so the
/// cofactor signs alternate and each minor skips the column its own
/// coefficient multiplies.
fn exactOrient4(comptime T: type, a: [4]T, b: [4]T, c: [4]T, d: [4]T) Sign {
    const rows: [4][4]f64 = .{ widenRow(T, a), widenRow(T, b), widenRow(T, c), widenRow(T, d) };
    var total = exactMinor3(rows, 1, 2, 3, 1, 2, 3).mul(Big.fromFloat(rows[0][0]));
    total = total.sub(exactMinor3(rows, 1, 2, 3, 0, 2, 3).mul(Big.fromFloat(rows[0][1])));
    total = total.add(exactMinor3(rows, 1, 2, 3, 0, 1, 3).mul(Big.fromFloat(rows[0][2])));
    total = total.sub(exactMinor3(rows, 1, 2, 3, 0, 1, 2).mul(Big.fromFloat(rows[0][3])));
    return total.sign();
}

fn widenRow(comptime T: type, a: [4]T) [4]f64 {
    return .{ widen(T, a[0]), widen(T, a[1]), widen(T, a[2]), widen(T, a[3]) };
}

/// The plain evaluation of the same four-by-four determinant, with the same
/// expansion along the first row and the same grouping inside each minor, so
/// that the two policies differ only in the arithmetic and never in the
/// algebra.
fn fastOrient4(comptime T: type, a: [4]T, b: [4]T, c: [4]T, d: [4]T) T {
    @setFloatMode(.strict);
    const rows = [4][4]T{ a, b, c, d };
    return minor3FastOf(T, rows, 1, 2, 3, 1, 2, 3) * rows[0][0] -
        minor3FastOf(T, rows, 1, 2, 3, 0, 2, 3) * rows[0][1] +
        minor3FastOf(T, rows, 1, 2, 3, 0, 1, 3) * rows[0][2] -
        minor3FastOf(T, rows, 1, 2, 3, 0, 1, 2) * rows[0][3];
}

fn minor3FastOf(comptime T: type, rows: [4][4]T, r0: usize, r1: usize, r2: usize, c0: usize, c1: usize, c2: usize) T {
    @setFloatMode(.strict);
    return rows[r0][c2] * (rows[r1][c0] * rows[r2][c1] - rows[r1][c1] * rows[r2][c0]) -
        rows[r1][c2] * (rows[r0][c0] * rows[r2][c1] - rows[r0][c1] * rows[r2][c0]) +
        rows[r2][c2] * (rows[r0][c0] * rows[r1][c1] - rows[r0][c1] * rows[r1][c0]);
}

fn orient4Sign(comptime T: type, a: [4]T, b: [4]T, c: [4]T, d: [4]T, policy: PredicatePolicy) Sign {
    return switch (policy) {
        .fast => signOfScalar(T, fastOrient4(T, a, b, c, d)),
        .adaptive_exact => exactOrient4(T, a, b, c, d),
    };
}

/// The exact sign of the three-by-three minor of the three-by-four matrix
/// whose columns are `u`, `v`, and `w`, over the three selected rows. Four of
/// these vanish exactly when the three vectors are linearly dependent, which
/// is the four-dimensional coplanarity test.
fn exactMinor3x4(comptime T: type, u: [4]T, v: [4]T, w: [4]T, rows: [3]usize) Sign {
    const columns = [3][4]f64{ widenRow(T, u), widenRow(T, v), widenRow(T, w) };
    const rows_d: [3][3]f64 = .{
        .{ columns[0][rows[0]], columns[1][rows[0]], columns[2][rows[0]] },
        .{ columns[0][rows[1]], columns[1][rows[1]], columns[2][rows[1]] },
        .{ columns[0][rows[2]], columns[1][rows[2]], columns[2][rows[2]] },
    };
    return minor3Plain(rows_d).sign();
}

/// The three-by-three determinant, expanded along the THIRD column into three
/// terms of two products each. The third column is the one the caller selected
/// rows against, so this is the cofactor expansion the four-by-four
/// determinant needs, and the exact and fast paths below use the same grouping.
fn minor3Plain(rows: [3][3]f64) Big {
    var total = Big.fromFloat(rows[0][2]).mul(
        Big.fromFloat(rows[1][0]).mul(Big.fromFloat(rows[2][1])).sub(Big.fromFloat(rows[1][1]).mul(Big.fromFloat(rows[2][0]))),
    );
    total = total.add(Big.fromFloat(rows[1][2]).mul(
        Big.fromFloat(rows[0][1]).mul(Big.fromFloat(rows[2][0])).sub(Big.fromFloat(rows[0][0]).mul(Big.fromFloat(rows[2][1]))),
    ));
    total = total.add(Big.fromFloat(rows[2][2]).mul(
        Big.fromFloat(rows[0][0]).mul(Big.fromFloat(rows[1][1])).sub(Big.fromFloat(rows[0][1]).mul(Big.fromFloat(rows[1][0]))),
    ));
    return total;
}

fn fastMinor3x4(comptime T: type, u: [4]T, v: [4]T, w: [4]T, rows: [3]usize) T {
    @setFloatMode(.strict);
    const columns = [3][4]T{ u, v, w };
    const m: [3][3]T = .{
        .{ columns[0][rows[0]], columns[1][rows[0]], columns[2][rows[0]] },
        .{ columns[0][rows[1]], columns[1][rows[1]], columns[2][rows[1]] },
        .{ columns[0][rows[2]], columns[1][rows[2]], columns[2][rows[2]] },
    };
    return m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]) +
        m[1][2] * (m[0][1] * m[2][0] - m[0][0] * m[2][1]) +
        m[2][2] * (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
}

/// The sign of `dot(a - o, n)` for a `dim`-wide pair, which is the signed
/// distance of `a` from the flat through `o` with unit normal `n`.
fn flatSideSign(comptime dim: usize, comptime T: type, a: [dim]T, o: [dim]T, n: [dim]T, policy: PredicatePolicy) Sign {
    var differences: [dim]f64 = undefined;
    var normals: [dim]f64 = undefined;
    inline for (0..dim) |index| {
        differences[index] = widen(T, a[index]) - widen(T, o[index]);
        normals[index] = widen(T, n[index]);
    }
    return switch (policy) {
        .fast => signOfScalar(T, flatDistanceFast(dim, T, a, o, n)),
        .adaptive_exact => exactSumOfProducts(dim, differences, normals).sign(),
    };
}

/// The signed distance of `a` from the flat through `o` with unit normal `n`,
/// the plain scalar evaluation whose sign `fast` reports. The terms are
/// summed left to right in declaration order.
fn flatDistanceFast(comptime dim: usize, comptime T: type, a: [dim]T, o: [dim]T, n: [dim]T) T {
    @setFloatMode(.strict);
    var sum: T = 0;
    inline for (0..dim) |index| {
        sum = sum + (a[index] - o[index]) * n[index];
    }
    return sum;
}

fn dotOfArray(comptime dim: usize, comptime T: type, a: [dim]T, b: [dim]T) T {
    @setFloatMode(.strict);
    var sum: T = 0;
    inline for (0..dim) |index| {
        sum = sum + a[index] * b[index];
    }
    return sum;
}

// ---------------------------------------------------------------------------
// Vector
// ---------------------------------------------------------------------------

/// The raw scalar kernels every arithmetic helper composes. Each sets strict
/// float mode, which is what stops the compiler from reassociating a sum or
/// contracting a product into the addition that follows it. The `fast`
/// predicate policy is DEFINED as the sign these produce, so a contraction
/// would silently change what `fast` means.
fn addRaw(comptime T: type, a: T, b: T) T {
    @setFloatMode(.strict);
    return a + b;
}

fn subRaw(comptime T: type, a: T, b: T) T {
    @setFloatMode(.strict);
    return a - b;
}

fn negRaw(comptime T: type, a: T) T {
    @setFloatMode(.strict);
    return -a;
}

fn scaleRaw(comptime T: type, a: T, b: T) T {
    @setFloatMode(.strict);
    return a * b;
}

fn divRaw(comptime T: type, a: T, b: T) T {
    @setFloatMode(.strict);
    return a / b;
}

/// Which componentwise operation a vector helper applies. A comptime tag
/// rather than a function value, because the raw kernels are generic over the
/// scalar type and a generic function does not coerce to a non-generic
/// function pointer.
const BinaryOp = enum { add, sub };
const UnaryOp = enum { neg };

fn applyBinary(comptime op: BinaryOp, comptime T: type, a: T, b: T) T {
    return switch (op) {
        .add => addRaw(T, a, b),
        .sub => subRaw(T, a, b),
    };
}

fn applyUnary(comptime op: UnaryOp, comptime T: type, a: T) T {
    return switch (op) {
        .neg => negRaw(T, a),
    };
}

fn vectorBinary(comptime V: type, a: V, b: V, comptime op: BinaryOp) error{NonFiniteInput}!V {
    try requireFiniteField(a);
    try requireFiniteField(b);
    const width = comptime fieldCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = comptime fieldName(index);
        @field(out, name) = applyBinary(op, componentType(V), @field(a, name), @field(b, name));
    }
    return out;
}

/// The same operation on two values the caller has already validated. Every
/// internal caller builds its operands from validated components, so the check
/// would be a second pass over data that cannot fail.
fn vectorBinaryQuiet(comptime V: type, a: V, b: V, comptime op: BinaryOp) V {
    return vectorBinary(V, a, b, op) catch unreachable;
}

fn vectorUnary(comptime V: type, a: V, comptime op: UnaryOp) error{NonFiniteInput}!V {
    try requireFiniteField(a);
    const width = comptime fieldCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = comptime fieldName(index);
        @field(out, name) = applyUnary(op, componentType(V), @field(a, name));
    }
    return out;
}

fn vectorUnaryQuiet(comptime V: type, a: V, comptime op: UnaryOp) V {
    return vectorUnary(V, a, op) catch unreachable;
}

fn vectorScaled(comptime V: type, a: V, factor: componentType(V)) error{NonFiniteInput}!V {
    try requireFiniteField(a);
    try requireFiniteScalar(componentType(V), factor);
    const width = comptime fieldCount(V);
    var out: V = undefined;
    inline for (0..width) |index| {
        const name = comptime fieldName(index);
        @field(out, name) = scaleRaw(componentType(V), @field(a, name), factor);
    }
    return out;
}

fn vectorScaledQuiet(comptime V: type, a: V, factor: componentType(V)) V {
    return vectorScaled(V, a, factor) catch unreachable;
}

/// The scalar product, summed left to right in declaration order. The order is
/// the contract: `fast` is defined to report the sign of exactly this
/// evaluation, so a caller that groups the terms differently gets a different
/// `fast` sign. That is why the order is written down and not left to the
/// optimizer.
fn vectorDot(comptime V: type, a: V, b: V) componentType(V) {
    @setFloatMode(.strict);
    const T = componentType(V);
    const width = comptime fieldCount(V);
    var sum: T = 0;
    inline for (0..width) |index| {
        const name = comptime fieldName(index);
        sum = sum + (@field(a, name) * @field(b, name));
    }
    return sum;
}

fn vectorLengthSquared(comptime V: type, a: V) componentType(V) {
    return vectorDot(V, a, a);
}

/// The length, folded through the accepted `hypot` kernel from a positive
/// zero. The fold cannot overflow while the components cannot, and it cannot
/// flush a representable norm to zero, which `sqrt(x * x + y * y)` does in
/// both directions. Nesting a two-argument kernel is not correctly rounded for
/// three or four components; it is the accepted kernel applied again, never a
/// re-derived polynomial.
fn vectorLength(comptime V: type, a: V) componentType(V) {
    const T = componentType(V);
    var accumulator: T = 0;
    inline for (toArrayOf(V, a)) |component| {
        accumulator = hypotOf(T, accumulator, component);
    }
    return accumulator;
}

/// The unit vector along `a`. A zero vector has no direction and an infinite
/// component has none either, so both are errors and neither is a NaN.
fn vectorNormalize(comptime V: type, a: V) error{ NonFiniteInput, Degenerate }!V {
    try requireFiniteField(a);
    const length = vectorLength(V, a);
    if (length == 0) return error.Degenerate;
    var out: V = undefined;
    const width = comptime fieldCount(V);
    inline for (0..width) |index| {
        @field(out, fieldName(index)) = divRaw(componentType(V), @field(a, fieldName(index)), length);
    }
    return out;
}

fn vectorDistanceSquared(comptime V: type, a: V, b: V) componentType(V) {
    const difference = vectorBinaryQuiet(V, a, b, .sub);
    return vectorDot(V, difference, difference);
}

fn vectorDistance(comptime V: type, a: V, b: V) componentType(V) {
    return vectorLength(V, vectorBinaryQuiet(V, a, b, .sub));
}

fn cross3Raw(comptime T: type, a: [3]T, b: [3]T) [3]T {
    @setFloatMode(.strict);
    return .{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

/// A two-component `f32` vector.
pub const Vector2f = Vector2(f32);
/// A two-component `f64` vector.
pub const Vector2d = Vector2(f64);
/// A three-component `f32` vector.
pub const Vector3f = Vector3(f32);
/// A three-component `f64` vector.
pub const Vector3d = Vector3(f64);
/// A four-component `f32` vector.
pub const Vector4f = Vector4(f32);
/// A four-component `f64` vector.
pub const Vector4d = Vector4(f64);

fn VectorOf(comptime dim: usize, comptime T: type) type {
    return switch (dim) {
        2 => Vector2(T),
        3 => Vector3(T),
        4 => Vector4(T),
        else => @compileError("zmath geometry dimensions are 2, 3, and 4"),
    };
}

fn Vector2(comptime T: type) type {
    return struct {
        /// The first component.
        x: T,
        /// The second component.
        y: T,

        const Self = @This();

        /// Builds a vector from named components.
        pub fn init(x: T, y: T) error{NonFiniteInput}!Self {
            const out = Self{ .x = x, .y = y };
            try requireFiniteField(out);
            return out;
        }

        /// Builds a vector from `[x, y]`.
        pub fn fromArray(a: [2]T) error{NonFiniteInput}!Self {
            return init(a[0], a[1]);
        }

        /// The components as `[x, y]`.
        pub fn toArray(self: Self) [2]T {
            return toArrayOf(Self, self);
        }

        /// Componentwise sum.
        pub fn add(self: Self, other: Self) error{NonFiniteInput}!Self {
            return vectorBinary(Self, self, other, .add);
        }

        /// Componentwise difference.
        pub fn sub(self: Self, other: Self) error{NonFiniteInput}!Self {
            return vectorBinary(Self, self, other, .sub);
        }

        /// Componentwise negation.
        pub fn neg(self: Self) error{NonFiniteInput}!Self {
            return vectorUnary(Self, self, .neg);
        }

        /// Every component times `factor`.
        pub fn scale(self: Self, factor: T) error{NonFiniteInput}!Self {
            return vectorScaled(Self, self, factor);
        }

        /// The scalar product, summed left to right.
        pub fn dot(self: Self, other: Self) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return vectorDot(Self, self, other);
        }

        /// The squared length, which is the dot product with itself.
        pub fn lengthSquared(self: Self) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            return vectorLengthSquared(Self, self);
        }

        /// The length, folded through the accepted `hypot` kernel.
        pub fn length(self: Self) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            return vectorLength(Self, self);
        }

        /// The unit vector along `self`.
        pub fn normalize(self: Self) error{ NonFiniteInput, Degenerate }!Self {
            return vectorNormalize(Self, self);
        }

        /// The squared distance to `point`.
        pub fn distanceSquared(self: Self, point: Point2(T)) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return vectorDistanceSquared(Self, self, toVectorOfPoint(Point2(T), point));
        }

        /// The distance to `point`, the length of the rounded difference.
        pub fn distance(self: Self, point: Point2(T)) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return vectorDistance(Self, self, toVectorOfPoint(Point2(T), point));
        }

        /// True when both components are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Self, self);
        }

        /// Bit equality, with every NaN equal to every other NaN.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(Self, self, other);
        }
    };
}

fn Vector3(comptime T: type) type {
    return struct {
        /// The first component.
        x: T,
        /// The second component.
        y: T,
        /// The third component.
        z: T,

        const Self = @This();

        /// Builds a vector from named components.
        pub fn init(x: T, y: T, z: T) error{NonFiniteInput}!Self {
            const out = Self{ .x = x, .y = y, .z = z };
            try requireFiniteField(out);
            return out;
        }

        /// Builds a vector from `[x, y, z]`.
        pub fn fromArray(a: [3]T) error{NonFiniteInput}!Self {
            return init(a[0], a[1], a[2]);
        }

        /// The components as `[x, y, z]`.
        pub fn toArray(self: Self) [3]T {
            return toArrayOf(Self, self);
        }

        /// Componentwise sum.
        pub fn add(self: Self, other: Self) error{NonFiniteInput}!Self {
            return vectorBinary(Self, self, other, .add);
        }

        /// Componentwise difference.
        pub fn sub(self: Self, other: Self) error{NonFiniteInput}!Self {
            return vectorBinary(Self, self, other, .sub);
        }

        /// Componentwise negation.
        pub fn neg(self: Self) error{NonFiniteInput}!Self {
            return vectorUnary(Self, self, .neg);
        }

        /// Every component times `factor`.
        pub fn scale(self: Self, factor: T) error{NonFiniteInput}!Self {
            return vectorScaled(Self, self, factor);
        }

        /// The scalar product, summed left to right.
        pub fn dot(self: Self, other: Self) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return vectorDot(Self, self, other);
        }

        /// The squared length, which is the dot product with itself.
        pub fn lengthSquared(self: Self) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            return vectorLengthSquared(Self, self);
        }

        /// The length, folded through the accepted `hypot` kernel.
        pub fn length(self: Self) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            return vectorLength(Self, self);
        }

        /// The unit vector along `self`.
        pub fn normalize(self: Self) error{ NonFiniteInput, Degenerate }!Self {
            return vectorNormalize(Self, self);
        }

        /// The right-handed cross product. It exists in three dimensions and
        /// nowhere else, because only there does the orthogonal complement of
        /// the plane two vectors span come out as a single vector.
        pub fn cross(self: Self, other: Self) error{NonFiniteInput}!Self {
            try requireFiniteField(self);
            try requireFiniteField(other);
            const product = cross3Raw(T, self.toArray(), other.toArray());
            return fromArrayOf(Self, product);
        }

        /// The squared distance to `point`.
        pub fn distanceSquared(self: Self, point: Point3(T)) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return vectorDistanceSquared(Self, self, toVectorOfPoint(Point3(T), point));
        }

        /// The distance to `point`, the length of the rounded difference.
        pub fn distance(self: Self, point: Point3(T)) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return vectorDistance(Self, self, toVectorOfPoint(Point3(T), point));
        }

        /// True when every component is finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Self, self);
        }

        /// Bit equality, with every NaN equal to every other NaN.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(Self, self, other);
        }
    };
}

fn Vector4(comptime T: type) type {
    return struct {
        /// The first component.
        x: T,
        /// The second component.
        y: T,
        /// The third component.
        z: T,
        /// The fourth component.
        w: T,

        const Self = @This();

        /// Builds a vector from named components.
        pub fn init(x: T, y: T, z: T, w: T) error{NonFiniteInput}!Self {
            const out = Self{ .x = x, .y = y, .z = z, .w = w };
            try requireFiniteField(out);
            return out;
        }

        /// Builds a vector from `[x, y, z, w]`.
        pub fn fromArray(a: [4]T) error{NonFiniteInput}!Self {
            return init(a[0], a[1], a[2], a[3]);
        }

        /// The components as `[x, y, z, w]`.
        pub fn toArray(self: Self) [4]T {
            return toArrayOf(Self, self);
        }

        /// Componentwise sum.
        pub fn add(self: Self, other: Self) error{NonFiniteInput}!Self {
            return vectorBinary(Self, self, other, .add);
        }

        /// Componentwise difference.
        pub fn sub(self: Self, other: Self) error{NonFiniteInput}!Self {
            return vectorBinary(Self, self, other, .sub);
        }

        /// Componentwise negation.
        pub fn neg(self: Self) error{NonFiniteInput}!Self {
            return vectorUnary(Self, self, .neg);
        }

        /// Every component times `factor`.
        pub fn scale(self: Self, factor: T) error{NonFiniteInput}!Self {
            return vectorScaled(Self, self, factor);
        }

        /// The scalar product, summed left to right.
        pub fn dot(self: Self, other: Self) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            try requireFiniteField(other);
            return vectorDot(Self, self, other);
        }

        /// The squared length, which is the dot product with itself.
        pub fn lengthSquared(self: Self) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            return vectorLengthSquared(Self, self);
        }

        /// The length, folded through the accepted `hypot` kernel.
        pub fn length(self: Self) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            return vectorLength(Self, self);
        }

        /// The unit vector along `self`.
        pub fn normalize(self: Self) error{ NonFiniteInput, Degenerate }!Self {
            return vectorNormalize(Self, self);
        }

        /// The exterior product `self ^ other`, which is the product four
        /// dimensions has. There is no cross product here: the orthogonal
        /// complement of the plane two vectors span is a whole plane, not a
        /// line, so no single vector completes the product.
        pub fn wedge(self: Self, other: Self) error{NonFiniteInput}!Bivector4(T) {
            try requireFiniteField(self);
            try requireFiniteField(other);
            const a = self.toArray();
            const b = other.toArray();
            const out = Bivector4(T){
                .xy = a[0] * b[1] - a[1] * b[0],
                .xz = a[0] * b[2] - a[2] * b[0],
                .xw = a[0] * b[3] - a[3] * b[0],
                .yz = a[1] * b[2] - a[2] * b[1],
                .yw = a[1] * b[3] - a[3] * b[1],
                .zw = a[2] * b[3] - a[3] * b[2],
            };
            try requireFiniteBivector(T, out);
            return out;
        }

        /// The squared distance to `point`.
        pub fn distanceSquared(self: Self, point: Point4(T)) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return vectorDistanceSquared(Self, self, toVectorOfPoint(Point4(T), point));
        }

        /// The distance to `point`, the length of the rounded difference.
        pub fn distance(self: Self, point: Point4(T)) error{NonFiniteInput}!T {
            try requireFiniteField(self);
            try requireFiniteField(point);
            return vectorDistance(Self, self, toVectorOfPoint(Point4(T), point));
        }

        /// True when every component is finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Self, self);
        }

        /// Bit equality, with every NaN equal to every other NaN.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(Self, self, other);
        }
    };
}

/// The exterior product of two four-dimensional vectors, in the six
/// independent components the order `xy, xz, xw, yz, yw, zw` gives.
pub const Bivector4f = Bivector4(f32);
/// The exterior product of two four-dimensional vectors, `f64`.
pub const Bivector4d = Bivector4(f64);

/// The six-component finiteness gate for a bivector. A bivector names its
/// components after the axis pairs they belong to rather than by an axis, so
/// it cannot share the axis-named helper the point and vector types use.
fn requireFiniteBivector(comptime T: type, value: Bivector4(T)) error{NonFiniteInput}!void {
    if (!isFiniteOf(T, value.xy) or !isFiniteOf(T, value.xz) or
        !isFiniteOf(T, value.xw) or !isFiniteOf(T, value.yz) or
        !isFiniteOf(T, value.yw) or !isFiniteOf(T, value.zw))
    {
        return error.NonFiniteInput;
    }
}

fn Bivector4(comptime T: type) type {
    return struct {
        /// The `x y` component.
        xy: T,
        /// The `x z` component.
        xz: T,
        /// The `x w` component.
        xw: T,
        /// The `y z` component.
        yz: T,
        /// The `y w` component.
        yw: T,
        /// The `z w` component.
        zw: T,

        const Self = @This();

        /// Builds a bivector from its six components.
        pub fn init(xy: T, xz: T, xw: T, yz: T, yw: T, zw: T) error{NonFiniteInput}!Self {
            const out = Self{ .xy = xy, .xz = xz, .xw = xw, .yz = yz, .yw = yw, .zw = zw };
            try requireFiniteBivector(T, out);
            return out;
        }

        /// Componentwise negation.
        pub fn neg(self: Self) Self {
            return .{
                .xy = negRaw(T, self.xy),
                .xz = negRaw(T, self.xz),
                .xw = negRaw(T, self.xw),
                .yz = negRaw(T, self.yz),
                .yw = negRaw(T, self.yw),
                .zw = negRaw(T, self.zw),
            };
        }

        /// True when every component is finite. The bivector's six components
        /// are named after the pairs they belong to rather than by an axis, so
        /// it does not share the axis-named helper the point and vector types
        /// use.
        pub fn isFinite(self: Self) bool {
            return isFiniteOf(T, self.xy) and isFiniteOf(T, self.xz) and
                isFiniteOf(T, self.xw) and isFiniteOf(T, self.yz) and
                isFiniteOf(T, self.yw) and isFiniteOf(T, self.zw);
        }
    };
}

// ---------------------------------------------------------------------------
// Point
// ---------------------------------------------------------------------------

fn toVectorOfPoint(comptime P: type, point: P) VectorOf(fieldCount(P), componentType(P)) {
    const V = VectorOf(fieldCount(P), componentType(P));
    return fromArrayOf(V, toArrayOf(P, point));
}

fn pointTranslate(comptime P: type, point: P, offset: VectorOf(fieldCount(P), componentType(P))) error{NonFiniteInput}!P {
    try requireFiniteField(point);
    try requireFiniteField(offset);
    return pointTranslateRaw(P, point, offset);
}

fn pointTranslateRaw(comptime P: type, point: P, offset: VectorOf(fieldCount(P), componentType(P))) P {
    const T = componentType(P);
    const components = toArrayOf(P, point);
    const deltas = toArrayOf(@TypeOf(offset), offset);
    var out: [fieldCount(P)]T = undefined;
    inline for (0..comptime fieldCount(P)) |index| {
        out[index] = addRaw(T, components[index], deltas[index]);
    }
    return fromArrayOf(P, out);
}

fn pointAtParameter(comptime P: type, origin: P, direction: VectorOf(fieldCount(P), componentType(P)), t: componentType(P)) P {
    return pointTranslateRaw(P, origin, vectorScaledQuiet(@TypeOf(direction), direction, t));
}

fn pointNegatedVector(comptime V: type, value: V) V {
    return vectorUnaryQuiet(V, value, .neg);
}

/// A two-dimensional `f32` point.
pub const Point2f = Point2(f32);
/// A two-dimensional `f64` point.
pub const Point2d = Point2(f64);
/// A three-dimensional `f32` point.
pub const Point3f = Point3(f32);
/// A three-dimensional `f64` point.
pub const Point3d = Point3(f64);
/// A four-dimensional `f32` point.
pub const Point4f = Point4(f32);
/// A four-dimensional `f64` point.
pub const Point4d = Point4(f64);

fn PointOf(comptime dim: usize, comptime T: type) type {
    return switch (dim) {
        2 => Point2(T),
        3 => Point3(T),
        4 => Point4(T),
        else => @compileError("zmath geometry dimensions are 2, 3, and 4"),
    };
}

fn Point2(comptime T: type) type {
    return struct {
        /// The first component.
        x: T,
        /// The second component.
        y: T,

        const Self = @This();

        /// Builds a point from named components.
        pub fn init(x: T, y: T) error{NonFiniteInput}!Self {
            const out = Self{ .x = x, .y = y };
            try requireFiniteField(out);
            return out;
        }

        /// Builds a point from `[x, y]`.
        pub fn fromArray(a: [2]T) error{NonFiniteInput}!Self {
            return init(a[0], a[1]);
        }

        /// The components as `[x, y]`.
        pub fn toArray(self: Self) [2]T {
            return toArrayOf(Self, self);
        }

        /// The point `offset` from `self`.
        pub fn add(self: Self, offset: Vector2(T)) error{NonFiniteInput}!Self {
            return pointTranslate(Self, self, offset);
        }

        /// The point `offset` back from `self`.
        pub fn sub(self: Self, offset: Vector2(T)) error{NonFiniteInput}!Self {
            return pointTranslate(Self, self, pointNegatedVector(Vector2(T), offset));
        }

        /// The same coordinates read as a vector.
        pub fn toVector(self: Self) Vector2(T) {
            return toVectorOfPoint(Self, self);
        }

        /// The sign of the oriented area of the triangle `a`, `b`, `c`: it is
        /// `negative`, `zero`, or `positive` as `c` lies to the right of, on,
        /// or to the left of the directed line from `a` to `b`.
        pub fn orientation(a: Self, b: Self, c: Self, policy: PredicatePolicy) error{NonFiniteInput}!Sign {
            try requireFiniteField(a);
            try requireFiniteField(b);
            try requireFiniteField(c);
            const av = a.toArray();
            const bv = b.toArray();
            const cv = c.toArray();
            return switch (policy) {
                .fast => signOfScalar(T, (bv[0] - av[0]) * (cv[1] - av[1]) -
                    (bv[1] - av[1]) * (cv[0] - av[0])),
                .adaptive_exact => exactOrient2(T, av, bv, cv),
            };
        }

        /// Which side of the directed line from `self` to `other` the `point`
        /// lies on. A two-dimensional line has two sides; a line in three or
        /// four dimensions does not, so this exists here only.
        pub fn side(self: Self, other: Self, point: Self, policy: PredicatePolicy) error{NonFiniteInput}!Sign {
            return orientation(self, other, point, policy);
        }

        /// True when both components are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Self, self);
        }

        /// Bit equality, with every NaN equal to every other NaN.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(Self, self, other);
        }
    };
}

fn Point3(comptime T: type) type {
    return struct {
        /// The first component.
        x: T,
        /// The second component.
        y: T,
        /// The third component.
        z: T,

        const Self = @This();

        /// Builds a point from named components.
        pub fn init(x: T, y: T, z: T) error{NonFiniteInput}!Self {
            const out = Self{ .x = x, .y = y, .z = z };
            try requireFiniteField(out);
            return out;
        }

        /// Builds a point from `[x, y, z]`.
        pub fn fromArray(a: [3]T) error{NonFiniteInput}!Self {
            return init(a[0], a[1], a[2]);
        }

        /// The components as `[x, y, z]`.
        pub fn toArray(self: Self) [3]T {
            return toArrayOf(Self, self);
        }

        /// The point `offset` from `self`.
        pub fn add(self: Self, offset: Vector3(T)) error{NonFiniteInput}!Self {
            return pointTranslate(Self, self, offset);
        }

        /// The point `offset` back from `self`.
        pub fn sub(self: Self, offset: Vector3(T)) error{NonFiniteInput}!Self {
            return pointTranslate(Self, self, pointNegatedVector(Vector3(T), offset));
        }

        /// The same coordinates read as a vector.
        pub fn toVector(self: Self) Vector3(T) {
            return toVectorOfPoint(Self, self);
        }

        /// The sign of the signed volume of the tetrahedron `a`, `b`, `c`,
        /// `d`, which is `dot(cross(b - a, c - a), d - a)`. It is zero
        /// exactly when the four points are coplanar, and that is the decision
        /// a plain determinant evaluation gets wrong when the terms cancel.
        pub fn orientation(a: Self, b: Self, c: Self, d: Self, policy: PredicatePolicy) error{NonFiniteInput}!Sign {
            try requireFiniteField(a);
            try requireFiniteField(b);
            try requireFiniteField(c);
            try requireFiniteField(d);
            const av = a.toArray();
            const bv = b.toArray();
            const cv = c.toArray();
            const dv = d.toArray();
            return switch (policy) {
                .fast => signOfScalar(T, crossDotFast(T, av, bv, cv, dv)),
                .adaptive_exact => exactOrient3(T, av, bv, cv, dv),
            };
        }

        /// True when every component is finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Self, self);
        }

        /// Bit equality, with every NaN equal to every other NaN.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(Self, self, other);
        }
    };
}

/// `dot(cross(b - a, c - a), d - a)`, summed left to right. The grouping is
/// the definition of `fast` for this predicate.
fn crossDotFast(comptime T: type, a: [3]T, b: [3]T, c: [3]T, d: [3]T) T {
    @setFloatMode(.strict);
    const n = cross3Raw(T, .{
        b[0] - a[0],
        b[1] - a[1],
        b[2] - a[2],
    }, .{
        c[0] - a[0],
        c[1] - a[1],
        c[2] - a[2],
    });
    const w = [3]T{ d[0] - a[0], d[1] - a[1], d[2] - a[2] };
    return n[0] * w[0] + n[1] * w[1] + n[2] * w[2];
}

fn Point4(comptime T: type) type {
    return struct {
        /// The first component.
        x: T,
        /// The second component.
        y: T,
        /// The third component.
        z: T,
        /// The fourth component.
        w: T,

        const Self = @This();

        /// Builds a point from named components.
        pub fn init(x: T, y: T, z: T, w: T) error{NonFiniteInput}!Self {
            const out = Self{ .x = x, .y = y, .z = z, .w = w };
            try requireFiniteField(out);
            return out;
        }

        /// Builds a point from `[x, y, z, w]`.
        pub fn fromArray(a: [4]T) error{NonFiniteInput}!Self {
            return init(a[0], a[1], a[2], a[3]);
        }

        /// The components as `[x, y, z, w]`.
        pub fn toArray(self: Self) [4]T {
            return toArrayOf(Self, self);
        }

        /// The point `offset` from `self`.
        pub fn add(self: Self, offset: Vector4(T)) error{NonFiniteInput}!Self {
            return pointTranslate(Self, self, offset);
        }

        /// The point `offset` back from `self`.
        pub fn sub(self: Self, offset: Vector4(T)) error{NonFiniteInput}!Self {
            return pointTranslate(Self, self, pointNegatedVector(Vector4(T), offset));
        }

        /// The same coordinates read as a vector.
        pub fn toVector(self: Self) Vector4(T) {
            return toVectorOfPoint(Self, self);
        }

        /// The sign of the four-by-four determinant whose rows are the four
        /// points. It is zero exactly when the four position vectors are
        /// linearly dependent, and its magnitude is the oriented
        /// four-dimensional volume they span from the origin. This is the
        /// degree-four predicate that sizes the exact arithmetic: its value
        /// needs every limb of the array.
        ///
        /// It is not the coplanarity test. Four four-dimensional points are
        /// coplanar through a different form, a rank test on three difference
        /// vectors, and the affine intersections use that one privately.
        pub fn orientation(a: Self, b: Self, c: Self, d: Self, policy: PredicatePolicy) error{NonFiniteInput}!Sign {
            try requireFiniteField(a);
            try requireFiniteField(b);
            try requireFiniteField(c);
            try requireFiniteField(d);
            return orient4Sign(T, a.toArray(), b.toArray(), c.toArray(), d.toArray(), policy);
        }

        /// True when every component is finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Self, self);
        }

        /// Bit equality, with every NaN equal to every other NaN.
        pub fn eq(self: Self, other: Self) bool {
            return eqOf(Self, self, other);
        }
    };
}

// ---------------------------------------------------------------------------
// Line, ray, segment
// ---------------------------------------------------------------------------

fn lineBuild(comptime P: type, point: P, direction: VectorOf(fieldCount(P), componentType(P))) error{ NonFiniteInput, Degenerate }!struct { origin: P, direction: VectorOf(fieldCount(P), componentType(P)) } {
    try requireFiniteField(point);
    const unit = try vectorNormalize(@TypeOf(direction), direction);
    return .{ .origin = point, .direction = unit };
}

fn projectOntoLine(comptime P: type, origin: P, direction: VectorOf(fieldCount(P), componentType(P)), point: P) error{NonFiniteInput}!P {
    try requireFiniteField(point);
    const V = @TypeOf(direction);
    const difference = vectorBinaryQuiet(V, toVectorOfPoint(P, point), toVectorOfPoint(P, origin), .sub);
    const t = vectorDot(V, difference, direction);
    return pointTranslateRaw(P, origin, vectorScaledQuiet(V, direction, t));
}

/// A two-dimensional `f32` line.
pub const Line2f = LineOf(2, f32);
/// A two-dimensional `f64` line.
pub const Line2d = LineOf(2, f64);
/// A three-dimensional `f32` line.
pub const Line3f = LineOf(3, f32);
/// A three-dimensional `f64` line.
pub const Line3d = LineOf(3, f64);
/// A four-dimensional `f32` line.
pub const Line4f = LineOf(4, f32);
/// A four-dimensional `f64` line.
pub const Line4d = LineOf(4, f64);

fn LineOf(comptime dim: usize, comptime T: type) type {
    const P = PointOf(dim, T);
    const V = VectorOf(dim, T);
    return struct {
        /// A point on the line.
        origin: P,
        /// The unit direction.
        direction: V,

        const Self = @This();

        /// Builds the line through `point` along `direction`, which is
        /// normalized. A zero direction has no line to build.
        pub fn init(point: P, direction: V) error{ NonFiniteInput, Degenerate }!Self {
            const built = try lineBuild(P, point, direction);
            return .{ .origin = built.origin, .direction = built.direction };
        }

        /// The point at `t`. There is no clamp: a line runs in both
        /// directions.
        pub fn at(self: Self, t: T) error{NonFiniteInput}!P {
            try requireFiniteScalar(T, t);
            return pointAtParameter(P, self.origin, self.direction, t);
        }

        /// The point of the line closest to `point`, which is the orthogonal
        /// projection.
        pub fn closestPointTo(self: Self, point: P) error{NonFiniteInput}!P {
            return projectOntoLine(P, self.origin, self.direction, point);
        }

        /// The set of points common to two lines. Two lines in two or more
        /// dimensions share at most one point unless they are the same line,
        /// so the answer is `none`, `point`, or `line`.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{NonFiniteInput}!IntersectionOf(dim, T) {
            try requireFiniteField(self.origin);
            try requireFiniteField(other.origin);
            return affineIntersectDispatch(dim, T, self.origin.toArray(), pointAtArray(dim, T, self.origin.toArray(), self.direction.toArray(), 1), -infinityOf(T), infinityOf(T), other.origin.toArray(), pointAtArray(dim, T, other.origin.toArray(), other.direction.toArray(), 1), -infinityOf(T), infinityOf(T), policy);
        }

        /// True when the origin and the direction are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(P, self.origin) and isFiniteQuiet(V, self.direction);
        }
    };
}

/// A two-dimensional `f32` ray.
pub const Ray2f = RayOf(2, f32);
/// A two-dimensional `f64` ray.
pub const Ray2d = RayOf(2, f64);
/// A three-dimensional `f32` ray.
pub const Ray3f = RayOf(3, f32);
/// A three-dimensional `f64` ray.
pub const Ray3d = RayOf(3, f64);
/// A four-dimensional `f32` ray.
pub const Ray4f = RayOf(4, f32);
/// A four-dimensional `f64` ray.
pub const Ray4d = RayOf(4, f64);

fn RayOf(comptime dim: usize, comptime T: type) type {
    const P = PointOf(dim, T);
    const V = VectorOf(dim, T);
    return struct {
        /// The end the ray starts from.
        origin: P,
        /// The unit direction the ray runs in.
        direction: V,

        const Self = @This();

        /// Builds the ray from `origin` along `direction`, which is
        /// normalized. A zero direction has no ray to build.
        pub fn init(origin: P, direction: V) error{ NonFiniteInput, Degenerate }!Self {
            const built = try lineBuild(P, origin, direction);
            return .{ .origin = built.origin, .direction = built.direction };
        }

        /// The point at `t`, clamped at `t = 0` so the ray never extends
        /// behind the end it starts from.
        pub fn at(self: Self, t: T) error{NonFiniteInput}!P {
            try requireFiniteScalar(T, t);
            if (t <= 0) return self.origin;
            return pointAtParameter(P, self.origin, self.direction, t);
        }

        /// The set of points common to two rays. Two rays share at most one
        /// point unless they leave the same end in the same direction, which
        /// is `contained` or `coincident`.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{NonFiniteInput}!IntersectionOf(dim, T) {
            try requireFiniteField(self.origin);
            try requireFiniteField(other.origin);
            return affineIntersectDispatch(dim, T, self.origin.toArray(), pointAtArray(dim, T, self.origin.toArray(), self.direction.toArray(), 1), 0, infinityOf(T), other.origin.toArray(), pointAtArray(dim, T, other.origin.toArray(), other.direction.toArray(), 1), 0, infinityOf(T), policy);
        }

        /// The set of points common to this ray and a line.
        pub fn intersectLine(self: Self, line: LineOf(dim, T), policy: PredicatePolicy) error{NonFiniteInput}!IntersectionOf(dim, T) {
            try requireFiniteField(self.origin);
            try requireFiniteField(line.origin);
            return affineIntersectDispatch(dim, T, self.origin.toArray(), pointAtArray(dim, T, self.origin.toArray(), self.direction.toArray(), 1), 0, infinityOf(T), line.origin.toArray(), pointAtArray(dim, T, line.origin.toArray(), line.direction.toArray(), 1), -infinityOf(T), infinityOf(T), policy);
        }

        /// The set of points common to this ray and a segment.
        pub fn intersectSegment(self: Self, segment: SegmentOf(dim, T), policy: PredicatePolicy) error{NonFiniteInput}!IntersectionOf(dim, T) {
            try requireFiniteField(self.origin);
            try requireFiniteField(segment);
            return affineIntersectDispatch(dim, T, self.origin.toArray(), pointAtArray(dim, T, self.origin.toArray(), self.direction.toArray(), 1), 0, infinityOf(T), segment.a.toArray(), segment.b.toArray(), 0, 1, policy);
        }

        /// True when the origin and the direction are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(P, self.origin) and isFiniteQuiet(V, self.direction);
        }
    };
}

/// A two-dimensional `f32` segment.
pub const Segment2f = SegmentOf(2, f32);
/// A two-dimensional `f64` segment.
pub const Segment2d = SegmentOf(2, f64);
/// A three-dimensional `f32` segment.
pub const Segment3f = SegmentOf(3, f32);
/// A three-dimensional `f64` segment.
pub const Segment3d = SegmentOf(3, f64);
/// A four-dimensional `f32` segment.
pub const Segment4f = SegmentOf(4, f32);
/// A four-dimensional `f64` segment.
pub const Segment4d = SegmentOf(4, f64);

fn SegmentOf(comptime dim: usize, comptime T: type) type {
    const P = PointOf(dim, T);
    const V = VectorOf(dim, T);
    return struct {
        /// The end at `t = 0`.
        a: P,
        /// The end at `t = 1`.
        b: P,

        const Self = @This();

        /// Builds the segment from `a` to `b`. The two ends may coincide,
        /// which gives the degenerate segment that is a single point.
        pub fn init(a: P, b: P) error{NonFiniteInput}!Self {
            try requireFiniteField(a);
            try requireFiniteField(b);
            return .{ .a = a, .b = b };
        }

        /// The point at `t`, clamped to `[0, 1]`, so the segment never
        /// extends past either end.
        pub fn at(self: Self, t: T) error{NonFiniteInput}!P {
            try requireFiniteScalar(T, t);
            const clamped = clampUnitOf(T, t);
            if (clamped == 0) return self.a;
            if (clamped == 1) return self.b;
            const difference = vectorBinaryQuiet(V, toVectorOfPoint(P, self.b), toVectorOfPoint(P, self.a), .sub);
            return pointTranslateRaw(P, self.a, vectorScaledQuiet(V, difference, clamped));
        }

        /// The distance between the two ends, which cannot overflow while the
        /// ends cannot.
        pub fn length(self: Self) error{NonFiniteInput}!T {
            return vectorDistance(V, toVectorOfPoint(P, self.a), toVectorOfPoint(P, self.b));
        }

        /// The point of the segment closest to `point`, which is the
        /// orthogonal projection clamped to the segment.
        pub fn closestPointTo(self: Self, point: P) error{NonFiniteInput}!P {
            try requireFiniteField(point);
            const from = toVectorOfPoint(P, self.a);
            const to = toVectorOfPoint(P, self.b);
            const difference = vectorBinaryQuiet(V, to, from, .sub);
            const denominator = vectorDot(V, difference, difference);
            if (denominator == 0) return self.a;
            const offset = vectorBinaryQuiet(V, toVectorOfPoint(P, point), from, .sub);
            const t = clampUnitOf(T, vectorDot(V, offset, difference) / denominator);
            if (t == 0) return self.a;
            if (t == 1) return self.b;
            return pointTranslateRaw(P, self.a, vectorScaledQuiet(V, difference, t));
        }

        /// The set of points common to two segments. The tags are `none` for
        /// disjoint segments, `point` for a crossing or a touching end,
        /// `segment` for a proper collinear overlap, `contained` when one
        /// segment lies inside the other, and `coincident` when the two are
        /// the same segment.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{NonFiniteInput}!IntersectionOf(dim, T) {
            try requireFiniteField(self.a);
            try requireFiniteField(other.a);
            return affineIntersectDispatch(dim, T, self.a.toArray(), self.b.toArray(), 0, 1, other.a.toArray(), other.b.toArray(), 0, 1, policy);
        }

        /// The set of points common to this segment and a ray.
        pub fn intersectRay(self: Self, ray: RayOf(dim, T), policy: PredicatePolicy) error{NonFiniteInput}!IntersectionOf(dim, T) {
            try requireFiniteField(self.a);
            try requireFiniteField(ray.origin);
            return affineIntersectDispatch(dim, T, self.a.toArray(), self.b.toArray(), 0, 1, ray.origin.toArray(), pointAtArray(dim, T, ray.origin.toArray(), ray.direction.toArray(), 1), 0, infinityOf(T), policy);
        }

        /// The set of points common to this segment and a line.
        pub fn intersectLine(self: Self, line: LineOf(dim, T), policy: PredicatePolicy) error{NonFiniteInput}!IntersectionOf(dim, T) {
            try requireFiniteField(self.a);
            try requireFiniteField(line.origin);
            return affineIntersectDispatch(dim, T, self.a.toArray(), self.b.toArray(), 0, 1, line.origin.toArray(), pointAtArray(dim, T, line.origin.toArray(), line.direction.toArray(), 1), -infinityOf(T), infinityOf(T), policy);
        }

        /// True when both ends are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(P, self.a) and isFiniteQuiet(P, self.b);
        }
    };
}

// ---------------------------------------------------------------------------
// Array helpers shared by the geometric construction paths
// ---------------------------------------------------------------------------

fn widenArray(comptime dim: usize, comptime T: type, a: [dim]T) [dim]f64 {
    var out: [dim]f64 = undefined;
    inline for (0..dim) |index| {
        out[index] = widen(T, a[index]);
    }
    return out;
}

fn arrayLength(comptime dim: usize, comptime T: type, a: [dim]T) T {
    var accumulator: T = 0;
    inline for (0..dim) |index| {
        accumulator = hypotOf(T, accumulator, a[index]);
    }
    return accumulator;
}

/// `a` scaled to unit length, or `a` unchanged when it has zero length. The
/// construction paths below call this on vectors that are provably nonzero,
/// so the fallback is a bound and not a result.
fn arrayNormalized(comptime dim: usize, comptime T: type, a: [dim]T) [dim]T {
    const length = arrayLength(dim, T, a);
    if (length == 0) return a;
    var out: [dim]T = undefined;
    inline for (0..dim) |index| {
        out[index] = divRaw(T, a[index], length);
    }
    return out;
}

fn arrayScale(comptime dim: usize, comptime T: type, a: [dim]T, factor: T) [dim]T {
    var out: [dim]T = undefined;
    inline for (0..dim) |index| {
        out[index] = scaleRaw(T, a[index], factor);
    }
    return out;
}

fn arrayAdd(comptime dim: usize, comptime T: type, a: [dim]T, b: [dim]T) [dim]T {
    var out: [dim]T = undefined;
    inline for (0..dim) |index| {
        out[index] = addRaw(T, a[index], b[index]);
    }
    return out;
}

fn arraySub(comptime dim: usize, comptime T: type, a: [dim]T, b: [dim]T) [dim]T {
    var out: [dim]T = undefined;
    inline for (0..dim) |index| {
        out[index] = subRaw(T, a[index], b[index]);
    }
    return out;
}

fn arrayDot(comptime dim: usize, comptime T: type, a: [dim]T, b: [dim]T) T {
    return dotOfArray(dim, T, a, b);
}

/// A unit vector orthogonal to `v`, built by projecting out the coordinate
/// axis least aligned with `v` and normalizing. Choosing the least aligned
/// axis is what keeps the residual far from zero, so the normalization is
/// well conditioned. The caller uses this to complete a frame, never to
/// decide a sign, so it needs no policy and no error.
fn anyPerpendicular(comptime dim: usize, comptime T: type, v: [dim]T) [dim]T {
    var best: usize = 0;
    var best_value: T = absOf(T, v[0]);
    inline for (1..dim) |index| {
        const value = absOf(T, v[index]);
        if (value < best_value) {
            best_value = value;
            best = index;
        }
    }
    var seed: [dim]T = undefined;
    inline for (0..dim) |index| {
        seed[index] = 0;
    }
    seed[best] = 1;
    const along = arrayDot(dim, T, seed, v);
    const squared = arrayDot(dim, T, v, v);
    return arrayNormalized(dim, T, arraySub(dim, T, seed, arrayScale(dim, T, v, divRaw(T, along, squared))));
}

/// `v` with the component along each vector of `against` removed, in the
/// order given, then normalized. A slice rather than a fixed-length array
/// because the number of vectors depends on the dimension, and an `inline
/// for` over a slice is not what a comptime frame needs.
fn orthogonalized(comptime dim: usize, comptime T: type, v: [dim]T, against: []const [dim]T) [dim]T {
    var current = v;
    for (against) |basis| {
        const projected = arrayDot(dim, T, current, basis);
        current = arraySub(dim, T, current, arrayScale(dim, T, basis, projected));
    }
    return arrayNormalized(dim, T, current);
}

fn pointAtArray(comptime dim: usize, comptime T: type, origin: [dim]T, direction: [dim]T, t: T) [dim]T {
    return arrayAdd(dim, T, origin, arrayScale(dim, T, direction, t));
}

// ---------------------------------------------------------------------------
// Degree-two predicate forms shared by the balls and the flats
// ---------------------------------------------------------------------------

/// The sign of `dot(v, v) - k * k`, degree two. This is the form behind ball
/// containment and behind every ball-to-ball comparison, so both policies
/// evaluate the same expression and differ only in the arithmetic.
fn squaredNormMinusSquareSign(comptime dim: usize, comptime T: type, v: [dim]T, k: T, policy: PredicatePolicy) Sign {
    return switch (policy) {
        .fast => signOfScalar(T, subRaw(T, dotOfArray(dim, T, v, v), scaleRaw(T, k, k))),
        .adaptive_exact => blk: {
            const wide = widenArray(dim, T, v);
            var total = exactSumOfProducts(dim, wide, wide);
            const radius = Big.fromFloat(widen(T, k));
            total = total.sub(radius.mul(radius));
            break :blk total.sign();
        },
    };
}

/// The sign of `hb * hb - (dot(w, w) - r * r)`, the ray-against-ball
/// discriminant, degree two.
fn rayBallDiscriminantSign(comptime dim: usize, comptime T: type, w: [dim]T, hb: T, r: T, policy: PredicatePolicy) Sign {
    return switch (policy) {
        .fast => signOfScalar(T, subRaw(T, scaleRaw(T, hb, hb), subRaw(T, dotOfArray(dim, T, w, w), scaleRaw(T, r, r)))),
        .adaptive_exact => blk: {
            const half = Big.fromFloat(widen(T, hb));
            const wide = widenArray(dim, T, w);
            const radius = Big.fromFloat(widen(T, r));
            const constant = exactSumOfProducts(dim, wide, wide).sub(radius.mul(radius));
            break :blk half.mul(half).sub(constant).sign();
        },
    };
}

/// The sign of `signed_distance * signed_distance - r * r`, the
/// ball-against-flat comparison, degree two. The distance is degree one and
/// the radius is degree zero, so the two squares carry the same scale and
/// need no adjustment.
fn flatBallSign(comptime T: type, distance: T, r: T, policy: PredicatePolicy) Sign {
    return switch (policy) {
        .fast => signOfScalar(T, subRaw(T, scaleRaw(T, distance, distance), scaleRaw(T, r, r))),
        .adaptive_exact => blk: {
            const d = Big.fromFloat(widen(T, distance));
            const radius = Big.fromFloat(widen(T, r));
            break :blk d.mul(d).sub(radius.mul(radius)).sign();
        },
    };
}

// ---------------------------------------------------------------------------
// A two-dimensional affine flat in four dimensions
// ---------------------------------------------------------------------------

/// A two-flat in four dimensions, which is what two hyperplanes meet in. The
/// `API.md` family list for four dimensions has no plane, but the four-
/// dimensional intersection union has to name the result of a hyperplane
/// meeting a hyperplane, and the name of that shape is a two-flat. It is
/// `origin` plus two orthonormal directions spanning it.
pub const Flat2f = struct {
    /// A point on the flat.
    origin: Point4f,
    /// The first unit direction spanning the flat.
    u: Vector4f,
    /// The second unit direction spanning the flat, orthogonal to `u`.
    v: Vector4f,
};

/// A two-flat in four dimensions, `f64`.
pub const Flat2d = struct {
    /// A point on the flat.
    origin: Point4d,
    /// The first unit direction spanning the flat.
    u: Vector4d,
    /// The second unit direction spanning the flat, orthogonal to `u`.
    v: Vector4d,
};

fn Flat2Of(comptime T: type) type {
    return if (T == f32) Flat2f else Flat2d;
}

// ---------------------------------------------------------------------------
// Intersection unions
// ---------------------------------------------------------------------------

/// The result of any two-dimensional intersection, in `f32`.
///
/// The tags, and why each is here:
///
/// - `none`: the two operands have no common point. Every pair can produce
///   this one, which is why it is a tag and not the absence of a tag.
/// - `point`: the common set is a single point. A proper crossing, a
///   touching endpoint, or tangency.
/// - `two_points`: two separate points. Two circles or two spheres that cut
///   each other in a circle of intersection reduce to the two points of the
///   line of centres that bound it; a segment or a ray through a ball that
///   goes right through it also lands here.
/// - `segment`: a proper sub-segment. Two collinear segments that overlap
///   without either containing the other.
/// - `line`: the common set is the whole line, which happens exactly when two
///   lines are the same line.
/// - `box`: the common set is an axis-aligned box that is neither operand.
///   Two boxes that overlap partially, including the touching case where the
///   box is degenerate.
/// - `contained`: the common set is exactly one operand. A segment inside a
///   flat, a segment inside a collinear segment, one circle strictly inside
///   another, one box inside another.
/// - `coincident`: the operands are the same set. Identical segments,
///   identical circles, identical boxes.
///
/// There is no `plane` tag here and no `hyperplane` tag: two dimensions have
/// neither, and a union that carried a tag it could never produce would be
/// a lie about the geometry. There is no `circle` tag either, because a flat
/// that could cut a circle needs a third dimension, and `API.md` puts the
/// circle in three dimensions.
pub const Intersection2f = Intersection2(f32);
/// The result of any two-dimensional intersection, in `f64`.
pub const Intersection2d = Intersection2(f64);
/// The result of any three-dimensional intersection, in `f32`.
///
/// The two-dimensional tags, and three more:
///
/// - `plane`: the common set is the whole plane, which happens exactly when
///   two planes are the same plane.
/// - `circle`: the common set is a circle. A plane that cuts a sphere not
///   tangent to it meets it in a circle, and no other pair in three dimensions
///   produces a curve.
///
/// There is no `hyperplane` tag: three dimensions have no four-flat, so a
/// hyperplane tag could never be produced.
pub const Intersection3f = Intersection3(f32);
/// The result of any three-dimensional intersection, in `f64`.
pub const Intersection3d = Intersection3(f64);
/// The result of any four-dimensional intersection, in `f32`.
///
/// The three-dimensional tags, and three more:
///
/// - `hyperplane`: the common set is a three-flat, which happens exactly when
///   two hyperplanes are the same hyperplane.
/// - `hypersphere`: the common set is a three-ball boundary. A hyperplane
///   that cuts a hypersphere not tangent to it meets it in one.
///
/// The `plane` tag here names a TWO-flat, not a three-dimensional plane: in
/// four dimensions the tag names the dimension of the common affine span, so
/// two hyperplanes that are not parallel meet in a two-flat and say so.
pub const Intersection4f = Intersection4(f32);
/// The result of any four-dimensional intersection, in `f64`.
pub const Intersection4d = Intersection4(f64);

fn IntersectionOf(comptime dim: usize, comptime T: type) type {
    return switch (dim) {
        2 => Intersection2(T),
        3 => Intersection3(T),
        4 => Intersection4(T),
        else => @compileError("zmath geometry dimensions are 2, 3, and 4"),
    };
}

fn Intersection2(comptime T: type) type {
    return union(enum) {
        /// No common point.
        none,
        /// One common point.
        point: Point2(T),
        /// Two common points.
        two_points: [2]Point2(T),
        /// A proper sub-segment of one of the operands.
        segment: SegmentOf(2, T),
        /// A whole line, which means the two lines are one line.
        line: LineOf(2, T),
        /// An axis-aligned box that is neither operand.
        box: Aabb2(T),
        /// One operand, whole.
        contained: ContainedSide,
        /// Both operands, the same set.
        coincident,
    };
}

fn Intersection3(comptime T: type) type {
    return union(enum) {
        /// No common point.
        none,
        /// One common point.
        point: Point3(T),
        /// Two common points.
        two_points: [2]Point3(T),
        /// A proper sub-segment of one of the operands.
        segment: SegmentOf(3, T),
        /// A whole line, which means the two lines are one line.
        line: LineOf(3, T),
        /// A whole plane, which means the two planes are one plane.
        plane: Plane3(T),
        /// A circle, from a plane cutting a sphere.
        circle: Sphere3(T),
        /// An axis-aligned box that is neither operand.
        box: Aabb3(T),
        /// One operand, whole.
        contained: ContainedSide,
        /// Both operands, the same set.
        coincident,
    };
}

fn Intersection4(comptime T: type) type {
    return union(enum) {
        /// No common point.
        none,
        /// One common point.
        point: Point4(T),
        /// Two common points.
        two_points: [2]Point4(T),
        /// A proper sub-segment of one of the operands.
        segment: SegmentOf(4, T),
        /// A whole line, which means the two lines are one line.
        line: LineOf(4, T),
        /// A two-flat, which is what two non-parallel hyperplanes meet in.
        plane: Flat2Of(T),
        /// A whole three-flat, which means the two hyperplanes are one.
        hyperplane: Hyperplane4(T),
        /// A three-ball boundary, from a hyperplane cutting a hypersphere.
        hypersphere: Hypersphere4(T),
        /// An axis-aligned box that is neither operand.
        box: Aabb4(T),
        /// One operand, whole.
        contained: ContainedSide,
        /// Both operands, the same set.
        coincident,
    };
}

// ---------------------------------------------------------------------------
// Circle, sphere, hypersphere
// ---------------------------------------------------------------------------

fn ballBuild(comptime P: type, center: P, radius: componentType(P)) error{ NonFiniteInput, NegativeRadius }!struct { center: P, radius: componentType(P) } {
    try requireFiniteField(center);
    const T = componentType(P);
    try requireFiniteScalar(T, radius);
    if (radius < 0) return error.NegativeRadius;
    return .{ .center = center, .radius = radius };
}

fn ballContainsSign(comptime dim: usize, comptime T: type, center: [dim]T, radius: T, point: [dim]T, policy: PredicatePolicy) Sign {
    return squaredNormMinusSquareSign(dim, T, arraySub(dim, T, point, center), radius, policy);
}

/// A two-dimensional `f32` circle.
pub const Circle2f = Circle2(f32);
/// A two-dimensional `f64` circle.
pub const Circle2d = Circle2(f64);
/// A three-dimensional `f32` sphere.
pub const Sphere3f = Sphere3(f32);
/// A three-dimensional `f64` sphere.
pub const Sphere3d = Sphere3(f64);
/// A four-dimensional `f32` hypersphere.
pub const Hypersphere4f = Hypersphere4(f32);
/// A four-dimensional `f64` hypersphere.
pub const Hypersphere4d = Hypersphere4(f64);

fn SphereOf(comptime dim: usize, comptime T: type) type {
    return switch (dim) {
        2 => Circle2(T),
        3 => Sphere3(T),
        4 => Hypersphere4(T),
        else => @compileError("zmath geometry dimensions are 2, 3, and 4"),
    };
}

fn Circle2(comptime T: type) type {
    return struct {
        /// The centre.
        center: Point2(T),
        /// The radius, finite and nonnegative.
        radius: T,

        const Self = @This();

        /// Builds the circle. The radius must be finite and nonnegative, and
        /// a negative one is an error rather than a sign flip.
        pub fn init(center: Point2(T), radius: T) error{ NonFiniteInput, NegativeRadius }!Self {
            const built = try ballBuild(Point2(T), center, radius);
            return .{ .center = built.center, .radius = built.radius };
        }

        /// True when `point` is inside or ON the circle. The boundary is
        /// inclusive, because a point at exactly the radius is on the set the
        /// circle names, and `adaptive_exact` is what makes that exact.
        pub fn contains(self: Self, point: Point2(T), policy: PredicatePolicy) error{NonFiniteInput}!bool {
            try requireFiniteField(self.center);
            try requireFiniteField(point);
            const sign = ballContainsSign(2, T, self.center.toArray(), self.radius, point.toArray(), policy);
            return sign == .zero or sign == .negative;
        }

        /// The set of points common to two circles.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{NonFiniteInput}!Intersection2(T) {
            try requireFiniteField(self.center);
            try requireFiniteField(other.center);
            return ballIntersect(2, T, self.center.toArray(), self.radius, other.center.toArray(), other.radius, policy);
        }

        /// The set of points common to this circle and a segment.
        pub fn intersectSegment(self: Self, segment: SegmentOf(2, T), policy: PredicatePolicy) error{NonFiniteInput}!Intersection2(T) {
            try requireFiniteField(self.center);
            try requireFiniteField(segment);
            return ballIntersectSegment(2, T, self.center.toArray(), self.radius, segment.a.toArray(), segment.b.toArray(), policy);
        }

        /// The set of points common to this circle and a ray.
        pub fn intersectRay(self: Self, ray: RayOf(2, T), policy: PredicatePolicy) error{NonFiniteInput}!Intersection2(T) {
            try requireFiniteField(self.center);
            try requireFiniteField(ray);
            return ballIntersectRay(2, T, self.center.toArray(), self.radius, ray.origin.toArray(), ray.direction.toArray(), policy);
        }

        /// True when the centre and the radius are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Point2(T), self.center) and isFiniteOf(T, self.radius);
        }
    };
}

fn Sphere3(comptime T: type) type {
    return struct {
        /// The centre.
        center: Point3(T),
        /// The radius, finite and nonnegative.
        radius: T,

        const Self = @This();

        /// Builds the sphere. The radius must be finite and nonnegative, and
        /// a negative one is an error rather than a sign flip.
        pub fn init(center: Point3(T), radius: T) error{ NonFiniteInput, NegativeRadius }!Self {
            const built = try ballBuild(Point3(T), center, radius);
            return .{ .center = built.center, .radius = built.radius };
        }

        /// True when `point` is inside or ON the sphere. The boundary is
        /// inclusive.
        pub fn contains(self: Self, point: Point3(T), policy: PredicatePolicy) error{NonFiniteInput}!bool {
            try requireFiniteField(self.center);
            try requireFiniteField(point);
            const sign = ballContainsSign(3, T, self.center.toArray(), self.radius, point.toArray(), policy);
            return sign == .zero or sign == .negative;
        }

        /// The set of points common to two spheres.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{NonFiniteInput}!Intersection3(T) {
            try requireFiniteField(self.center);
            try requireFiniteField(other.center);
            return ballIntersect(3, T, self.center.toArray(), self.radius, other.center.toArray(), other.radius, policy);
        }

        /// The set of points common to this sphere and a segment.
        pub fn intersectSegment(self: Self, segment: SegmentOf(3, T), policy: PredicatePolicy) error{NonFiniteInput}!Intersection3(T) {
            try requireFiniteField(self.center);
            try requireFiniteField(segment);
            return ballIntersectSegment(3, T, self.center.toArray(), self.radius, segment.a.toArray(), segment.b.toArray(), policy);
        }

        /// The set of points common to this sphere and a ray.
        pub fn intersectRay(self: Self, ray: RayOf(3, T), policy: PredicatePolicy) error{NonFiniteInput}!Intersection3(T) {
            try requireFiniteField(self.center);
            try requireFiniteField(ray);
            return ballIntersectRay(3, T, self.center.toArray(), self.radius, ray.origin.toArray(), ray.direction.toArray(), policy);
        }

        /// The set of points common to this sphere and a plane: nothing, the
        /// one tangent point, or the circle the plane cuts.
        pub fn intersectPlane(self: Self, plane: Plane3(T), policy: PredicatePolicy) error{NonFiniteInput}!Intersection3(T) {
            try requireFiniteField(self.center);
            try requireFiniteField(plane);
            return ballIntersectFlat(3, T, self.center.toArray(), self.radius, plane.origin.toArray(), plane.normal.toArray(), policy);
        }

        /// True when the centre and the radius are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Point3(T), self.center) and isFiniteOf(T, self.radius);
        }
    };
}

fn Hypersphere4(comptime T: type) type {
    return struct {
        /// The centre.
        center: Point4(T),
        /// The radius, finite and nonnegative.
        radius: T,

        const Self = @This();

        /// Builds the hypersphere. The radius must be finite and nonnegative,
        /// and a negative one is an error rather than a sign flip.
        pub fn init(center: Point4(T), radius: T) error{ NonFiniteInput, NegativeRadius }!Self {
            const built = try ballBuild(Point4(T), center, radius);
            return .{ .center = built.center, .radius = built.radius };
        }

        /// True when `point` is inside or ON the hypersphere. The boundary is
        /// inclusive.
        pub fn contains(self: Self, point: Point4(T), policy: PredicatePolicy) error{NonFiniteInput}!bool {
            try requireFiniteField(self.center);
            try requireFiniteField(point);
            const sign = ballContainsSign(4, T, self.center.toArray(), self.radius, point.toArray(), policy);
            return sign == .zero or sign == .negative;
        }

        /// The set of points common to two hyperspheres.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{NonFiniteInput}!Intersection4(T) {
            try requireFiniteField(self.center);
            try requireFiniteField(other.center);
            return ballIntersect(4, T, self.center.toArray(), self.radius, other.center.toArray(), other.radius, policy);
        }

        /// The set of points common to this hypersphere and a segment.
        pub fn intersectSegment(self: Self, segment: SegmentOf(4, T), policy: PredicatePolicy) error{NonFiniteInput}!Intersection4(T) {
            try requireFiniteField(self.center);
            try requireFiniteField(segment);
            return ballIntersectSegment(4, T, self.center.toArray(), self.radius, segment.a.toArray(), segment.b.toArray(), policy);
        }

        /// The set of points common to this hypersphere and a ray.
        pub fn intersectRay(self: Self, ray: RayOf(4, T), policy: PredicatePolicy) error{NonFiniteInput}!Intersection4(T) {
            try requireFiniteField(self.center);
            try requireFiniteField(ray);
            return ballIntersectRay(4, T, self.center.toArray(), self.radius, ray.origin.toArray(), ray.direction.toArray(), policy);
        }

        /// The set of points common to this hypersphere and a hyperplane:
        /// nothing, the one tangent point, or the three-ball boundary the
        /// hyperplane cuts.
        pub fn intersectHyperplane(self: Self, plane: Hyperplane4(T), policy: PredicatePolicy) error{NonFiniteInput}!Intersection4(T) {
            try requireFiniteField(self.center);
            try requireFiniteField(plane);
            return ballIntersectFlat(4, T, self.center.toArray(), self.radius, plane.origin.toArray(), plane.normal.toArray(), policy);
        }

        /// True when the centre and the radius are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Point4(T), self.center) and isFiniteOf(T, self.radius);
        }
    };
}

// ---------------------------------------------------------------------------
// Plane, hyperplane
// ---------------------------------------------------------------------------

/// A three-dimensional `f32` plane.
pub const Plane3f = Plane3(f32);
/// A three-dimensional `f64` plane.
pub const Plane3d = Plane3(f64);
/// A four-dimensional `f32` hyperplane.
pub const Hyperplane4f = Hyperplane4(f32);
/// A four-dimensional `f64` hyperplane.
pub const Hyperplane4d = Hyperplane4(f64);

fn PlaneOf(comptime dim: usize, comptime T: type) type {
    return switch (dim) {
        3 => Plane3(T),
        4 => Hyperplane4(T),
        else => @compileError("zmath flats are three- and four-dimensional"),
    };
}

fn flatBuild(comptime P: type, point: P, normal: VectorOf(fieldCount(P), componentType(P))) error{ NonFiniteInput, Degenerate }!struct { origin: P, normal: VectorOf(fieldCount(P), componentType(P)) } {
    const built = try lineBuild(P, point, normal);
    return .{ .origin = built.origin, .normal = built.direction };
}

fn Plane3(comptime T: type) type {
    return struct {
        /// A point on the plane.
        origin: Point3(T),
        /// The unit normal.
        normal: Vector3(T),

        const Self = @This();

        /// Builds the plane through `point` with `normal`, which is
        /// normalized. A zero normal has no plane to build.
        pub fn init(point: Point3(T), normal: Vector3(T)) error{ NonFiniteInput, Degenerate }!Self {
            const built = try flatBuild(Point3(T), point, normal);
            return .{ .origin = built.origin, .normal = built.normal };
        }

        /// The signed distance from `point` to the plane. It is negative on
        /// the side the normal points away from, zero exactly on the plane,
        /// and positive on the other side.
        pub fn signedDistance(self: Self, point: Point3(T)) error{NonFiniteInput}!T {
            try requireFiniteScalar(T, point.x);
            try requireFiniteScalar(T, point.y);
            try requireFiniteScalar(T, point.z);
            return flatDistanceFast(3, T, point.toArray(), self.origin.toArray(), self.normal.toArray());
        }

        /// The set of points common to two planes: nothing, a whole line when
        /// they are the same line's plane, or the line they meet in.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{NonFiniteInput}!Intersection3(T) {
            try requireFiniteField(self.origin);
            try requireFiniteField(other.origin);
            return flatPair3(T, self.origin.toArray(), self.normal.toArray(), other.origin.toArray(), other.normal.toArray(), policy);
        }

        /// The set of points common to this plane and a segment: nothing, the
        /// one point the plane cuts the segment at, or the whole segment when
        /// it lies in the plane.
        pub fn intersectSegment(self: Self, segment: SegmentOf(3, T), policy: PredicatePolicy) error{NonFiniteInput}!Intersection3(T) {
            try requireFiniteField(self.origin);
            try requireFiniteField(segment);
            return flatIntersectSegment(3, T, self.origin.toArray(), self.normal.toArray(), segment.a.toArray(), segment.b.toArray(), policy);
        }

        /// The set of points common to this plane and a sphere: nothing, the
        /// one tangent point, or the circle.
        pub fn intersectSphere(self: Self, sphere: Sphere3(T), policy: PredicatePolicy) error{NonFiniteInput}!Intersection3(T) {
            try requireFiniteField(self.origin);
            try requireFiniteField(sphere);
            return ballIntersectFlat(3, T, sphere.center.toArray(), sphere.radius, self.origin.toArray(), self.normal.toArray(), policy);
        }

        /// True when the origin and the normal are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Point3(T), self.origin) and isFiniteQuiet(Vector3(T), self.normal);
        }
    };
}

fn Hyperplane4(comptime T: type) type {
    return struct {
        /// A point on the hyperplane.
        origin: Point4(T),
        /// The unit normal.
        normal: Vector4(T),

        const Self = @This();

        /// Builds the hyperplane through `point` with `normal`, which is
        /// normalized. A zero normal has no hyperplane to build.
        pub fn init(point: Point4(T), normal: Vector4(T)) error{ NonFiniteInput, Degenerate }!Self {
            const built = try flatBuild(Point4(T), point, normal);
            return .{ .origin = built.origin, .normal = built.normal };
        }

        /// The signed distance from `point` to the hyperplane.
        pub fn signedDistance(self: Self, point: Point4(T)) error{NonFiniteInput}!T {
            try requireFiniteScalar(T, point.x);
            try requireFiniteScalar(T, point.y);
            try requireFiniteScalar(T, point.z);
            try requireFiniteScalar(T, point.w);
            return flatDistanceFast(4, T, point.toArray(), self.origin.toArray(), self.normal.toArray());
        }

        /// The set of points common to two hyperplanes: nothing, the whole
        /// three-flat when they are one, or the two-flat they meet in.
        pub fn intersect(self: Self, other: Self, policy: PredicatePolicy) error{NonFiniteInput}!Intersection4(T) {
            try requireFiniteField(self.origin);
            try requireFiniteField(other.origin);
            return flatPair4(T, self.origin.toArray(), self.normal.toArray(), other.origin.toArray(), other.normal.toArray(), policy);
        }

        /// The set of points common to this hyperplane and a segment.
        pub fn intersectSegment(self: Self, segment: SegmentOf(4, T), policy: PredicatePolicy) error{NonFiniteInput}!Intersection4(T) {
            try requireFiniteField(self.origin);
            try requireFiniteField(segment);
            return flatIntersectSegment(4, T, self.origin.toArray(), self.normal.toArray(), segment.a.toArray(), segment.b.toArray(), policy);
        }

        /// The set of points common to this hyperplane and a hypersphere:
        /// nothing, the one tangent point, or the three-ball boundary.
        pub fn intersectSphere(self: Self, sphere: Hypersphere4(T), policy: PredicatePolicy) error{NonFiniteInput}!Intersection4(T) {
            try requireFiniteField(self.origin);
            try requireFiniteField(sphere);
            return ballIntersectFlat(4, T, sphere.center.toArray(), sphere.radius, self.origin.toArray(), self.normal.toArray(), policy);
        }

        /// True when the origin and the normal are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(Point4(T), self.origin) and isFiniteQuiet(Vector4(T), self.normal);
        }
    };
}

// ---------------------------------------------------------------------------
// Axis-aligned bounding box
// ---------------------------------------------------------------------------

fn aabbBuild(comptime P: type, minimum: P, maximum: P) error{ NonFiniteInput, InvertedBounds }!P {
    try requireFiniteField(minimum);
    try requireFiniteField(maximum);
    const lows = toArrayOf(P, minimum);
    const highs = toArrayOf(P, maximum);
    inline for (0..comptime fieldCount(P)) |index| {
        if (lows[index] > highs[index]) return error.InvertedBounds;
    }
    return minimum;
}

/// The finiteness gate for a box CONSUMED by an operation, as opposed to the
/// one a constructor ran. It has to cover BOTH corners: the fields are
/// public, so a caller can poison `maximum` after `init` accepted it, and
/// every operation below reads `maximum` in its arithmetic. Validating
/// `minimum` alone let an infinite `maximum` through and produced a finite
/// answer from it, which is exactly the propagation the module comment
/// forbids.
fn requireFiniteBox(comptime P: type, minimum: P, maximum: P) error{NonFiniteInput}!void {
    try requireFiniteField(minimum);
    try requireFiniteField(maximum);
}

/// A two-dimensional `f32` axis-aligned box.
pub const Aabb2f = Aabb2(f32);
/// A two-dimensional `f64` axis-aligned box.
pub const Aabb2d = Aabb2(f64);
/// A three-dimensional `f32` axis-aligned box.
pub const Aabb3f = Aabb3(f32);
/// A three-dimensional `f64` axis-aligned box.
pub const Aabb3d = Aabb3(f64);
/// A four-dimensional `f32` axis-aligned box.
pub const Aabb4f = Aabb4(f32);
/// A four-dimensional `f64` axis-aligned box.
pub const Aabb4d = Aabb4(f64);

fn AabbOf(comptime dim: usize, comptime T: type) type {
    return switch (dim) {
        2 => Aabb2(T),
        3 => Aabb3(T),
        4 => Aabb4(T),
        else => @compileError("zmath geometry dimensions are 2, 3, and 4"),
    };
}

fn aabbContainsOf(comptime P: type, minimum: P, maximum: P, point: P) bool {
    const lows = toArrayOf(P, minimum);
    const highs = toArrayOf(P, maximum);
    const values = toArrayOf(P, point);
    var inside = true;
    inline for (0..comptime fieldCount(P)) |index| {
        inside = inside and values[index] >= lows[index] and values[index] <= highs[index];
    }
    return inside;
}

fn aabbOverlapsOf(comptime P: type, a_minimum: P, a_maximum: P, b_minimum: P, b_maximum: P) bool {
    const a_low = toArrayOf(P, a_minimum);
    const a_high = toArrayOf(P, a_maximum);
    const b_low = toArrayOf(P, b_minimum);
    const b_high = toArrayOf(P, b_maximum);
    var overlaps = true;
    inline for (0..comptime fieldCount(P)) |index| {
        overlaps = overlaps and a_low[index] <= b_high[index] and b_low[index] <= a_high[index];
    }
    return overlaps;
}

/// The product of the extents: the area in two dimensions, the volume in
/// three, and the hypervolume in four. One expression covers all three because
/// the measure a box "defines" is the product of its side lengths in every
/// dimension.
fn aabbMeasureOf(comptime P: type, minimum: P, maximum: P) componentType(P) {
    const lows = toArrayOf(P, minimum);
    const highs = toArrayOf(P, maximum);
    var total = multRaw(componentType(P), highs[0] - lows[0], 1);
    inline for (1..comptime fieldCount(P)) |index| {
        total = multRaw(componentType(P), total, highs[index] - lows[index]);
    }
    return total;
}

fn aabbMeasuresOf(comptime P: type, minimum: P, maximum: P) [fieldCount(P)]componentType(P) {
    const lows = toArrayOf(P, minimum);
    const highs = toArrayOf(P, maximum);
    var out: [fieldCount(P)]componentType(P) = undefined;
    inline for (0..comptime fieldCount(P)) |index| {
        out[index] = highs[index] - lows[index];
    }
    return out;
}

/// The smallest box that holds `minimum`, `maximum`, and `point` on both
/// corners. `fmin` and `fmax` are the accepted C23 kernels, so a point that
/// is already inside leaves both corners exactly as they were.
fn aabbExpandedOf(comptime P: type, minimum: P, maximum: P, point: P) error{NonFiniteInput}!struct { minimum: P, maximum: P } {
    try requireFiniteField(point);
    const lows = toArrayOf(P, minimum);
    const highs = toArrayOf(P, maximum);
    const values = toArrayOf(P, point);
    var new_low: [fieldCount(P)]componentType(P) = undefined;
    var new_high: [fieldCount(P)]componentType(P) = undefined;
    inline for (0..comptime fieldCount(P)) |index| {
        new_low[index] = minOf(componentType(P), lows[index], values[index]);
        new_high[index] = maxOf(componentType(P), highs[index], values[index]);
    }
    return .{ .minimum = fromArrayOf(P, new_low), .maximum = fromArrayOf(P, new_high) };
}

fn aabbCenterOf(comptime P: type, minimum: P, maximum: P) P {
    const lows = toArrayOf(P, minimum);
    const highs = toArrayOf(P, maximum);
    var out: [fieldCount(P)]componentType(P) = undefined;
    inline for (0..comptime fieldCount(P)) |index| {
        out[index] = divRaw(componentType(P), addRaw(componentType(P), lows[index], highs[index]), 2);
    }
    return fromArrayOf(P, out);
}

fn Aabb2(comptime T: type) type {
    const P = Point2(T);
    return struct {
        /// The lower corner, which is not above `maximum` on any axis.
        minimum: P,
        /// The upper corner.
        maximum: P,

        const Self = @This();

        /// Builds the box. Every `minimum` component must be at or below the
        /// matching `maximum` component, and an inverted box is an error
        /// rather than a silently swapped pair.
        pub fn init(minimum: P, maximum: P) error{ NonFiniteInput, InvertedBounds }!Self {
            _ = try aabbBuild(P, minimum, maximum);
            return .{ .minimum = minimum, .maximum = maximum };
        }

        /// True when `point` is inside or ON the box. The boundary is
        /// INCLUSIVE: a component exactly on a face is inside, because the
        /// face is part of the box and the range comparison rounds nothing.
        /// This is why `contains` takes no policy: there is no rounding here
        /// for a policy to change.
        pub fn contains(self: Self, point: P) error{NonFiniteInput}!bool {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteField(point);
            return aabbContainsOf(P, self.minimum, self.maximum, point);
        }

        /// The smallest box that holds both this one and `point`, on both
        /// corners.
        pub fn expand(self: Self, point: P) error{NonFiniteInput}!Self {
            try requireFiniteBox(P, self.minimum, self.maximum);
            const grown = try aabbExpandedOf(P, self.minimum, self.maximum, point);
            return .{ .minimum = grown.minimum, .maximum = grown.maximum };
        }

        /// True when the two boxes have a common point. Touching counts:
        /// two boxes that share only a face do have that face in common. The
        /// comparison rounds nothing, so this takes no policy.
        pub fn intersects(self: Self, other: Self) error{NonFiniteInput}!bool {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteBox(P, other.minimum, other.maximum);
            return aabbOverlapsOf(P, self.minimum, self.maximum, other.minimum, other.maximum);
        }

        /// The midpoint of the box.
        pub fn center(self: Self) error{NonFiniteInput}!P {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbCenterOf(P, self.minimum, self.maximum);
        }

        /// The size of the box along each axis, that is `maximum - minimum`.
        /// This is the full extent, not the half extent from the centre.
        pub fn extent(self: Self) error{NonFiniteInput}![2]T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbMeasuresOf(P, self.minimum, self.maximum);
        }

        /// The area of the box.
        pub fn area(self: Self) error{NonFiniteInput}!T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbMeasureOf(P, self.minimum, self.maximum);
        }

        /// The set of points common to two boxes: `none`, a `box` that is
        /// neither operand, the `contained` operand, or `coincident`.
        pub fn intersect(self: Self, other: Self) error{NonFiniteInput}!Intersection2(T) {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteBox(P, other.minimum, other.maximum);
            return aabbIntersect(2, T, self.minimum, self.maximum, other.minimum, other.maximum);
        }

        /// True when both corners are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(P, self.minimum) and isFiniteQuiet(P, self.maximum);
        }
    };
}

fn multRaw(comptime T: type, a: T, b: T) T {
    @setFloatMode(.strict);
    return a * b;
}

fn Aabb3(comptime T: type) type {
    const P = Point3(T);
    return struct {
        /// The lower corner, which is not above `maximum` on any axis.
        minimum: P,
        /// The upper corner.
        maximum: P,

        const Self = @This();

        /// Builds the box, rejecting an inverted one.
        pub fn init(minimum: P, maximum: P) error{ NonFiniteInput, InvertedBounds }!Self {
            _ = try aabbBuild(P, minimum, maximum);
            return .{ .minimum = minimum, .maximum = maximum };
        }

        /// True when `point` is inside or ON the box. The boundary is
        /// INCLUSIVE.
        pub fn contains(self: Self, point: P) error{NonFiniteInput}!bool {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteField(point);
            return aabbContainsOf(P, self.minimum, self.maximum, point);
        }

        /// The smallest box that holds both this one and `point`, on both
        /// corners.
        pub fn expand(self: Self, point: P) error{NonFiniteInput}!Self {
            try requireFiniteBox(P, self.minimum, self.maximum);
            const grown = try aabbExpandedOf(P, self.minimum, self.maximum, point);
            return .{ .minimum = grown.minimum, .maximum = grown.maximum };
        }

        /// True when the two boxes have a common point. Touching counts.
        pub fn intersects(self: Self, other: Self) error{NonFiniteInput}!bool {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteBox(P, other.minimum, other.maximum);
            return aabbOverlapsOf(P, self.minimum, self.maximum, other.minimum, other.maximum);
        }

        /// The midpoint of the box.
        pub fn center(self: Self) error{NonFiniteInput}!P {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbCenterOf(P, self.minimum, self.maximum);
        }

        /// The size of the box along each axis, that is `maximum - minimum`.
        pub fn extent(self: Self) error{NonFiniteInput}![3]T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbMeasuresOf(P, self.minimum, self.maximum);
        }

        /// The volume of the box.
        pub fn volume(self: Self) error{NonFiniteInput}!T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbMeasureOf(P, self.minimum, self.maximum);
        }

        /// Twice the area of the box, `2 * (xy + xz + yz)`. The factor of
        /// two keeps the expression a product of the extents, so it is exact
        /// wherever the extents are; a caller who wants the area divides.
        pub fn doubledSurfaceArea(self: Self) error{NonFiniteInput}!T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            const sides = aabbMeasuresOf(P, self.minimum, self.maximum);
            return scaleRaw(T, addRaw(T, addRaw(T, multRaw(T, sides[0], sides[1]), multRaw(T, sides[0], sides[2])), multRaw(T, sides[1], sides[2])), 2);
        }

        /// The set of points common to two boxes.
        pub fn intersect(self: Self, other: Self) error{NonFiniteInput}!Intersection3(T) {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteBox(P, other.minimum, other.maximum);
            return aabbIntersect(3, T, self.minimum, self.maximum, other.minimum, other.maximum);
        }

        /// True when both corners are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(P, self.minimum) and isFiniteQuiet(P, self.maximum);
        }
    };
}

fn Aabb4(comptime T: type) type {
    const P = Point4(T);
    return struct {
        /// The lower corner, which is not above `maximum` on any axis.
        minimum: P,
        /// The upper corner.
        maximum: P,

        const Self = @This();

        /// Builds the box, rejecting an inverted one.
        pub fn init(minimum: P, maximum: P) error{ NonFiniteInput, InvertedBounds }!Self {
            _ = try aabbBuild(P, minimum, maximum);
            return .{ .minimum = minimum, .maximum = maximum };
        }

        /// True when `point` is inside or ON the box. The boundary is
        /// INCLUSIVE.
        pub fn contains(self: Self, point: P) error{NonFiniteInput}!bool {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteField(point);
            return aabbContainsOf(P, self.minimum, self.maximum, point);
        }

        /// The smallest box that holds both this one and `point`, on both
        /// corners.
        pub fn expand(self: Self, point: P) error{NonFiniteInput}!Self {
            try requireFiniteBox(P, self.minimum, self.maximum);
            const grown = try aabbExpandedOf(P, self.minimum, self.maximum, point);
            return .{ .minimum = grown.minimum, .maximum = grown.maximum };
        }

        /// True when the two boxes have a common point. Touching counts.
        pub fn intersects(self: Self, other: Self) error{NonFiniteInput}!bool {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteBox(P, other.minimum, other.maximum);
            return aabbOverlapsOf(P, self.minimum, self.maximum, other.minimum, other.maximum);
        }

        /// The midpoint of the box.
        pub fn center(self: Self) error{NonFiniteInput}!P {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbCenterOf(P, self.minimum, self.maximum);
        }

        /// The size of the box along each axis, that is `maximum - minimum`.
        pub fn extent(self: Self) error{NonFiniteInput}![4]T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbMeasuresOf(P, self.minimum, self.maximum);
        }

        /// The four-dimensional hypervolume of the box.
        pub fn hyperVolume(self: Self) error{NonFiniteInput}!T {
            try requireFiniteBox(P, self.minimum, self.maximum);
            return aabbMeasureOf(P, self.minimum, self.maximum);
        }

        /// The set of points common to two boxes.
        pub fn intersect(self: Self, other: Self) error{NonFiniteInput}!Intersection4(T) {
            try requireFiniteBox(P, self.minimum, self.maximum);
            try requireFiniteBox(P, other.minimum, other.maximum);
            return aabbIntersect(4, T, self.minimum, self.maximum, other.minimum, other.maximum);
        }

        /// True when both corners are finite.
        pub fn isFinite(self: Self) bool {
            return isFiniteQuiet(P, self.minimum) and isFiniteQuiet(P, self.maximum);
        }
    };
}

fn aabbIntersect(comptime dim: usize, comptime T: type, a_minimum: PointOf(dim, T), a_maximum: PointOf(dim, T), b_minimum: PointOf(dim, T), b_maximum: PointOf(dim, T)) IntersectionOf(dim, T) {
    const P = PointOf(dim, T);
    const a_low = toArrayOf(P, a_minimum);
    const a_high = toArrayOf(P, a_maximum);
    const b_low = toArrayOf(P, b_minimum);
    const b_high = toArrayOf(P, b_maximum);
    var low: [dim]T = undefined;
    var high: [dim]T = undefined;
    var disjoint = false;
    inline for (0..dim) |index| {
        low[index] = maxOf(T, a_low[index], b_low[index]);
        high[index] = minOf(T, a_high[index], b_high[index]);
        disjoint = disjoint or low[index] > high[index];
    }
    if (disjoint) return .none;
    var equal_a = true;
    var equal_b = true;
    var equal_both = true;
    inline for (0..dim) |index| {
        equal_a = equal_a and low[index] == a_low[index] and high[index] == a_high[index];
        equal_b = equal_b and low[index] == b_low[index] and high[index] == b_high[index];
        equal_both = equal_both and low[index] == a_low[index] and high[index] == a_high[index] and
            low[index] == b_low[index] and high[index] == b_high[index];
    }
    if (equal_both) return .coincident;
    if (equal_a) return .{ .contained = .first };
    if (equal_b) return .{ .contained = .second };
    return .{ .box = .{ .minimum = fromArrayOf(P, low), .maximum = fromArrayOf(P, high) } };
}

// ---------------------------------------------------------------------------
// The affine intersection engine
//
// Every pair of lines, rays, and segments reduces to the same question, in
// any dimension: do two coplanar objects meet, and where. In two dimensions
// the question is asked directly. In three and four dimensions it is asked
// after an exact coplanarity gate and a projection into a plane, because two
// objects in a higher dimension meet only if they share a plane first.
// ---------------------------------------------------------------------------

const AffineKind = enum { none, point, segment, line, coincident, contained };

fn AffineResultOf(comptime T: type) type {
    return struct {
        kind: AffineKind,
        /// The near parameter along the first object, and for `.point` the
        /// only one.
        first: T,
        /// The far parameter along the first object, used by `.segment`.
        second: T,
        /// Which operand the whole intersection is, for `.contained`.
        side: ContainedSide,
    };
}

fn noneAffine(comptime T: type) AffineResultOf(T) {
    return .{ .kind = .none, .first = 0, .second = 0, .side = .first };
}

/// One affine object as the engine sees it: the two points that define its
/// line, the difference, and the closed range of parameters along that
/// difference. A segment is `[0, 1]`, a ray is `[0, +infinity)`, and a line
/// is the whole real line, so all three share one code path.
fn AffineObjectOf(comptime T: type) type {
    return struct {
        start: [2]T,
        finish: [2]T,
        direction: [2]T,
        lower: T,
        upper: T,
    };
}

fn affineObjectOf(comptime T: type, start: [2]T, finish: [2]T, lower: T, upper: T) AffineObjectOf(T) {
    return .{
        .start = start,
        .finish = finish,
        .direction = arraySub(2, T, finish, start),
        .lower = lower,
        .upper = upper,
    };
}

fn inRangeOf(comptime T: type, object: AffineObjectOf(T), t: T) bool {
    return t >= object.lower and t <= object.upper;
}

/// The sign of `cross(direction, point - start)`, degree two.
fn lineSideSign2(comptime T: type, object: AffineObjectOf(T), point: [2]T, policy: PredicatePolicy) Sign {
    const offset = arraySub(2, T, point, object.start);
    return switch (policy) {
        .fast => signOfScalar(T, subRaw(T, scaleRaw(T, object.direction[0], offset[1]), scaleRaw(T, object.direction[1], offset[0]))),
        .adaptive_exact => exactDet2(widen(T, object.direction[0]), widen(T, offset[1]), widen(T, object.direction[1]), widen(T, offset[0])),
    };
}

/// True when `point` lies on the supporting line of `object` and its
/// parameter is inside the object's range. The collinearity test is the
/// exact predicate; the range test is a comparison of computed values.
fn pointOnAffine2(comptime T: type, object: AffineObjectOf(T), point: [2]T, policy: PredicatePolicy) bool {
    if (lineSideSign2(T, object, point, policy) != .zero) return false;
    const denominator = arrayDot(2, T, object.direction, object.direction);
    if (denominator == 0) return point[0] == object.start[0] and point[1] == object.start[1];
    const offset = arraySub(2, T, point, object.start);
    const parameter = divRaw(T, arrayDot(2, T, offset, object.direction), denominator);
    return inRangeOf(T, object, parameter);
}

/// The overlap of two collinear objects, projected onto the axis along which
/// the first object's direction is largest. That axis is the one the shared
/// direction is largest in, so the projection is monotone on both objects
/// and the one-dimensional overlap is the whole answer.
///
/// The projection goes through each object's PARAMETER RANGE, not through its
/// two defining points. A ray's two defining points are one unit apart, so
/// reading its extent off them would clip a ray to a segment, and a ray
/// overlapping a segment on the same line is the case a caller reaches first.
/// The range bounds are what make an infinite end an infinite coordinate.
fn collinearOverlap(comptime T: type, a: AffineObjectOf(T), b: AffineObjectOf(T)) AffineResultOf(T) {
    var axis: usize = 0;
    var widest: T = absOf(T, a.direction[0]);
    inline for (1..2) |index| {
        const candidate = absOf(T, a.direction[index]);
        if (candidate > widest) {
            widest = candidate;
            axis = index;
        }
    }
    const a_at_lower = addRaw(T, a.start[axis], scaleRaw(T, a.direction[axis], a.lower));
    const a_at_upper = addRaw(T, a.start[axis], scaleRaw(T, a.direction[axis], a.upper));
    const b_at_lower = addRaw(T, b.start[axis], scaleRaw(T, b.direction[axis], b.lower));
    const b_at_upper = addRaw(T, b.start[axis], scaleRaw(T, b.direction[axis], b.upper));
    const a_low = minOf(T, a_at_lower, a_at_upper);
    const a_high = maxOf(T, a_at_lower, a_at_upper);
    const b_low = minOf(T, b_at_lower, b_at_upper);
    const b_high = maxOf(T, b_at_lower, b_at_upper);
    const low = maxOf(T, a_low, b_low);
    const high = minOf(T, a_high, b_high);
    if (low > high) return noneAffine(T);
    const span = a.direction[axis];
    const first = divRaw(T, subRaw(T, low, a.start[axis]), span);
    const second = divRaw(T, subRaw(T, high, a.start[axis]), span);
    if (low == high) return .{ .kind = .point, .first = first, .second = first, .side = .first };
    const whole_a = low == a_low and high == a_high;
    const whole_b = low == b_low and high == b_high;
    if (whole_a and whole_b) return .{ .kind = .coincident, .first = first, .second = second, .side = .first };
    if (whole_a) return .{ .kind = .contained, .first = first, .second = second, .side = .first };
    if (whole_b) return .{ .kind = .contained, .first = first, .second = second, .side = .second };
    return .{ .kind = .segment, .first = first, .second = second, .side = .first };
}

/// The intersection of two coplanar two-dimensional affine objects.
fn affine2Intersect(comptime T: type, a: AffineObjectOf(T), b: AffineObjectOf(T), policy: PredicatePolicy) AffineResultOf(T) {
    // A degenerate segment is a single point, and then the whole question is
    // whether that point is on the other object.
    if (arrayLength(2, T, a.direction) == 0) {
        if (pointOnAffine2(T, b, a.start, policy)) {
            return .{ .kind = .point, .first = 0, .second = 0, .side = .first };
        }
        return noneAffine(T);
    }
    if (arrayLength(2, T, b.direction) == 0) {
        if (pointOnAffine2(T, a, b.start, policy)) {
            return .{ .kind = .point, .first = 1, .second = 1, .side = .second };
        }
        return noneAffine(T);
    }

    const a_first = lineSideSign2(T, a, b.start, policy);
    const a_second = lineSideSign2(T, a, b.finish, policy);
    const b_first = lineSideSign2(T, b, a.start, policy);
    const b_second = lineSideSign2(T, b, a.finish, policy);
    const a_collinear = a_first == .zero and a_second == .zero;
    const b_collinear = b_first == .zero and b_second == .zero;
    if (a_collinear or b_collinear) return collinearOverlap(T, a, b);

    // The supporting lines are not parallel, so they meet at exactly one
    // parameter pair, found by Cramer's rule on the two-by-two system.
    const denominator = subRaw(T, scaleRaw(T, a.direction[0], b.direction[1]), scaleRaw(T, a.direction[1], b.direction[0]));
    if (denominator == 0) return noneAffine(T);
    const gap = arraySub(2, T, b.start, a.start);
    const first = divRaw(T, subRaw(T, scaleRaw(T, gap[0], b.direction[1]), scaleRaw(T, gap[1], b.direction[0])), denominator);
    const second = divRaw(T, subRaw(T, scaleRaw(T, gap[0], a.direction[1]), scaleRaw(T, gap[1], a.direction[0])), denominator);
    if (inRangeOf(T, a, first) and inRangeOf(T, b, second)) {
        return .{ .kind = .point, .first = first, .second = second, .side = .first };
    }

    // The crossing is off the end of at least one object. The only place they
    // can still meet is a boundary of one of them, so each boundary is tested
    // against the other object.
    if (first < a.lower) {
        const candidate = pointAtArray(2, T, a.start, a.direction, a.lower);
        if (pointOnAffine2(T, b, candidate, policy)) {
            return .{ .kind = .point, .first = a.lower, .second = a.lower, .side = .first };
        }
    }
    if (first > a.upper) {
        const candidate = pointAtArray(2, T, a.start, a.direction, a.upper);
        if (pointOnAffine2(T, b, candidate, policy)) {
            return .{ .kind = .point, .first = a.upper, .second = a.upper, .side = .first };
        }
    }
    if (second < b.lower) {
        const candidate = pointAtArray(2, T, b.start, b.direction, b.lower);
        if (pointOnAffine2(T, a, candidate, policy)) {
            return .{ .kind = .point, .first = 0, .second = b.lower, .side = .second };
        }
    }
    if (second > b.upper) {
        const candidate = pointAtArray(2, T, b.start, b.direction, b.upper);
        if (pointOnAffine2(T, a, candidate, policy)) {
            return .{ .kind = .point, .first = 1, .second = b.upper, .side = .second };
        }
    }
    return noneAffine(T);
}

/// The sign of the three-by-four minor over the three given rows of the
/// matrix whose columns are `u`, `v`, and `w`. Four of these vanish exactly
/// when the three vectors are linearly dependent, which is the
/// four-dimensional coplanarity test, and each is degree two.
fn minor3x4Sign(comptime T: type, u: [4]T, v: [4]T, w: [4]T, rows: [3]usize, policy: PredicatePolicy) Sign {
    return switch (policy) {
        .fast => signOfScalar(T, fastMinor3x4(T, u, v, w, rows)),
        .adaptive_exact => exactMinor3x4(T, u, v, w, rows),
    };
}

/// The sign of the three-dimensional signed volume of the tetrahedron
/// `a`, `b`, `c`, `d`, which is zero exactly when the four points are
/// coplanar. That coplanarity is the gate the three-dimensional affine
/// intersections open on.
fn orient3Sign(comptime T: type, a: [3]T, b: [3]T, c: [3]T, d: [3]T, policy: PredicatePolicy) Sign {
    return switch (policy) {
        .fast => signOfScalar(T, crossDotFast(T, a, b, c, d)),
        .adaptive_exact => exactOrient3(T, a, b, c, d),
    };
}

/// A coordinate axis whose projection out of the whole basis is largest, so
/// the residual is far from zero and the normalization is well conditioned.
/// The axis is chosen against EVERY basis vector, not the first one alone:
/// an axis that is least aligned with the first vector can still lie wholly
/// inside the span of the others, and a seed inside the span projects to
/// exactly zero, which normalizes to the zero vector rather than a direction.
/// `basis` must be shorter than the dimension, so some axis is outside it.
fn seedOutsideBasis(comptime dim: usize, comptime T: type, basis: []const [dim]T) [dim]T {
    var seed: [dim]T = undefined;
    inline for (0..dim) |index| {
        seed[index] = 0;
    }
    var best: usize = 0;
    var best_weight: T = std.math.inf(T);
    inline for (0..dim) |index| {
        // The squared projection of the axis onto the span, so the SMALLEST
        // weight names the axis that lies most outside it.
        var weight: T = 0;
        for (basis) |vector| {
            weight = addRaw(T, weight, scaleRaw(T, vector[index], vector[index]));
        }
        if (weight < best_weight) {
            best_weight = weight;
            best = index;
        }
    }
    seed[best] = 1;
    return orthogonalized(dim, T, seed, basis);
}

///
/// The intersection of two coplanar affine objects in three or four
/// dimensions. The coplanarity gate is the exact predicate on the ORIGINAL
/// coordinates; the classification after the projection is the exact
/// two-dimensional predicate on the PROJECTED coordinates, which is what the
/// module comment promises and no more.
fn affineIntersectHigh(
    comptime dim: usize,
    comptime T: type,
    a_start: [dim]T,
    a_finish: [dim]T,
    a_lower: T,
    a_upper: T,
    b_start: [dim]T,
    b_finish: [dim]T,
    b_lower: T,
    b_upper: T,
    policy: PredicatePolicy,
) AffineResultOf(T) {
    const u = arraySub(dim, T, a_finish, a_start);
    const v = arraySub(dim, T, b_finish, b_start);
    const w = arraySub(dim, T, b_start, a_start);

    if (dim == 3) {
        if (orient3Sign(T, a_start, a_finish, b_start, b_finish, policy) != .zero) return noneAffine(T);
    } else {
        if (minor3x4Sign(T, u, v, w, .{ 0, 1, 2 }, policy) != .zero) return noneAffine(T);
        if (minor3x4Sign(T, u, v, w, .{ 0, 1, 3 }, policy) != .zero) return noneAffine(T);
        if (minor3x4Sign(T, u, v, w, .{ 0, 2, 3 }, policy) != .zero) return noneAffine(T);
        if (minor3x4Sign(T, u, v, w, .{ 1, 2, 3 }, policy) != .zero) return noneAffine(T);
    }

    const first = arrayNormalized(dim, T, u);
    if (arrayLength(dim, T, first) == 0) return noneAffine(T);
    const residual = arraySub(dim, T, v, arrayScale(dim, T, first, arrayDot(dim, T, v, first)));
    const second = if (arrayLength(dim, T, residual) > 0)
        arrayNormalized(dim, T, residual)
    else
        anyPerpendicular(dim, T, first);

    const first_object = affineObjectOf(T, [2]T{
        arrayDot(dim, T, a_start, first),
        arrayDot(dim, T, a_start, second),
    }, [2]T{
        arrayDot(dim, T, a_finish, first),
        arrayDot(dim, T, a_finish, second),
    }, a_lower, a_upper);
    const second_object = affineObjectOf(T, [2]T{
        arrayDot(dim, T, b_start, first),
        arrayDot(dim, T, b_start, second),
    }, [2]T{
        arrayDot(dim, T, b_finish, first),
        arrayDot(dim, T, b_finish, second),
    }, b_lower, b_upper);
    return affine2Intersect(T, first_object, second_object, policy);
}

/// The concrete union for one classification. `a_lower` and `a_upper` say how
/// far the first object runs, and a `coincident` result over an unbounded
/// range is a whole LINE, not two objects that happen to match: naming the
/// span tells the caller what it got.
fn affineToIntersection(
    comptime dim: usize,
    comptime T: type,
    a_start: [dim]T,
    a_direction: [dim]T,
    a_lower: T,
    a_upper: T,
    result: AffineResultOf(T),
) IntersectionOf(dim, T) {
    const P = PointOf(dim, T);
    switch (result.kind) {
        .none => return .none,
        .point => return .{ .point = fromArrayOf(P, pointAtArray(dim, T, a_start, a_direction, result.first)) },
        .segment => return .{ .segment = .{
            .a = fromArrayOf(P, pointAtArray(dim, T, a_start, a_direction, result.first)),
            .b = fromArrayOf(P, pointAtArray(dim, T, a_start, a_direction, result.second)),
        } },
        .line => {
            if (arrayLength(dim, T, a_direction) == 0) return .none;
            return .{ .line = .{
                .origin = fromArrayOf(P, a_start),
                .direction = fromArrayOf(VectorOf(dim, T), arrayNormalized(dim, T, a_direction)),
            } };
        },
        .coincident => {
            if (!(std.math.isFinite(a_lower) and std.math.isFinite(a_upper))) {
                if (arrayLength(dim, T, a_direction) == 0) return .none;
                return .{ .line = .{
                    .origin = fromArrayOf(P, a_start),
                    .direction = fromArrayOf(VectorOf(dim, T), arrayNormalized(dim, T, a_direction)),
                } };
            }
            return .coincident;
        },
        .contained => return .{ .contained = result.side },
    }
}

fn affineIntersectDispatch(
    comptime dim: usize,
    comptime T: type,
    a_start: [dim]T,
    a_finish: [dim]T,
    a_lower: T,
    a_upper: T,
    b_start: [dim]T,
    b_finish: [dim]T,
    b_lower: T,
    b_upper: T,
    policy: PredicatePolicy,
) IntersectionOf(dim, T) {
    const direction = arraySub(dim, T, a_finish, a_start);
    const result = if (dim == 2)
        affine2Intersect(T, affineObjectOf(T, a_start, a_finish, a_lower, a_upper), affineObjectOf(T, b_start, b_finish, b_lower, b_upper), policy)
    else
        affineIntersectHigh(dim, T, a_start, a_finish, a_lower, a_upper, b_start, b_finish, b_lower, b_upper, policy);
    return affineToIntersection(dim, T, a_start, direction, a_lower, a_upper, result);
}

// ---------------------------------------------------------------------------
// Ball intersections
// ---------------------------------------------------------------------------

fn ballIntersect(comptime dim: usize, comptime T: type, first_center: [dim]T, first_radius: T, second_center: [dim]T, second_radius: T, policy: PredicatePolicy) IntersectionOf(dim, T) {
    const P = PointOf(dim, T);
    const separation = arraySub(dim, T, second_center, first_center);
    // The two comparisons are the two tangency conditions, and their SIGNS run
    // the way the geometry does: the centres are farther apart than the radii
    // add when `outer` is positive, and closer than they differ when `inner`
    // is negative. Getting either sign backwards turns every classification
    // into its opposite, which is the mistake the tests below exist to catch.
    if (arrayLength(dim, T, separation) == 0 and first_radius == second_radius) return .coincident;
    const outer = squaredNormMinusSquareSign(dim, T, separation, addRaw(T, first_radius, second_radius), policy);
    const inner = squaredNormMinusSquareSign(dim, T, separation, subRaw(T, first_radius, second_radius), policy);
    if (outer == .positive) return .none;
    if (outer == .zero or inner == .zero) {
        // Tangency: one point on the line of centres, at distance `first`
        // from the first centre along the separation.
        const unit = arrayNormalized(dim, T, separation);
        const touch = arrayAdd(dim, T, first_center, arrayScale(dim, T, unit, first_radius));
        return .{ .point = fromArrayOf(P, touch) };
    }
    if (inner == .negative) return .{ .contained = if (first_radius >= second_radius) .first else .second };
    return .{ .two_points = ballCrossingPoints(dim, T, first_center, first_radius, second_center, second_radius) };
}

/// The two points where two balls cut each other. The geometry is fixed: the
/// chord lies in the plane through both centres, perpendicular to the line of
/// centres, at the distance the law of cosines puts it.
fn ballCrossingPoints(comptime dim: usize, comptime T: type, first_center: [dim]T, first_radius: T, second_center: [dim]T, second_radius: T) [2]PointOf(dim, T) {
    const P = PointOf(dim, T);
    const separation = arraySub(dim, T, second_center, first_center);
    const distance = arrayLength(dim, T, separation);
    const unit = arrayScale(dim, T, separation, divRaw(T, 1, distance));
    const along = divRaw(T, addRaw(T, subRaw(T, scaleRaw(T, first_radius, first_radius), scaleRaw(T, second_radius, second_radius)), scaleRaw(T, distance, distance)), scaleRaw(T, 2, distance));
    const height_squared = subRaw(T, scaleRaw(T, first_radius, first_radius), scaleRaw(T, along, along));
    const height = sqrtOf(T, absOf(T, height_squared));
    const base = arrayAdd(dim, T, first_center, arrayScale(dim, T, unit, along));
    const side = arrayScale(dim, T, anyPerpendicular(dim, T, unit), height);
    return .{
        fromArrayOf(P, arrayAdd(dim, T, base, side)),
        fromArrayOf(P, arraySub(dim, T, base, side)),
    };
}

fn ballIntersectSegment(comptime dim: usize, comptime T: type, center: [dim]T, radius: T, a: [dim]T, b: [dim]T, policy: PredicatePolicy) IntersectionOf(dim, T) {
    const P = PointOf(dim, T);
    const from_a = arraySub(dim, T, a, center);
    const from_b = arraySub(dim, T, b, center);
    const start_inside = squaredNormMinusSquareSign(dim, T, from_a, radius, policy);
    const end_inside = squaredNormMinusSquareSign(dim, T, from_b, radius, policy);
    // The segment is the second operand of `intersectSegment`, and it is the
    // whole of the common set, so the side is `second`.
    if (start_inside != .positive and end_inside != .positive) return .{ .contained = .second };
    const span = arrayLength(dim, T, arraySub(dim, T, b, a));
    if (span == 0) return .none;
    const unit = arrayScale(dim, T, arraySub(dim, T, b, a), divRaw(T, 1, span));
    const along = arrayDot(dim, T, from_a, unit);
    const discriminant = rayBallDiscriminantSign(dim, T, from_a, along, radius, policy);
    if (discriminant == .negative) return .none;
    const root = sqrtOf(T, absOf(T, subRaw(T, scaleRaw(T, along, along), subRaw(T, arrayDot(dim, T, from_a, from_a), scaleRaw(T, radius, radius)))));
    // The two roots of the ray quadratic are `-along -+ root`. The smaller one
    // enters the ball first, so it is the `low` end of the overlap.
    const near = divRaw(T, subRaw(T, negRaw(T, along), root), span);
    const far = divRaw(T, addRaw(T, negRaw(T, along), root), span);
    const low = maxOf(T, near, 0);
    const high = minOf(T, far, 1);
    if (low > high) return .none;
    if (low == high) return .{ .point = fromArrayOf(P, pointAtArray(dim, T, a, arraySub(dim, T, b, a), low)) };
    return .{ .segment = .{
        .a = fromArrayOf(P, pointAtArray(dim, T, a, arraySub(dim, T, b, a), low)),
        .b = fromArrayOf(P, pointAtArray(dim, T, a, arraySub(dim, T, b, a), high)),
    } };
}

/// The common set of a ray and a ball. The ray quadratic is
/// `t^2 + 2 * along * t + (dot(offset, offset) - r * r)` with the direction
/// already unit, so the roots are `-along -+ root`.
///
/// The three shapes follow from where the ORIGIN sits, because the ball is
/// the CLOSED solid ball and the ray is `[0, +infinity)`:
///
/// - `near >= 0`: the origin is outside, so the common set is the two
///   separate points where the ray enters and leaves, or the single tangent
///   point when they coincide.
/// - `near < 0`: the origin is inside, so the common set is the whole run
///   from the origin to `far`, which is a SEGMENT and not a point. The
///   origin itself is a common point, so reporting only the exit would drop
///   it, and reporting only the origin would drop the run between them.
///   `far == 0` is the boundary case, where that run has zero length.
fn ballIntersectRay(comptime dim: usize, comptime T: type, center: [dim]T, radius: T, origin: [dim]T, direction: [dim]T, policy: PredicatePolicy) IntersectionOf(dim, T) {
    const P = PointOf(dim, T);
    const offset = arraySub(dim, T, origin, center);
    const along = arrayDot(dim, T, offset, direction);
    const discriminant = rayBallDiscriminantSign(dim, T, offset, along, radius, policy);
    if (discriminant == .negative) return .none;
    const root = sqrtOf(T, absOf(T, subRaw(T, scaleRaw(T, along, along), subRaw(T, arrayDot(dim, T, offset, offset), scaleRaw(T, radius, radius)))));
    // Both roots carry `-along`. The smaller is the near root, so it is
    // `-along - root`; the larger is `-along + root`.
    const near = negRaw(T, addRaw(T, along, root));
    const far = addRaw(T, negRaw(T, along), root);
    // The whole chord lies behind the origin, which no point of the ray
    // reaches.
    if (far < 0) return .none;
    if (near >= 0) {
        if (near == far) return .{ .point = fromArrayOf(P, pointAtArray(dim, T, origin, direction, near)) };
        return .{ .two_points = .{
            fromArrayOf(P, pointAtArray(dim, T, origin, direction, near)),
            fromArrayOf(P, pointAtArray(dim, T, origin, direction, far)),
        } };
    }
    // The origin is inside the ball, so the run starts at the origin. A run
    // of zero length is a single point, which is the case of a ray that
    // starts on the surface and points away from the ball.
    if (far == 0) return .{ .point = fromArrayOf(P, origin) };
    return .{ .segment = .{
        .a = fromArrayOf(P, origin),
        .b = fromArrayOf(P, pointAtArray(dim, T, origin, direction, far)),
    } };
}

fn ballIntersectFlat(comptime dim: usize, comptime T: type, center: [dim]T, radius: T, flat_point: [dim]T, flat_normal: [dim]T, policy: PredicatePolicy) IntersectionOf(dim, T) {
    const distance = arrayDot(dim, T, arraySub(dim, T, center, flat_point), flat_normal);
    const sign = flatBallSign(T, distance, radius, policy);
    if (sign == .positive) return .none;
    const foot = arraySub(dim, T, center, arrayScale(dim, T, flat_normal, distance));
    if (sign == .zero) return .{ .point = fromArrayOf(PointOf(dim, T), foot) };
    const cut_radius = sqrtOf(T, absOf(T, subRaw(T, scaleRaw(T, radius, radius), scaleRaw(T, distance, distance))));
    const ball = SphereOf(dim, T){
        .center = fromArrayOf(PointOf(dim, T), foot),
        .radius = cut_radius,
    };
    return if (dim == 3) .{ .circle = ball } else .{ .hypersphere = ball };
}

// ---------------------------------------------------------------------------
// Flat intersections
// ---------------------------------------------------------------------------

fn flatIntersectSegment(comptime dim: usize, comptime T: type, flat_point: [dim]T, flat_normal: [dim]T, a: [dim]T, b: [dim]T, policy: PredicatePolicy) IntersectionOf(dim, T) {
    const P = PointOf(dim, T);
    const first = flatSideSign(dim, T, a, flat_point, flat_normal, policy);
    const second = flatSideSign(dim, T, b, flat_point, flat_normal, policy);
    if (first == .zero and second == .zero) return .{ .contained = .second };
    if (first == .zero) return .{ .point = fromArrayOf(P, a) };
    if (second == .zero) return .{ .point = fromArrayOf(P, b) };
    if ((first == .negative) == (second == .negative)) return .none;
    // The flat cuts the segment strictly inside, so the parameter is the
    // share of the two signed distances.
    const left = flatDistanceFast(dim, T, a, flat_point, flat_normal);
    const right = flatDistanceFast(dim, T, b, flat_point, flat_normal);
    const t = divRaw(T, left, subRaw(T, left, right));
    return .{ .point = fromArrayOf(P, pointAtArray(dim, T, a, arraySub(dim, T, b, a), t)) };
}

/// True when the two normals are parallel, decided on the two-by-two minors
/// of the matrix whose ROWS are the normals. A minor uses one component of
/// each normal against a DIFFERENT component of the other, which is what
/// `minor2Sign(first[i], second[j], first[j], second[i])` reads. Passing the
/// same index twice on each normal instead tests a relation between two
/// components of one normal and two of the other, which vanishes for every
/// coordinate-aligned pair and would call the two axis normals parallel.
fn normalsParallel3(comptime T: type, first: [3]T, second: [3]T, policy: PredicatePolicy) bool {
    return minor2Sign(T, first[0], second[1], first[1], second[0], policy) == .zero and
        minor2Sign(T, first[1], second[2], first[2], second[1], policy) == .zero and
        minor2Sign(T, first[2], second[0], first[0], second[2], policy) == .zero;
}

/// The four-dimensional gate, over all six index pairs. The same argument
/// order as the three-dimensional one, and all six minors vanish exactly when
/// the two rows are linearly dependent.
fn normalsParallel4(comptime T: type, first: [4]T, second: [4]T, policy: PredicatePolicy) bool {
    var parallel = true;
    inline for (0..4) |row| {
        inline for (row + 1..4) |column| {
            if (minor2Sign(T, first[row], second[column], first[column], second[row], policy) != .zero) {
                parallel = false;
            }
        }
    }
    return parallel;
}

fn minor2Sign(comptime T: type, a: T, b: T, c: T, d: T, policy: PredicatePolicy) Sign {
    return switch (policy) {
        .fast => signOfScalar(T, subRaw(T, scaleRaw(T, a, b), scaleRaw(T, c, d))),
        .adaptive_exact => exactDet2(widen(T, a), widen(T, b), widen(T, c), widen(T, d)),
    };
}

/// The intersection of two planes in three dimensions.
///
/// When the normals are not parallel the meeting set is a line whose
/// direction is the cross product of the normals. The point on that line is
/// found by sliding along the component of the second normal that is
/// orthogonal to the first. That component is the right transversal to use
/// precisely because it is NOT orthogonal to the second normal: its dot
/// product with the second normal is `1 - (n1 . n2)^2`, which a non-parallel
/// pair of unit normals never makes zero. A transversal chosen to be
/// orthogonal to both normals instead would divide by zero on every pair
/// whose normals are coordinate-aligned, which is most of them.
fn flatPair3(comptime T: type, first_point: [3]T, first_normal: [3]T, second_point: [3]T, second_normal: [3]T, policy: PredicatePolicy) Intersection3(T) {
    if (normalsParallel3(T, first_normal, second_normal, policy)) {
        if (flatSideSign(3, T, second_point, first_point, first_normal, policy) != .zero) return .none;
        return .{ .plane = .{ .origin = fromArrayOf(Point3(T), first_point), .normal = fromArrayOf(Vector3(T), first_normal) } };
    }
    const direction = arrayNormalized(3, T, cross3Raw(T, first_normal, second_normal));
    const transversal = arrayNormalized(3, T, arraySub(3, T, second_normal, arrayScale(3, T, first_normal, arrayDot(3, T, second_normal, first_normal))));
    const gap = arrayDot(3, T, arraySub(3, T, second_point, first_point), second_normal);
    const step = divRaw(T, gap, arrayDot(3, T, transversal, second_normal));
    const meeting = arrayAdd(3, T, first_point, arrayScale(3, T, transversal, step));
    return .{ .line = .{
        .origin = fromArrayOf(Point3(T), meeting),
        .direction = fromArrayOf(Vector3(T), direction),
    } };
}

fn flatPair4(comptime T: type, first_point: [4]T, first_normal: [4]T, second_point: [4]T, second_normal: [4]T, policy: PredicatePolicy) Intersection4(T) {
    if (normalsParallel4(T, first_normal, second_normal, policy)) {
        if (flatSideSign(4, T, second_point, first_point, first_normal, policy) != .zero) return .none;
        return .{ .hyperplane = .{ .origin = fromArrayOf(Point4(T), first_point), .normal = fromArrayOf(Vector4(T), first_normal) } };
    }
    // `in_plane_normal` is the part of the second normal orthogonal to the
    // first. It slides between the two hyperplanes, and it is orthogonal to
    // the first by construction, so the step lands on both.
    const in_plane_normal = arrayNormalized(4, T, arraySub(4, T, second_normal, arrayScale(4, T, first_normal, arrayDot(4, T, second_normal, first_normal))));
    const gap = arrayDot(4, T, arraySub(4, T, second_point, first_point), second_normal);
    const step = divRaw(T, gap, arrayDot(4, T, in_plane_normal, second_normal));
    const meeting = arrayAdd(4, T, first_point, arrayScale(4, T, in_plane_normal, step));
    // The meeting set has direction span the two-dimensional orthogonal
    // complement of both normals, and in four dimensions that complement is
    // exactly two-dimensional, so two directions describe it.
    const first_direction = seedOutsideBasis(4, T, &[_][4]T{ first_normal, in_plane_normal });
    const second_direction = seedOutsideBasis(4, T, &[_][4]T{ first_normal, in_plane_normal, first_direction });
    return .{ .plane = .{
        .origin = fromArrayOf(Point4(T), meeting),
        .u = fromArrayOf(Vector4(T), first_direction),
        .v = fromArrayOf(Vector4(T), second_direction),
    } };
}
