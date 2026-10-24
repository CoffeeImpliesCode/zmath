//! Concrete geometric algebras.
//!
//! Twelve types and no generic algebra. Every public name below takes a
//! concrete type and a concrete element type; there is no `Algebra(p, q, r, T)`,
//! no generic vector, no generic quaternion, and no `anytype` product. The one
//! generic helper in this file, `Multivector`, is private: it is instantiated
//! once per algebra and each instantiation is a distinct concrete type with
//! its own methods.
//!
//! Basis order, metric, and coefficient count
//!
//! | Module | Ordered basis | Metric | Coefficients | Storage |
//! | --- | --- | --- | --- | --- |
//! | `Cl3` | e1, e2, e3 | +1, +1, +1 | 8 | 32 B (f32), 64 B (f64) |
//! | `Cl4` | e1, e2, e3, e4 | +1, +1, +1, +1 | 16 | 64 B (f32), 128 B (f64) |
//! | `Pga3` | e0, e1, e2, e3 | 0, +1, +1, +1 | 16 | 64 B (f32), 128 B (f64) |
//! | `Pga4` | e0, e1, e2, e3, e4 | 0, +1, +1, +1, +1 | 32 | 128 B (f32), 256 B (f64) |
//! | `Cga3` | e1, e2, e3, e_plus, e_minus | +1, +1, +1, +1, -1 | 32 | 128 B (f32), 256 B (f64) |
//! | `Cga4` | e1, e2, e3, e4, e_plus, e_minus | +1, +1, +1, +1, +1, -1 | 64 | 256 B (f32), 512 B (f64) |
//!
//! The basis order assigns one bit to each basis vector: basis vector `i` is
//! `1 << i`, the scalar is mask `0`, and the pseudoscalar is mask
//! `coefficientCount - 1`. A full multivector stores its coefficients in
//! ascending mask order, so the coefficient at index `mask` multiplies the
//! ascending canonical blade that `mask` represents. Grade is the population
//! count of the mask.
//!
//! The largest value type is `Cga4d`: 512 bytes on the stack. A multivector is
//! a value type and this file never nests one in a loop body that recurses.
//! Every operation uses a fixed `[coefficientCount]T` result array and no
//! dynamic allocation. The product rules are built once at comptime from the
//! metric; nothing is built at runtime.
//!
//! The duals
//!
//! `Cl` and `Cga` use the right dual `A * inverse(I)`, where `I` is the
//! pseudoscalar. Its square is `I * I = (-1)^(n * (n - 1) / 2) * product(metric)`,
//! which is `+1` for `Cl4` and `Cga4` and `-1` for `Cl3` and `Cga3`, so the
//! right dual satisfies
//!
//! ```text
//! dual(dual(A)) == I * I * A == dualSquare * A
//! ```
//!
//! with `dualSquare = (-1)^n` for all four of those algebras. `Pga` uses the
//! signed complementary Poincare dual: a grade `k` blade of mask `m` maps to
//! the blade of mask `top & ~m` with the sign `(-1)^(k * (n - 1 - k))`. Applying
//! it twice gives the sign
//!
//! ```text
//! (-1)^(k * (n - 1 - k) + (n - k) * (k - 1)) = (-1)^(2k * (n - k) - n) = (-1)^n
//! ```
//!
//! because `2k * (n - k)` is even. The same `dualSquare = (-1)^n` therefore
//! holds for the signed Poincare dual: it is an involution exactly when the
//! basis has an even number of vectors, and it negates the multivector exactly
//! when the count is odd. Every one of the twelve types publishes `dualSquare`
//! and satisfies `dual(dual(A)) == dualSquare * A`.
//!
//! Special values
//!
//! These operations are exact linear algebra on coefficients. `add`, `sub`,
//! `neg`, `scale`, `wedge`, `inner`, `meet`, `dot`, `grade`, `bladeValue`, and
//! `reverse` propagate a NaN or an infinity into every coefficient they touch.
//! `geometricProduct` skips a row whose left coefficient is exactly zero, so
//! `zero() * A` is exactly `zero()` even when `A` holds a non-finite entry: the
//! algebraic result wins there, not the IEEE one. No operation reads,
//! validates, sets, or restores `fenv`, so the caller must keep
//! round-to-nearest-ties-to-even active. `normalized`, `inverseChecked`, and
//! the `unit` constructors reject invalid scale values with `null`. They do not
//! guarantee finite output coefficients; final scaling is not checked.
//!
//! `inverse`, `inverseChecked`, and `normalized` divide a nonzero finite scale.
//! A binary32 quotient of a nonzero divisor is normal whenever it is nonzero,
//! so binary32 keeps the hardware divide. A binary64 reciprocal of a scale near
//! the top of the range is subnormal, and a soft-float binary64 divide flushes
//! that to a signed zero, so the binary64 reciprocal goes through the shared
//! exact quotient packer of `scalar_div.zig` and keeps its subnormal value.
//!
//! Normalization table
//!
//! Every row is implemented. Scale is the factor the operation applies, sign is
//! the grade sign it uses, and failure is the only way it can fail: no other
//! row reports failure, and every row propagates a non-finite coefficient from
//! its inputs.
//!
//! | Operation | Output | Scale and sign | Precondition | Failure |
//! | --- | --- | --- | --- | --- |
//! | `zero` | all coefficients `+0` | exact `+0` | none | none |
//! | `one` | mask `0` coefficient `+1` | exact `+1` | none | none |
//! | `pseudoscalar` | top mask coefficient `+1` | exact `+1` | none | none |
//! | `blade` | one basis vector | `+1`, no metric sign | index below `basisCount` | none |
//! | `fromCoefficients`, `setCoefficients`, `coefficients` | as supplied | none | none | none |
//! | `add`, `sub`, `neg` | coefficientwise | exact, sign follows the operands | none | none |
//! | `scale` | coefficientwise | sign follows the scalar | none | none |
//! | `geometricProduct` | all grades | metric signs from the comptime table | none | none |
//! | `wedge` | grade `r + s` | metric signs, disjoint blades only | none | none |
//! | `inner` | grade `r - s` | metric signs, `b` inside `a` only | none | none |
//! | `meet` | grade `r + s - n` | via the duals of this type | none | none |
//! | `dot` | scalar | grade-0 of `A * B(reverse)` | none | none |
//! | `grade` | scalar | plain sum over one grade | none | none |
//! | `bladeValue` | scalar | the stored coefficient | none | none |
//! | `scalarPart` | scalar | the mask `0` coefficient | none | none |
//! | `vectorPart`, `euclideanPart` | `[dimension]T` | basis order, no sign | none | none |
//! | `reverse` | same grades | `(-1)^(k * (k - 1) / 2)` | none | none |
//! | `conjugate` | same grades | `(-1)^(k * (k + 1) / 2)`: negative on grades 1, 2, 5, 6 | none | none |
//! | `dual` | grade `n - k` | right dual, or signed Poincare | none | none |
//! | `dualComplement` | grade `n - k` | unsigned complement | none | none |
//! | `evenPart` | all even grades | exact | none | none |
//! | `norm` | scalar | grade-0 of `A * A(reverse)` | none | none |
//! | `normalized` | unit multivector | divide by the square root of `norm` | `norm` positive and finite | `null` |
//! | `inverse` | candidate inverse | even conjugate over scalar sandwich | finite nonzero scalar sandwich | unchecked |
//! | `inverseChecked` | as `inverse` | checks input and sandwich, not final coefficients | as `inverse` | `null` |
//! | `isFinite` | bool | exact | none | none |
//! | `unitRotor` | normalized even part of `1 + b` | normalize | caller establishes rotor membership | `null` |
//! | `unitVersor` | normalized even part | normalize | caller establishes versor membership | `null` |
//! | `unitFlat` | normalized single grade | normalize | caller establishes blade decomposability | `null` |
//! | `motorPart` | grade `4` | exact | none | none |
//! | `point` (Pga) | grade `2` | `e0 wedge x`, `e0` coefficient `+1` | none | none |
//! | `n_o` (Cga) | grade `1`, `(e_minus - e_plus) / 2` | exact, fixed | none | none |
//! | `n_inf` (Cga) | grade `1`, `e_minus + e_plus` | exact, fixed | none | none |
//! | `origin` (Cga) | grade `1`, and exactly `n_o` | exact, fixed | none | none |
//! | `embed` (Cga) | grade `1`: `n_o + x + (x dot x) * n_inf / 2` | exact, `+1` | none | none |
//! | `unitRound` (Cga) | unit grade `n - 1` | normalize that grade | that grade is present | `null` |
//!
//! The normalization table does not name a PGA line or plane operation, and
//! that is a property of the algebra rather than a gap. In `Pga` a point is
//! the grade-2 blade `e0 wedge x`, two point blades multiply to exactly zero,
//! every grade-2 blade has an identically zero scalar part against a point
//! blade, and the grade-2 part of `L * P` is one linear functional of `x` whose
//! kernel is the direction of `L`. No single `Pga` blade therefore carries a
//! point-incidence test, which is why `zmath.geometry` owns the line and plane
//! types.
//!
//! The twelve types share one code path, so the operations that only a family
//! has are declared on every type and reject the wrong family at compile time
//! rather than being absent. `point` is projective only, and `n_o`, `n_inf`,
//! `origin`, `embed`, and `unitRound` are conformal only. Every other name in
//! the table exists on all twelve types, and `motorPart` and `vectorPart` are
//! meaningful on all three families.
//!
//! ## A shortest use
//!
//! ```zig
//! const std = @import("std");
//! const zmath = @import("zmath");
//!
//! pub fn rotorAndRefusals() !void {
//!     const A = zmath.algebra.Cl3d;
//!
//!     // A pure bivector, stored as the coefficient of the blade mask.
//!     var coefficients: [A.coefficient_count]f64 = @splat(0);
//!     coefficients[0b011] = 1.0; // e1 wedge e2
//!     const rotor = A.unitRotor(A.fromCoefficients(coefficients)) orelse
//!         return error.NoRotor;
//!     // The inverse of a unit versor is its conjugate, so the product is one.
//!     const product = rotor.geometricProduct(rotor.inverse());
//!     std.debug.assert(@abs(product.scalarPart() - 1.0) < 1e-12);
//!
//!     // A zero multivector has neither a unit form nor an inverse, and both
//!     // report the reason as `null` rather than writing a non-finite answer.
//!     std.debug.assert(A.zero().normalized() == null);
//!     std.debug.assert(A.zero().inverseChecked() == null);
//! }
//! ```

const std = @import("std");
const scalar = @import("scalar.zig");
const scalar_div = @import("scalar_div.zig");

/// The three families this file instantiates.
pub const Family = enum { clifford, projective, conformal };

/// The widest basis any of the twelve types uses, `Cga4`.
const max_basis = 6;
const MaskShift = std.math.Log2Int(usize);

/// The blade mask of the basis vector at `index`.
fn bit(index: usize) usize {
    return @as(usize, 1) << @as(MaskShift, @intCast(index));
}

/// The grade of the blade at `mask`, its population count.
fn gradeOf(mask: usize) usize {
    return @popCount(mask);
}

/// `+1` for an even exponent and `-1` for an odd one, for any exponent.
fn paritySign(exponent: i32) i8 {
    return if (@mod(exponent, 2) == 0) @as(i8, 1) else @as(i8, -1);
}

/// The metric sign of one basis vector.
fn metricSign(comptime T: type, value: T) i8 {
    if (value > 0) return 1;
    if (value < 0) return -1;
    return 0;
}

/// The accepted square root kernel for the element type.
fn Sqrt(comptime T: type, x: T) T {
    return switch (T) {
        f32 => scalar.sqrtf(x),
        f64 => scalar.sqrtd(x),
        else => @compileError("zmath.algebra supports f32 and f64 only"),
    };
}

/// The exact finite test for the element type.
fn IsFinite(comptime T: type, x: T) bool {
    return switch (T) {
        f32 => std.math.isFinite(x),
        f64 => std.math.isFinite(x),
        else => @compileError("zmath.algebra supports f32 and f64 only"),
    };
}

/// The reciprocal of a nonzero coefficient. The binary32 quotient of a
/// nonzero divisor is normal whenever it is nonzero, so it keeps the
/// hardware divide. The binary64 reciprocal of a coefficient near the top of
/// the range is subnormal, and a soft-float binary64 divide flushes that to a
/// signed zero, so it goes through the shared exact quotient packer.
fn Reciprocal(comptime T: type, x: T) T {
    return switch (T) {
        f32 => 1 / x,
        f64 => scalar_div.div64(1.0, x),
        else => @compileError("zmath.algebra supports f32 and f64 only"),
    };
}

/// One basis-vector product, folded to a mask and a `+1`, `-1`, or `0` scale.
const ProductEntry = struct {
    mask: u8,
    coeff: i8,
};

/// One dual image, folded to a mask and a scale.
fn DualEntry(comptime T: type) type {
    return struct { mask: u8, scale: T };
}

/// Builds one concrete geometric algebra.
///
/// This helper is private. Each call returns a distinct concrete type, so the
/// twelve published types have no shared generic surface.
fn Multivector(
    comptime n: usize,
    comptime T: type,
    comptime kind: Family,
    comptime basis_metric: [n]T,
    comptime euclidean_dimension: usize,
) type {
    return struct {
        const Self = @This();
        /// The coefficients, indexed by ascending blade mask.
        c: [coefficient_count]T,

        /// The element type of every coefficient.
        pub const element = T;
        /// The number of basis vectors, and so the log2 of the coefficient count.
        pub const basis_count = n;
        /// The number of coefficients, one per blade mask.
        pub const coefficient_count = 1 << n;
        /// The mask of the scalar.
        pub const scalar_mask = 0;
        /// The mask of the top blade, the pseudoscalar.
        pub const pseudoscalar_mask = coefficient_count - 1;

        // `coefficient_count` is a left shift of one, so it is never zero: the
        // narrowest multivector still has the scalar blade, and every wider one
        // has the pseudoscalar above it. `pseudoscalar_mask` is that count less
        // one, so it is a real index into every `coefficient_count` wide table.
        comptime {
            if (coefficient_count < 1) @compileError("a multivector has at least the scalar blade");
        }

        /// The number of distinct grades, from `0` to `n`.
        pub const grade_count = n + 1;
        /// The Euclidean dimension, which is the basis count less the null vector
        /// and, for a conformal algebra, the two extra directions.
        pub const dimension = euclidean_dimension;
        /// The family this algebra belongs to.
        pub const family = kind;
        /// The metric of the ordered basis.
        pub const metric = basis_metric;
        /// The signed square of the pseudoscalar.
        pub const pseudoscalar_square: T = pseudoscalarSquare();
        /// The sign of `dual(dual(A))` relative to `A`, which is `(-1)^n` for
        /// every one of the twelve types.
        pub const dual_square: T = if (@mod(n, 2) == 0) 1 else -1;
        /// The size of the coefficient array, in bytes.
        pub const storage_bytes = @sizeOf([coefficient_count]T);

        /// The bit index of the `i`th Euclidean direction. `Pga` leads with the
        /// null vector, so its directions start at bit 1.
        fn spatialBit(index: usize) usize {
            return if (kind == .projective) index + 1 else index;
        }

        /// The mask of the grade-1 blade for the `i`th Euclidean direction.
        fn spatialMask(index: usize) usize {
            return bit(spatialBit(index));
        }

        /// `I * I` for the ordered basis, the signed square of the top blade.
        fn pseudoscalarSquare() T {
            var magnitude: i32 = 1;
            if (@mod(n * (n - 1) / 2, 2) != 0) magnitude = -1;
            for (0..n) |index| {
                magnitude *= metricSign(T, basis_metric[index]);
            }
            return @as(T, @floatFromInt(magnitude));
        }

        /// The product of two blades, folded with the metric signs.
        ///
        /// The two basis index lists concatenate, sort into ascending order with
        /// the inversion count folded into the sign, and then collapse every
        /// adjacent equal pair into its metric sign. Sorting first is valid
        /// because `e_i * e_j = -e_j * e_i` for `i != j`, and an equal pair on
        /// both sides is the same basis vector, so it collapses to its metric.
        fn bladeProduct(a_mask: usize, b_mask: usize) ProductEntry {
            var list: [2 * max_basis]usize = undefined;
            var length: usize = 0;
            for (0..n) |index| {
                if (a_mask & bit(index) != 0) {
                    list[length] = index;
                    length += 1;
                }
            }
            for (0..n) |index| {
                if (b_mask & bit(index) != 0) {
                    list[length] = index;
                    length += 1;
                }
            }
            var sign: i8 = 1;
            var index: usize = 1;
            while (index < length) : (index += 1) {
                var cursor = index;
                while (cursor > 0 and list[cursor - 1] > list[cursor]) : (cursor -= 1) {
                    const held = list[cursor - 1];
                    list[cursor - 1] = list[cursor];
                    list[cursor] = held;
                    sign = -sign;
                }
            }
            var mask: usize = 0;
            var magnitude: i8 = sign;
            var collapse: usize = 0;
            while (collapse < length) {
                if (collapse + 1 < length and list[collapse] == list[collapse + 1]) {
                    magnitude *%= metricSign(T, basis_metric[list[collapse]]);
                    collapse += 2;
                } else {
                    mask |= bit(list[collapse]);
                    collapse += 1;
                }
            }
            return .{ .mask = @intCast(mask), .coeff = magnitude };
        }

        const product_table: [coefficient_count][coefficient_count]ProductEntry =
            buildProductTable();

        fn buildProductTable() [coefficient_count][coefficient_count]ProductEntry {
            // The widest table is 64 by 64 pairs, each a short sort and fold.
            @setEvalBranchQuota(200_000);
            var table: [coefficient_count][coefficient_count]ProductEntry = undefined;
            for (0..coefficient_count) |a_mask| {
                for (0..coefficient_count) |b_mask| {
                    table[a_mask][b_mask] = bladeProduct(a_mask, b_mask);
                }
            }
            return table;
        }

        /// The dual image of one blade, folded to a mask and a scale.
        ///
        /// `Cl` and `Cga` multiply by the multivector inverse of the
        /// pseudoscalar, whose coefficient is the sign the product table already
        /// holds for the top blade divided by the signed square of that
        /// pseudoscalar. `Pga` takes the sign `(-1)^(k * (n - 1 - k))` of the
        /// signed complementary Poincare dual.
        fn dualEntry(mask: usize) DualEntry(T) {
            const complement = (coefficient_count - 1) & ~mask;
            const coefficient: T = switch (kind) {
                .projective => blk: {
                    const basis: i32 = @intCast(n);
                    const k: i32 = @intCast(gradeOf(mask));
                    break :blk @floatFromInt(paritySign(k * (basis - 1 - k)));
                },
                else => blk: {
                    const sign: i32 = product_table[mask][pseudoscalar_mask].coeff;
                    break :blk @as(T, @floatFromInt(sign)) / pseudoscalarSquare();
                },
            };
            return .{ .mask = @intCast(complement), .scale = coefficient };
        }

        const dual_table: [coefficient_count]DualEntry(T) = buildDualTable();

        fn buildDualTable() [coefficient_count]DualEntry(T) {
            @setEvalBranchQuota(200_000);
            var table: [coefficient_count]DualEntry(T) = undefined;
            for (0..coefficient_count) |mask| {
                table[mask] = dualEntry(mask);
            }
            return table;
        }

        /// The grade-0 part of the geometric product. This is the sandwich
        /// scalar of the even grades, and it is not `dot`, whose second argument
        /// is reversed.
        ///
        /// A blade pair lands on the scalar only when it multiplies with itself,
        /// so this is the diagonal of the product table summed in ascending
        /// mask order, with the same `+0` seed, the same zero left skip, the
        /// same zero scale skip, and the same operand order as the full
        /// product. Every other pair lands above grade `0` and is dropped, so
        /// this computes the same value as `a.geometricProduct(b).c[scalar_mask]`
        /// without the `coefficient_count` squared terms.
        fn scalarProduct(a: Self, b: Self) T {
            const table = &product_table;
            var sum: T = 0;
            for (0..coefficient_count) |mask| {
                const value = a.c[mask];
                if (value == 0) continue;
                const entry = table[mask][mask];
                if (entry.coeff == 0) continue;
                sum += @as(T, @floatFromInt(entry.coeff)) * value * b.c[mask];
            }
            return sum;
        }

        /// The versor conjugate, which is the Clifford conjugate restricted to
        /// the even grades. It keeps grade 0 and negates grades 2 and 4, so the
        /// sandwich of an even multivector with it has a scalar grade-0 part.
        fn sandwichConjugate(self: Self) Self {
            return self.conjugate().evenPart();
        }

        /// The blade of one grade, or `zero()` for a grade above the top one.
        fn gradedPart(self: Self, k: usize) Self {
            var result: [coefficient_count]T = @splat(0);
            if (k > n) return .{ .c = result };
            for (0..coefficient_count) |mask| {
                if (gradeOf(mask) == k) result[mask] = self.c[mask];
            }
            return .{ .c = result };
        }

        // ------------------------------------------------------------------
        // Construction and access
        // ------------------------------------------------------------------

        /// The multivector with every coefficient `+0`.
        pub fn zero() Self {
            return .{ .c = @splat(0) };
        }

        /// The multivector with the scalar coefficient `+1` and no other.
        pub fn one() Self {
            var result: [coefficient_count]T = @splat(0);
            result[scalar_mask] = 1;
            return .{ .c = result };
        }

        /// The multivector with the top blade coefficient `+1` and no other.
        pub fn pseudoscalar() Self {
            var result: [coefficient_count]T = @splat(0);
            result[pseudoscalar_mask] = 1;
            return .{ .c = result };
        }

        /// The basis vector at `index`, in the published basis order, with
        /// coefficient `+1` and no metric sign. An index at or above the basis
        /// count gives `zero()`.
        pub fn blade(index: usize) Self {
            var result: [coefficient_count]T = @splat(0);
            if (index < basis_count) result[bit(index)] = 1;
            return .{ .c = result };
        }

        /// The multivector with the supplied coefficients, in ascending mask
        /// order, without any normalization.
        pub fn fromCoefficients(values: [coefficient_count]T) Self {
            return .{ .c = values };
        }

        /// The coefficients, in ascending mask order.
        pub fn coefficients(self: Self) [coefficient_count]T {
            return self.c;
        }

        /// The multivector with the supplied coefficients, replacing every one.
        /// Nothing is mutated: the receiver arrives by value and the answer is a
        /// new multivector, so a caller keeps whatever `self` held.
        pub fn setCoefficients(self: Self, values: [coefficient_count]T) Self {
            _ = self;
            return .{ .c = values };
        }

        /// The coefficient of the blade at `mask`, which is `+0` for a mask at
        /// or above the coefficient count.
        pub fn bladeValue(self: Self, mask: usize) T {
            if (mask >= coefficient_count) return 0;
            return self.c[mask];
        }

        /// The sum of the coefficients of grade `k`, which is `+0` for a grade
        /// above the top one.
        pub fn grade(self: Self, k: usize) T {
            var sum: T = 0;
            if (k > n) return sum;
            for (0..coefficient_count) |mask| {
                if (gradeOf(mask) == k) sum += self.c[mask];
            }
            return sum;
        }

        /// The scalar coefficient.
        pub fn scalarPart(self: Self) T {
            return self.c[scalar_mask];
        }

        /// True when every coefficient is finite.
        pub fn isFinite(self: Self) bool {
            for (self.c) |value| {
                if (!IsFinite(T, value)) return false;
            }
            return true;
        }

        // ------------------------------------------------------------------
        // Linear algebra
        // ------------------------------------------------------------------

        /// The coefficientwise sum.
        pub fn add(a: Self, b: Self) Self {
            var result: [coefficient_count]T = @splat(0);
            for (0..coefficient_count) |mask| {
                result[mask] = a.c[mask] + b.c[mask];
            }
            return .{ .c = result };
        }

        /// The coefficientwise difference.
        pub fn sub(a: Self, b: Self) Self {
            var result: [coefficient_count]T = @splat(0);
            for (0..coefficient_count) |mask| {
                result[mask] = a.c[mask] - b.c[mask];
            }
            return .{ .c = result };
        }

        /// The coefficientwise negation.
        pub fn neg(self: Self) Self {
            var result: [coefficient_count]T = @splat(0);
            for (0..coefficient_count) |mask| {
                result[mask] = -self.c[mask];
            }
            return .{ .c = result };
        }

        /// The coefficientwise product with `factor`.
        pub fn scale(self: Self, factor: T) Self {
            var result: [coefficient_count]T = @splat(0);
            for (0..coefficient_count) |mask| {
                result[mask] = self.c[mask] * factor;
            }
            return .{ .c = result };
        }

        /// Retains every even grade, including grade 6 where supported.
        pub fn evenPart(self: Self) Self {
            var result: [coefficient_count]T = @splat(0);
            for (0..coefficient_count) |mask| {
                if (@mod(gradeOf(mask), 2) == 0) result[mask] = self.c[mask];
            }
            return .{ .c = result };
        }

        // ------------------------------------------------------------------
        // Products
        // ------------------------------------------------------------------

        /// The geometric product, the sum over every blade pair with the
        /// comptime sign of that pair.
        ///
        /// A row whose left coefficient is exactly zero is skipped, so a zero
        /// factor annihilates that term instead of producing a NaN.
        pub fn geometricProduct(a: Self, b: Self) Self {
            const table = &product_table;
            var result: [coefficient_count]T = @splat(0);
            for (0..coefficient_count) |left| {
                const value = a.c[left];
                if (value == 0) continue;
                for (0..coefficient_count) |right| {
                    const entry = table[left][right];
                    if (entry.coeff == 0) continue;
                    result[entry.mask] += @as(T, @floatFromInt(entry.coeff)) * value * b.c[right];
                }
            }
            return .{ .c = result };
        }

        /// The exterior product, the grade `r + s` part of the geometric
        /// product. A pair that shares a basis vector contributes nothing.
        pub fn wedge(a: Self, b: Self) Self {
            const table = &product_table;
            var result: [coefficient_count]T = @splat(0);
            for (0..coefficient_count) |left| {
                const value = a.c[left];
                if (value == 0) continue;
                for (0..coefficient_count) |right| {
                    if (left & right != 0) continue;
                    const entry = table[left][right];
                    if (entry.coeff == 0) continue;
                    result[entry.mask] += @as(T, @floatFromInt(entry.coeff)) * value * b.c[right];
                }
            }
            return .{ .c = result };
        }

        /// The interior product, the part of the geometric product whose right
        /// blade lies inside its left blade. It has grade `r - s` and it
        /// includes the scalar part.
        pub fn inner(a: Self, b: Self) Self {
            const table = &product_table;
            var result: [coefficient_count]T = @splat(0);
            for (0..coefficient_count) |left| {
                const value = a.c[left];
                if (value == 0) continue;
                for (0..coefficient_count) |right| {
                    if (right & ~left != 0) continue;
                    const entry = table[left][right];
                    if (entry.coeff == 0) continue;
                    result[entry.mask] += @as(T, @floatFromInt(entry.coeff)) * value * b.c[right];
                }
            }
            return .{ .c = result };
        }

        /// The meet, the dual of the wedge of the two duals. It has grade
        /// `r + s - n`, so it is `zero()` for a grade sum below the basis count,
        /// which is the usual degeneracy of the `Pga` flats.
        pub fn meet(a: Self, b: Self) Self {
            return dual(a).wedge(dual(b)).dual();
        }

        /// The scalar form, the grade-0 part of `A * B(reverse)`. It is
        /// symmetric, and on the vector grade it is the ordinary dot product.
        pub fn dot(a: Self, b: Self) T {
            return scalarProduct(a, b.reverse());
        }

        // ------------------------------------------------------------------
        // Involutions
        // ------------------------------------------------------------------

        /// The multivector with the grade sign `(-1)^(k * (k - 1) / 2)`.
        pub fn reverse(self: Self) Self {
            var result: [coefficient_count]T = @splat(0);
            for (0..coefficient_count) |mask| {
                const k: i32 = @intCast(gradeOf(mask));
                const sign: i32 = paritySign(@divTrunc(k * (k - 1), 2));
                result[mask] = @as(T, @floatFromInt(sign)) * self.c[mask];
            }
            return .{ .c = result };
        }

        /// The Clifford conjugate, with grade sign `(-1)^(k * (k + 1) / 2)`.
        /// In the supported dimensions it negates grades 1, 2, 5, and 6,
        /// and preserves grades 0, 3, and 4.
        ///
        /// This is distinct from `reverse`; it is not a general versor inverse.
        pub fn conjugate(self: Self) Self {
            var result: [coefficient_count]T = @splat(0);
            for (0..coefficient_count) |mask| {
                const k: i32 = @intCast(gradeOf(mask));
                const sign: T = @floatFromInt(paritySign(@divTrunc(k * (k + 1), 2)));
                result[mask] = sign * self.c[mask];
            }
            return .{ .c = result };
        }

        /// The dual of this algebra: the right dual `A * inverse(I)` for `Cl`
        /// and `Cga`, the signed complementary Poincare dual for `Pga`. It
        /// satisfies `dual(dual(A)) == dual_square * A`.
        pub fn dual(self: Self) Self {
            const table = &dual_table;
            var result: [coefficient_count]T = @splat(0);
            for (0..coefficient_count) |mask| {
                const entry = table[mask];
                result[entry.mask] += entry.scale * self.c[mask];
            }
            return .{ .c = result };
        }

        /// The unsigned complementary dual, which drops the dual sign.
        pub fn dualComplement(self: Self) Self {
            var result: [coefficient_count]T = @splat(0);
            for (0..coefficient_count) |mask| {
                result[(coefficient_count - 1) & ~mask] += self.c[mask];
            }
            return .{ .c = result };
        }

        /// Scales the even Clifford conjugate by the reciprocal sandwich scalar.
        ///
        /// Precondition: its product with `self` is a finite, nonzero scalar.
        /// Only then is the scaled value an inverse. This is not a general
        /// multivector inverse; an odd-only blade fails that precondition.
        /// Final scaling can overflow. `inverseChecked` checks the input and
        /// sandwich, but does not check the final scaled coefficients.
        pub fn inverse(self: Self) Self {
            const sandwich = self.sandwichConjugate();
            return sandwich.scale(Reciprocal(T, scalarProduct(sandwich, self)));
        }

        /// The same candidate as `inverse`, or `null` for non-finite input,
        /// a zero/non-finite sandwich scalar, or a non-scalar sandwich.
        /// An odd-only blade returns `null`; general invertibility is not
        /// tested. Final scaled coefficients can still be non-finite.
        pub fn inverseChecked(self: Self) ?Self {
            if (!self.isFinite()) return null;
            const conjugate_part = self.sandwichConjugate();
            const product = conjugate_part.geometricProduct(self);
            // The scan accumulates a flag instead of returning from inside the
            // loop. The accepted set is identical, but with the early return in
            // the loop the big-endian AArch64 backend builds this optional
            // aggregate return wrongly: `Pga4f` and `Cga3f`, the two 128-byte
            // binary32 payloads, then answer `null` for an inverse they accept,
            // while `normalized`, which returns the same optional shape and has
            // no loop, stays correct.
            var residue = false;
            for (0..coefficient_count) |mask| {
                if (mask != scalar_mask and product.c[mask] != 0) residue = true;
            }
            // The residue scan reads the whole product, and its grade-0 part is
            // the sandwich scalar, so the scale is the one `inverse` reads from
            // that same product instead of building a second one.
            const sandwich = product.c[scalar_mask];
            if (residue or sandwich == 0 or !IsFinite(T, sandwich)) return null;
            return conjugate_part.scale(Reciprocal(T, sandwich));
        }

        /// The squared norm, the grade-0 part of `A * A(reverse)`. It is
        /// nonnegative for `Cl` and `Cga`, and it may be negative or zero for
        /// `Pga`.
        pub fn norm(self: Self) T {
            return dot(self, self);
        }

        /// The multivector divided by the square root of its norm, or `null`
        /// when the norm is zero, negative, or non-finite, or when any
        /// coefficient is non-finite.
        pub fn normalized(self: Self) ?Self {
            if (!self.isFinite()) return null;
            const squared = self.norm();
            if (squared <= 0 or !IsFinite(T, squared)) return null;
            const magnitude = Sqrt(T, squared);
            if (magnitude == 0 or !IsFinite(T, magnitude)) return null;
            return self.scale(Reciprocal(T, magnitude));
        }

        // ------------------------------------------------------------------
        // Named multivector operations
        // ------------------------------------------------------------------

        /// Normalizes the even part of `1 + b`. For a rotor result, the caller
        /// must establish rotor membership; this does not check it. Returns
        /// `null` when the normalizer rejects its norm.
        pub fn unitRotor(b: Self) ?Self {
            return one().add(b).evenPart().normalized();
        }

        /// Normalizes the even part of `self`. The caller must establish that
        /// it is a versor; even grade alone does not establish membership.
        /// Returns `null` when the normalizer rejects its norm.
        pub fn unitVersor(self: Self) ?Self {
            return self.evenPart().normalized();
        }

        /// Normalizes a nonzero single-grade value. The caller must establish
        /// blade decomposability to use it as a flat; that is not checked.
        /// Returns `null` for mixed grades, non-finite coefficients, zero, or
        /// a norm the normalizer rejects.
        pub fn unitFlat(self: Self) ?Self {
            if (!self.isFinite()) return null;
            var seen: usize = 0;
            var seen_grade: usize = 0;
            for (0..coefficient_count) |mask| {
                if (self.c[mask] == 0) continue;
                const k = gradeOf(mask);
                if (seen != 0 and k != seen_grade) return null;
                seen += 1;
                seen_grade = k;
            }
            if (seen == 0) return null;
            return self.normalized();
        }

        /// The grade-4 part of `self`, the translation part of a motor. It is
        /// `zero()` for a basis that has no grade 4.
        pub fn motorPart(self: Self) Self {
            return self.gradedPart(4);
        }

        /// The Euclidean components, in the published basis order from `e1` to
        /// `e_dimension`. For `Pga` the leading null direction is not one of
        /// them, and for `Cga` the two conformal directions are not.
        pub fn vectorPart(self: Self) [dimension]T {
            var result: [dimension]T = undefined;
            for (0..dimension) |index| {
                result[index] = self.c[spatialMask(index)];
            }
            return result;
        }

        /// The `euclideanPart` spelling of `vectorPart`, which is the name the
        /// conformal row of the normalization table uses. The two are the same
        /// map in every family, because the spatial bits are the Euclidean
        /// directions in all three.
        pub const euclideanPart = vectorPart;

        /// The homogeneous point at the Euclidean position `x`, the grade-2
        /// blade `e0 wedge x` with the `e0` coefficient `+1`.
        ///
        /// Only `Pga` leads its basis with a null direction, so this operation
        /// exists only there and is a compile error elsewhere.
        pub fn point(x: [dimension]T) Self {
            if (kind != .projective) @compileError("zmath.algebra.point exists only for Pga");
            var result: [coefficient_count]T = @splat(0);
            for (0..dimension) |index| {
                result[bit(0) | spatialMask(index)] = x[index];
            }
            return .{ .c = result };
        }

        /// The conformal origin, the grade-1 blade `(e_minus - e_plus) / 2`: the
        /// coefficient `+0.5` on `e_minus` and `-0.5` on `e_plus`, which are the
        /// last two basis vectors, and nothing else. It is a sum of two grade-1
        /// basis coefficients and NOT a wedge, so it is not a grade-2 blade. Its
        /// inner product with `n_inf` is `-1`.
        ///
        /// Only `Cga` ends its basis with the two conformal directions, so this
        /// operation exists only there.
        pub fn n_o() Self {
            if (kind != .conformal) @compileError("zmath.algebra.n_o exists only for Cga");
            var result: [coefficient_count]T = @splat(0);
            result[bit(n - 1)] = @as(T, 0.5);
            result[bit(n - 2)] = -@as(T, 0.5);
            return .{ .c = result };
        }

        /// The conformal direction, the grade-1 blade `e_minus + e_plus`: the
        /// coefficient `+1` on each of the last two basis vectors and nothing
        /// else. It squares to zero, because the two conformal directions carry
        /// opposite metric signs, so `norm` of it is `0`.
        ///
        /// Only `Cga` ends its basis with the two conformal directions, so this
        /// operation exists only there.
        pub fn n_inf() Self {
            if (kind != .conformal) @compileError("zmath.algebra.n_inf exists only for Cga");
            var result: [coefficient_count]T = @splat(0);
            result[bit(n - 1)] = 1;
            result[bit(n - 2)] = 1;
            return .{ .c = result };
        }

        /// The point at the Euclidean origin, which is `n_o` itself: the grade-1
        /// blade with `+0.5` on `e_minus` and `-0.5` on `e_plus`, and no
        /// scalar part, no Euclidean part, and no `n_inf` part. It is the
        /// embedding of the zero vector, so `origin()` and `embed(@splat(0))`
        /// are the same multivector.
        ///
        /// Only `Cga` has a conformal origin, so this operation exists only
        /// there.
        pub fn origin() Self {
            return n_o();
        }

        /// The multivector with only the Euclidean components, used by `embed`.
        fn vectorOf(x: [dimension]T) Self {
            var result: [coefficient_count]T = @splat(0);
            for (0..dimension) |index| {
                result[spatialMask(index)] = x[index];
            }
            return .{ .c = result };
        }

        /// The embedding of the Euclidean vector `x`, which is
        /// `n_o + x + (x dot x) * n_inf / 2`.
        ///
        /// Only `Cga` has the conformal pair this embedding needs, so this
        /// operation exists only there.
        pub fn embed(x: [dimension]T) Self {
            if (kind != .conformal) @compileError("zmath.algebra.embed exists only for Cga");
            var squared: T = 0;
            for (x) |value| squared += value * value;
            return n_o().add(vectorOf(x)).add(n_inf().scale(squared / 2));
        }

        /// The unit round, the grade `n - 1` blade normalized.
        ///
        /// A round in the conformal model is the grade `n - 1` blade, so this
        /// operation exists only for `Cga`. Precondition: that grade carries
        /// the whole round. Failure: `null` when it is zero or non-finite.
        pub fn unitRound(self: Self) ?Self {
            if (kind != .conformal) @compileError("zmath.algebra.unitRound exists only for Cga");
            return self.gradedPart(n - 1).normalized();
        }
    };
}

/// The Euclidean geometric algebra on three vectors with `f32` coefficients.
pub const Cl3f = Multivector(3, f32, .clifford, .{ 1, 1, 1 }, 3);
/// The Euclidean geometric algebra on three vectors with `f64` coefficients.
pub const Cl3d = Multivector(3, f64, .clifford, .{ 1, 1, 1 }, 3);
/// The Euclidean geometric algebra on four vectors with `f32` coefficients.
pub const Cl4f = Multivector(4, f32, .clifford, .{ 1, 1, 1, 1 }, 4);
/// The Euclidean geometric algebra on four vectors with `f64` coefficients.
pub const Cl4d = Multivector(4, f64, .clifford, .{ 1, 1, 1, 1 }, 4);
/// The projective geometric algebra on three dimensions with `f32`
/// coefficients.
pub const Pga3f = Multivector(4, f32, .projective, .{ 0, 1, 1, 1 }, 3);
/// The projective geometric algebra on three dimensions with `f64`
/// coefficients.
pub const Pga3d = Multivector(4, f64, .projective, .{ 0, 1, 1, 1 }, 3);
/// The projective geometric algebra on four dimensions with `f32`
/// coefficients.
pub const Pga4f = Multivector(5, f32, .projective, .{ 0, 1, 1, 1, 1 }, 4);
/// The projective geometric algebra on four dimensions with `f64`
/// coefficients.
pub const Pga4d = Multivector(5, f64, .projective, .{ 0, 1, 1, 1, 1 }, 4);
/// The conformal geometric algebra on three dimensions with `f32`
/// coefficients.
pub const Cga3f = Multivector(5, f32, .conformal, .{ 1, 1, 1, 1, -1 }, 3);
/// The conformal geometric algebra on three dimensions with `f64`
/// coefficients.
pub const Cga3d = Multivector(5, f64, .conformal, .{ 1, 1, 1, 1, -1 }, 3);
/// The conformal geometric algebra on four dimensions with `f32`
/// coefficients.
pub const Cga4f = Multivector(6, f32, .conformal, .{ 1, 1, 1, 1, 1, -1 }, 4);
/// The conformal geometric algebra on four dimensions with `f64`
/// coefficients.
pub const Cga4d = Multivector(6, f64, .conformal, .{ 1, 1, 1, 1, 1, -1 }, 4);
