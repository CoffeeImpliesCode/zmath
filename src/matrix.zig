//! Fixed-size matrices: the one structural generic in the zmath API
//! (`API.md#scope`) and the module `API.md#matrices` specifies.
//!
//! ## Storage and access
//!
//! `Mat(rows, columns, Scalar)` stores `data: [columns][rows]Scalar`, so
//! `data[column][row]`, while `at(row, column)` takes the row first. The
//! transposition between the two is deliberate on both sides: a matrix is
//! WRITTEN row-major, and a column vector is APPLIED by walking down columns,
//! which is what this layout makes contiguous. It is also the one sharp edge in
//! this module, so `init` takes a row-major literal and `data` stays public
//! for a caller already holding column-major data.
//!
//! ## Composition order
//!
//! Vectors are columns and composition is right-to-left: `(A * B) * v` applies
//! `B` first and `A` second, so `A * B` reads "A after B". `mulVec` applies a
//! matrix to a column vector in the same sense.
//!
//! ## Accumulation order
//!
//! Every contraction sums its products in increasing contraction index,
//! starting from the index-zero product rather than from an explicit zero, and
//! every arithmetic body is `@setFloatMode(.strict)`, so the compiler may
//! neither reassociate a sum nor contract a multiply and an add into an FMA.
//! The order is fixed in the source, which makes a result reproducible instead
//! of merely correct. Starting from the first product is why `M * I` does not
//! preserve a negative zero entry: the leading `0 * M[r][k]` terms sum to
//! `+0`, and `+0 + (-0)` is `+0`. No other entry is affected.
//!
//! ## Determinant
//!
//! `determinant` expands along the first row by cofactors, recursively down to
//! a closed form at size two and three. It does NOT use Gaussian elimination.
//! Elimination is faster and better conditioned, but it makes the result depend
//! on a pivot order that a single rounding difference can change, and it pays
//! for divisions where cofactor expansion has none. Cofactor expansion has a
//! fixed operation order that depends on nothing but the entries, which is what
//! makes a determinant checkable against a hand derivation. Its cost is real
//! and is not hidden here: the size-4 form forms 40 multiplications and the
//! size-5 form 205, against the 2 and 9 of the closed forms, and the relative
//! error still grows with the condition number just as it does everywhere
//! else.
//!
//! ## Inverse
//!
//! The inverses use Gauss-Jordan elimination with PARTIAL PIVOTING, which is a
//! different algorithm on purpose: it is more stable than the adjugate at size
//! four, and its singularity test is a genuine pivot test rather than a test on
//! a computed determinant. The limit is the ordinary one, and no doc comment
//! here pretends otherwise. The inverse of an ill-conditioned matrix is not
//! accurate, and the bits `inverseMat4d` returns for a near-singular input are
//! set by the condition number rather than by this implementation. A caller who
//! needs a bound has to measure the condition number, which this module does
//! not provide.
//!
//! ## Special values
//!
//! Arithmetic on a real entry follows the encodings and the NaN rules of
//! `spec/core-math-2026-08-22.md#special-values`: a NaN operand propagates
//! with its sign and diagnostic payload intact and its quiet bit set, and an
//! invalid operation on non-NaN operands gives the POSITIVE canonical quiet
//! NaN of the format, which a host unit does not produce on its own. A complex
//! entry is `zmath.complex` arithmetic, which already meets the same rules, so
//! this module composes that kernel rather than restating it. zmath never
//! reads, validates, sets, or restores `fenv`; the caller must keep
//! round-to-nearest, ties-to-even active.
//!
//! ## Sizes
//!
//! `Mat4f` and `Mat4d` represent homogeneous 3D transforms and linear 4D
//! transforms. `Mat5f` and `Mat5d` represent homogeneous 4D affine or
//! projective transforms. That is why those two sizes exist. The arithmetic
//! treats every row alike and knows nothing about the interpretation.
//!
//! No operation here allocates. `Mat` is a comptime-bounded array.
//!
//! ## A shortest use
//!
//! ```zig
//! const std = @import("std");
//! const zmath = @import("zmath");
//!
//! pub fn applyAndReject() void {
//!     const rect = zmath.matrix.Mat(2, 4, f32).init(.{
//!         .{ 1, 3, 5, 7 },
//!         .{ 2, 4, 6, 8 },
//!     });
//!     // `init` takes a ROW-MAJOR literal and `data` is COLUMN-MAJOR, so
//!     // `data[1]` is the second COLUMN and holds row 0 above row 1.
//!     std.debug.assert(rect.data[1][0] == 3);
//!
//!     // A two-row matrix takes a four-lane COLUMN and returns a two-lane one.
//!     const applied: zmath.vector.Vec2f =
//!         rect.mulVec(zmath.vector.Vec4f.init(1, 1, 1, 1));
//!     std.debug.assert(applied.eq(zmath.vector.Vec2f.init(16, 20)));
//!
//!     // A zero matrix has no pivot above the threshold, so the inverse is
//!     // reported as `null` rather than as a zero or a non-finite matrix.
//!     const inverse = zmath.matrix.inverseMat2f(zmath.matrix.Mat2f.fill(0));
//!     std.debug.assert(inverse == null);
//! }
//! ```

const std = @import("std");
const complex = @import("complex.zig");
const vector = @import("vector.zig");

const Complexf = complex.Complexf;
const Complexd = complex.Complexd;
const Vec2f = vector.Vec2f;
const Vec2d = vector.Vec2d;
const Vec3f = vector.Vec3f;
const Vec3d = vector.Vec3d;
const Vec4f = vector.Vec4f;
const Vec4d = vector.Vec4d;

// ---------------------------------------------------------------------------
// Entry types
//
// `Component` is the private entry-type guard. Unsupported types are checked
// at the actual `Mat` constructor boundary, not through a public test probe.
// ---------------------------------------------------------------------------

/// The real type each component of a `Scalar` entry has. It is `Scalar` itself
/// for a real entry. A determinant keeps the entry type, so a complex matrix
/// has a complex determinant and this is not that type.
fn Component(comptime Scalar: type) type {
    return switch (Scalar) {
        f32, f64 => Scalar,
        Complexf, Complexd => RealOf(Scalar),
        else => @compileError(
            "zmath.matrix: unsupported entry type '" ++ @typeName(Scalar) ++
                "'. Expected f32, f64, Complexf, or Complexd.",
        ),
    };
}

/// The real type inside `Complexf` or `Complexd`.
fn RealOf(comptime Scalar: type) type {
    if (Scalar == Complexf) return f32;
    if (Scalar == Complexd) return f64;
    @compileError("zmath.matrix: '" ++ @typeName(Scalar) ++ "' is not a complex entry type");
}

/// The additive identity of `Scalar`, with positive zero components.
fn zeroOf(comptime Scalar: type) Scalar {
    return switch (Scalar) {
        f32, f64 => 0,
        Complexf, Complexd => .{ .re = 0, .im = 0 },
        else => @compileError("zmath.matrix: '" ++ @typeName(Scalar) ++ "' is not an entry type"),
    };
}

/// The multiplicative identity of `Scalar`, with a positive zero imaginary
/// part.
fn oneOf(comptime Scalar: type) Scalar {
    return switch (Scalar) {
        f32, f64 => 1,
        Complexf, Complexd => .{ .re = 1, .im = 0 },
        else => @compileError("zmath.matrix: '" ++ @typeName(Scalar) ++ "' is not an entry type"),
    };
}

// ---------------------------------------------------------------------------
// Entry arithmetic
//
// A real entry is IEEE arithmetic with the NaN rules of
// `spec/core-math-2026-08-22.md#special-values`; a complex entry is the matching
// `zmath.complex` kernel, which already implements those rules for both of its
// components. Neither path calls a host libm.
// ---------------------------------------------------------------------------

/// The unsigned integer of the same width as `T`, for the NaN masks.
fn Bits(comptime T: type) type {
    return switch (@bitSizeOf(T)) {
        32 => u32,
        64 => u64,
        else => @compileError("zmath.matrix: an entry component is f32 or f64"),
    };
}

/// The positive canonical quiet NaN of `T`, from the special-value table.
fn canonicalNaN(comptime T: type) T {
    return @bitCast(@as(Bits(T), if (T == f32) 0x7fc0_0000 else 0x7ff8_0000_0000_0000));
}

/// `value` with the quiet bit set, which promotes a signaling NaN and leaves a
/// quiet one alone. The sign and the diagnostic payload survive.
fn quietNaN(comptime T: type, value: T) T {
    const mask: Bits(T) = if (T == f32) 0x0040_0000 else 0x0008_0000_0000_0000;
    return @bitCast(@as(Bits(T), @bitCast(value)) | mask);
}

/// The result of a binary entry operation, given what the host unit produced.
///
/// The substitution is not decoration. An invalid SSE operation yields the
/// NEGATIVE indefinite, 0xffc00000, while the specification pins the positive
/// encoding, so without this the same program would emit a different NaN
/// payload on a different target.
fn binaryResult(comptime T: type, a: T, b: T, raw: T) T {
    if (!std.math.isNan(raw)) return raw;
    if (std.math.isNan(a)) return quietNaN(T, a);
    if (std.math.isNan(b)) return quietNaN(T, b);
    return canonicalNaN(T);
}

/// The result of a unary entry operation, given what the host unit produced.
fn unaryResult(comptime T: type, a: T, raw: T) T {
    if (!std.math.isNan(raw)) return raw;
    if (std.math.isNan(a)) return quietNaN(T, a);
    return canonicalNaN(T);
}

/// `a + b`.
fn entryAdd(comptime Scalar: type, a: Scalar, b: Scalar) Scalar {
    return switch (Scalar) {
        f32 => binaryResult(f32, a, b, a + b),
        f64 => binaryResult(f64, a, b, a + b),
        Complexf => complex.addf(a, b),
        Complexd => complex.addd(a, b),
        else => @compileError("zmath.matrix: '" ++ @typeName(Scalar) ++ "' is not an entry type"),
    };
}

/// `a - b`.
fn entrySub(comptime Scalar: type, a: Scalar, b: Scalar) Scalar {
    return switch (Scalar) {
        f32 => binaryResult(f32, a, b, a - b),
        f64 => binaryResult(f64, a, b, a - b),
        Complexf => complex.subf(a, b),
        Complexd => complex.subd(a, b),
        else => @compileError("zmath.matrix: '" ++ @typeName(Scalar) ++ "' is not an entry type"),
    };
}

/// `a * b`.
fn entryMul(comptime Scalar: type, a: Scalar, b: Scalar) Scalar {
    return switch (Scalar) {
        f32 => binaryResult(f32, a, b, a * b),
        f64 => binaryResult(f64, a, b, a * b),
        Complexf => complex.mulf(a, b),
        Complexd => complex.muld(a, b),
        else => @compileError("zmath.matrix: '" ++ @typeName(Scalar) ++ "' is not an entry type"),
    };
}

/// `-a`. A complex entry negates BOTH components, which is not the same as its
/// conjugate: `-conj(a)` would be right, `conj(a)` alone would be wrong, and
/// the difference is the sign of the real part. `zmath.complex` has no `neg`,
/// so the negation is `0 - a` and inherits the module's NaN rules.
fn entryNeg(comptime Scalar: type, a: Scalar) Scalar {
    return switch (Scalar) {
        f32 => unaryResult(f32, a, -a),
        f64 => unaryResult(f64, a, -a),
        Complexf => complex.subf(.{ .re = 0, .im = 0 }, a),
        Complexd => complex.subd(.{ .re = 0, .im = 0 }, a),
        else => @compileError("zmath.matrix: '" ++ @typeName(Scalar) ++ "' is not an entry type"),
    };
}

/// True when every component of `value` is finite.
fn isFiniteEntry(comptime Scalar: type, value: Scalar) bool {
    return switch (Scalar) {
        f32, f64 => std.math.isFinite(value),
        Complexf => std.math.isFinite(value.re) and std.math.isFinite(value.im),
        Complexd => std.math.isFinite(value.re) and std.math.isFinite(value.im),
        else => @compileError("zmath.matrix: '" ++ @typeName(Scalar) ++ "' is not an entry type"),
    };
}

// ---------------------------------------------------------------------------
// Vectors
//
// A matrix uses a named `Vec*` value when one matches the lane count and the
// entry type, and a plain array otherwise. Always taking a plain array would
// give a caller holding a `Vec3f` a second undocumented spelling of the same
// three numbers at the call site; inventing `Vec1f` and `Vec5f` would
// contradict `API.md#scope`.
// ---------------------------------------------------------------------------

/// The component names of the `Vec*` types, in declaration order.
const lane_names = [4][]const u8{ "x", "y", "z", "w" };

/// The input or output vector representation for the given matrix extent.
fn VecOf(comptime n_rows: usize, comptime Scalar: type) type {
    return switch (Scalar) {
        f32 => switch (n_rows) {
            2 => Vec2f,
            3 => Vec3f,
            4 => Vec4f,
            else => [n_rows]f32,
        },
        f64 => switch (n_rows) {
            2 => Vec2d,
            3 => Vec3d,
            4 => Vec4d,
            else => [n_rows]f64,
        },
        Complexf, Complexd => [n_rows]Scalar,
        else => @compileError("zmath.matrix: '" ++ @typeName(Scalar) ++ "' is not an entry type"),
    };
}

/// The component type of a `Vec*` or a plain vector array.
fn laneType(comptime V: type) type {
    return switch (@typeInfo(V)) {
        .array => |array| array.child,
        .@"struct" => |shape| shape.fields[0].type,
        else => @compileError("zmath.matrix: a vector is a Vec* value type or a plain array"),
    };
}

/// The component of `v` at `index`, which is comptime known.
fn laneAt(comptime V: type, v: V, comptime index: usize) laneType(V) {
    return switch (@typeInfo(V)) {
        .array => v[index],
        else => @field(v, lane_names[index]),
    };
}

/// Writes `value` into the component of `v` at `index`, which is comptime
/// known.
fn laneSet(comptime V: type, v: *V, comptime index: usize, value: laneType(V)) void {
    switch (@typeInfo(V)) {
        .array => v.*[index] = value,
        else => @field(v.*, lane_names[index]) = value,
    }
}

// ---------------------------------------------------------------------------
// Product typing
//
// The contract writes `mul` with undeclared `rows`, `columns`, and `Scalar`
// (`API.md#matrices`), which is shorthand for "the operands carry them". Zig has
// no syntax for inferring a trailing type parameter, so the operands are taken
// as `anytype` and the product type is computed and CHECKED here. That keeps
// the contracted call shape `mul(output_columns, left, right)` with nothing
// dynamic in it.
// ---------------------------------------------------------------------------

/// True when `T` declares a `pub` member called `name`.
fn hasDecl(comptime T: type, comptime name: []const u8) bool {
    for (@typeInfo(T).@"struct".decls) |decl| {
        if (std.mem.eql(u8, decl.name, name)) return true;
    }
    return false;
}

/// True when `T` is an instantiation of `Mat`.
fn isMat(comptime T: type) bool {
    if (@typeInfo(T) != .@"struct") return false;
    return hasDecl(T, "rows") and hasDecl(T, "columns") and hasDecl(T, "Element");
}

/// The product type of two `Mat` instantiations, or a compile error naming what
/// does not line up.
fn Product(comptime Left: type, comptime Right: type, comptime n_out_columns: usize) type {
    if (!isMat(Left)) {
        @compileError("zmath.matrix.mul: the left operand is not a Mat, it is '" ++ @typeName(Left) ++ "'");
    }
    if (!isMat(Right)) {
        @compileError("zmath.matrix.mul: the right operand is not a Mat, it is '" ++ @typeName(Right) ++ "'");
    }
    if (@as(usize, @field(Right, "rows")) != @as(usize, @field(Left, "columns"))) {
        @compileError(
            "zmath.matrix.mul: the operands do not chain. The right operand has " ++
                std.fmt.comptimePrint("{d}", .{@field(Right, "rows")}) ++
                " rows but the left operand has " ++
                std.fmt.comptimePrint("{d}", .{@field(Left, "columns")}) ++
                " columns.",
        );
    }
    if (@as(usize, @field(Right, "columns")) != n_out_columns) {
        @compileError(
            "zmath.matrix.mul: the output width is " ++
                std.fmt.comptimePrint("{d}", .{n_out_columns}) ++
                " but the right operand has " ++
                std.fmt.comptimePrint("{d}", .{@field(Right, "columns")}) ++
                " columns.",
        );
    }
    if (@field(Left, "Element") != @field(Right, "Element")) {
        @compileError(
            "zmath.matrix.mul: the operands have different entry types, '" ++
                @typeName(@field(Left, "Element")) ++ "' and '" ++
                @typeName(@field(Right, "Element")) ++ "'.",
        );
    }
    return Mat(@field(Left, "rows"), @field(Right, "columns"), @field(Left, "Element"));
}

// ---------------------------------------------------------------------------
// Mat
// ---------------------------------------------------------------------------

/// A `n_rows` by `n_columns` matrix of `Scalar` entries, stored column-major.
///
/// `Scalar` must be `f32`, `f64`, `Complexf`, or
/// `Complexd`, and any other entry type is a compile error that names it.
/// `n_rows` and `n_columns` are each at least one, so the narrowest matrix this
/// factory builds is 1 by 1 and a zero extent is a compile error too.
pub fn Mat(
    comptime n_rows: usize,
    comptime n_columns: usize,
    comptime Scalar: type,
) type {
    if (comptime n_rows == 0 or n_columns == 0) {
        @compileError("zmath.matrix.Mat: a matrix has at least one row and at least one column");
    }
    // Force the entry-type guard here, where the message can name the type,
    // rather than deep inside the first unrolled loop that uses a component.
    const component = Component(Scalar);
    const lane = VecOf(n_rows, Scalar);
    const input_lane = VecOf(n_columns, Scalar);

    return struct {
        const Self = @This();

        /// The matrix entry type.
        pub const Element = Scalar;
        /// The real type each component of an entry has. It is `Element` itself
        /// when the entry is real. A determinant keeps the entry type, so a
        /// complex matrix has a complex determinant and this decl is not that
        /// type; it is here for a caller that needs the component format.
        pub const element_component = component;
        /// The number of rows. `at` and `set` index it first, `data` second.
        pub const rows: usize = n_rows;
        /// The number of columns. `data` indexes it first, `at` and `set`
        /// second.
        pub const columns: usize = n_columns;
        /// The vector type `mulVec` returns, and it is ROW-sized: `Vec3f` for a
        /// three-row `f32` matrix, and `[rows]Element` wherever no named vector
        /// type matches, which is every complex matrix and every row count
        /// outside two to four.
        pub const vector_type = lane;

        /// The entries, COLUMN-MAJOR: `data[column][row]`.
        ///
        /// Public because the contract states the layout as the type's
        /// definition, and because a caller holding column-major data should be
        /// able to store it without a transposing loop.
        data: [n_columns][n_rows]Scalar,

        /// Builds a matrix from a ROW-MAJOR literal.
        ///
        /// `literal[row][column]` is the entry at ROW `row`, COLUMN `column`, so
        /// the initializer reads exactly the way the matrix is written down:
        ///
        ///     Mat2f.init(.{ .{ a, b },
        ///                   .{ c, d } })
        ///
        /// puts `a` at (0, 0), `b` at (0, 1), `c` at (1, 0), and `d` at (1, 1).
        ///
        /// The literal order is the OPPOSITE of `data`, which is column-major.
        /// That asymmetry is the point of the two orders: a literal is read by a
        /// human row by row, and storage is walked column by column during a
        /// contraction.
        pub fn init(literal: [n_rows][n_columns]Scalar) Self {
            var out: Self = undefined;
            inline for (0..n_rows) |row| {
                inline for (0..n_columns) |column| {
                    out.data[column][row] = literal[row][column];
                }
            }
            return out;
        }

        /// A matrix with every entry set to `value`.
        pub fn fill(value: Scalar) Self {
            var out: Self = undefined;
            inline for (0..n_columns) |column| {
                inline for (0..n_rows) |row| {
                    out.data[column][row] = value;
                }
            }
            return out;
        }

        /// The identity: `at(row, column)` is one when the indices agree and
        /// zero otherwise.
        ///
        /// It is available for a NON-SQUARE instantiation too, and there it
        /// means exactly that: `Mat(2, 3, f32).identity()` is the 2 by 3 matrix
        /// with a one at (0, 0) and at (1, 1) and a zero elsewhere. The
        /// rectangular identity is a left and right unit for every shape that
        /// chains with it, so there is nothing for a square-only guard to
        /// protect. The one thing it does not give is bit-exact equality back
        /// through a product, because the leading zero terms of a contraction
        /// cannot carry a negative zero; see the module comment.
        pub fn identity() Self {
            var out: Self = undefined;
            inline for (0..n_columns) |column| {
                inline for (0..n_rows) |row| {
                    out.data[column][row] = if (row == column) oneOf(Scalar) else zeroOf(Scalar);
                }
            }
            return out;
        }

        /// The entry at `row` and `column`, that is `data[column][row]`.
        ///
        /// Both indices are runtime values, so the bound is checked here and an
        /// index outside the matrix PANICS in every optimize mode, ReleaseFast
        /// and ReleaseSmall included: this is a `std.debug.panic`, not an
        /// `assert`, so the check is never compiled away. A compile error is
        /// impossible for a runtime index, and the message names the bound that
        /// was violated.
        pub fn at(self: Self, row: usize, column: usize) Scalar {
            if (row >= n_rows) {
                std.debug.panic(
                    "zmath.matrix.at: row {d} is out of range, the matrix has {d} rows",
                    .{ row, n_rows },
                );
            }
            if (column >= n_columns) {
                std.debug.panic(
                    "zmath.matrix.at: column {d} is out of range, the matrix has {d} columns",
                    .{ column, n_columns },
                );
            }
            return self.data[column][row];
        }

        /// Writes `value` at `row` and `column`, checked exactly as `at` is: an
        /// index outside the matrix panics in every optimize mode, and nothing is
        /// written on that path.
        pub fn set(self: *Self, row: usize, column: usize, value: Scalar) void {
            if (row >= n_rows) {
                std.debug.panic(
                    "zmath.matrix.set: row {d} is out of range, the matrix has {d} rows",
                    .{ row, n_rows },
                );
            }
            if (column >= n_columns) {
                std.debug.panic(
                    "zmath.matrix.set: column {d} is out of range, the matrix has {d} columns",
                    .{ column, n_columns },
                );
            }
            self.data[column][row] = value;
        }

        /// The entrywise sum.
        pub fn add(self: Self, other: Self) Self {
            @setFloatMode(.strict);
            var out: Self = undefined;
            inline for (0..n_columns) |column| {
                inline for (0..n_rows) |row| {
                    out.data[column][row] =
                        entryAdd(Scalar, self.data[column][row], other.data[column][row]);
                }
            }
            return out;
        }

        /// The entrywise difference.
        pub fn sub(self: Self, other: Self) Self {
            @setFloatMode(.strict);
            var out: Self = undefined;
            inline for (0..n_columns) |column| {
                inline for (0..n_rows) |row| {
                    out.data[column][row] =
                        entrySub(Scalar, self.data[column][row], other.data[column][row]);
                }
            }
            return out;
        }

        /// The entrywise negation, which maps `-0` to `+0` for a real entry and
        /// negates BOTH components of a complex one. A complex negation is not
        /// the conjugate: `-(1 + 2i)` is `-1 - 2i`, not `1 - 2i`.
        pub fn neg(self: Self) Self {
            @setFloatMode(.strict);
            var out: Self = undefined;
            inline for (0..n_columns) |column| {
                inline for (0..n_rows) |row| {
                    out.data[column][row] = entryNeg(Scalar, self.data[column][row]);
                }
            }
            return out;
        }

        /// Every entry times `factor`.
        pub fn scale(self: Self, factor: Scalar) Self {
            @setFloatMode(.strict);
            var out: Self = undefined;
            inline for (0..n_columns) |column| {
                inline for (0..n_rows) |row| {
                    out.data[column][row] = entryMul(Scalar, self.data[column][row], factor);
                }
            }
            return out;
        }

        /// The transpose, a `Mat(columns, rows, Element)`.
        pub fn transpose(self: Self) Mat(n_columns, n_rows, Scalar) {
            const Out = Mat(n_columns, n_rows, Scalar);
            var out: Out = undefined;
            inline for (0..n_columns) |column| {
                inline for (0..n_rows) |row| {
                    out.data[row][column] = self.data[column][row];
                }
            }
            return out;
        }

        /// Applies the matrix to a COLUMN vector with `columns` components,
        /// returning a `rows`-component vector of `vector_type`.
        ///
        /// The argument type is fixed at comptime and is NOT `vector_type` when
        /// the matrix is not square: it is the named vector `Vec2f`, `Vec3f`, or
        /// `Vec4f`, and the `d` counterparts, when the matrix is real and has
        /// exactly that many columns, and the plain array `[columns]Element` in
        /// every other case, which is every `Complexf` and `Complexd` matrix and
        /// every column count outside two to four. The result stays row-sized,
        /// so `Mat(2, 4, f32)` takes a `Vec4f` and returns a `Vec2f` while
        /// `Mat(4, 2, f32)` takes a `Vec2f` and returns a `Vec4f`. A vector of
        /// the wrong width or the wrong entry type does not convert: that call
        /// is a compile error.
        ///
        /// Each row of the result sums its products in increasing contraction
        /// index, starting from the index-zero product, and strictly: no
        /// reassociation and no fused multiply-add. See the module comment for
        /// what that costs at a negative zero.
        pub fn mulVec(self: Self, v: input_lane) vector_type {
            @setFloatMode(.strict);
            var lanes: [n_columns]Scalar = undefined;
            inline for (0..n_columns) |index| lanes[index] = laneAt(input_lane, v, index);
            var out: vector_type = undefined;
            inline for (0..n_rows) |row| {
                var sum = entryMul(Scalar, self.data[0][row], lanes[0]);
                for (1..n_columns) |column| {
                    sum = entryAdd(
                        Scalar,
                        sum,
                        entryMul(Scalar, self.data[column][row], lanes[column]),
                    );
                }
                laneSet(vector_type, &out, row, sum);
            }
            return out;
        }

        /// True when every entry is finite.
        pub fn isFinite(self: Self) bool {
            for (self.data) |column| {
                for (column) |entry| {
                    if (!isFiniteEntry(Scalar, entry)) return false;
                }
            }
            return true;
        }

        /// The determinant, by cofactor expansion along the first row.
        ///
        /// The result is a real number for a real matrix and a complex one for a
        /// complex matrix, so this is the only way to reach a complex
        /// determinant. A non-square instantiation is a compile error, not a
        /// zero and not a `null`. The named `determinantMat*` functions are the
        /// `f32` and `f64` surface the contract asks for, and a complex-suffixed
        /// name would be a naming convention `API.md#scope` does not have.
        pub fn determinant(self: Self) Element {
            if (comptime n_rows != n_columns) {
                @compileError(
                    "zmath.matrix.determinant: a determinant needs a square matrix, this one is " ++
                        std.fmt.comptimePrint("{d} by {d}", .{ n_rows, n_columns }),
                );
            }
            return cofactorDet(Self, &self);
        }
    };
}

/// The entry of `value` at `row` and `column`, without a bounds check.
///
/// Every kernel below indexes a pair its own loop bounds have already proved,
/// so the check in `at` would be dead work repeated throughout a contraction.
/// The matrix arrives by pointer because these kernels are generic over an
/// UNBOUNDED size: by value, every one of them would copy the whole matrix into
/// its own frame, once per entry read, and a `Mat` is a comptime-bounded array
/// the caller is not paying to have copied.
fn entryAt(comptime M: type, value: *const M, row: usize, column: usize) M.Element {
    return value.data[column][row];
}

/// `value` with `dropped_row` and `dropped_column` removed, a
/// `Mat(M.rows - 1, M.columns - 1, M.Element)`.
///
/// The copy is straight-line and indexed by comptime constants, which is what
/// lets the scalar replacement pass keep a minor in registers instead of
/// writing it out and reading it back. A pair of runtime counters stepping
/// through the surviving rows and columns computes the same addresses but
/// gives that pass nothing to fold.
fn minorOf(
    comptime M: type,
    value: *const M,
    comptime dropped_row: usize,
    comptime dropped_column: usize,
) Mat(M.rows - 1, M.columns - 1, M.Element) {
    if (comptime dropped_row >= M.rows or dropped_column >= M.columns) {
        @compileError(
            "zmath.matrix: a minor drops a row and a column of the matrix it is taken from",
        );
    }
    const Out = Mat(M.rows - 1, M.columns - 1, M.Element);
    var out: Out = undefined;
    inline for (0..(M.rows - 1)) |target_row| {
        inline for (0..(M.columns - 1)) |target_column| {
            // Where a copy lands maps back to the parent by stepping over the
            // dropped index, and both indices are comptime.
            const row = target_row + (if (target_row >= dropped_row) 1 else 0);
            const column = target_column + (if (target_column >= dropped_column) 1 else 0);
            out.data[target_column][target_row] = entryAt(M, value, row, column);
        }
    }
    return out;
}

/// The determinant of a 2 by 2, with the two products formed in this order.
fn det2(comptime M: type, value: *const M) M.Element {
    @setFloatMode(.strict);
    const S = M.Element;
    const leading = entryMul(S, entryAt(M, value, 0, 0), entryAt(M, value, 1, 1));
    const trailing = entryMul(S, entryAt(M, value, 0, 1), entryAt(M, value, 1, 0));
    return entrySub(S, leading, trailing);
}

/// The determinant of a 3 by 3, in closed form.
///
/// The three parenthesized differences are the 2 by 2 minors of the first row
/// taken in turn, so this is that row's cofactor expansion with its own
/// sub-determinants already evaluated. Stopping here is what keeps the
/// recursion in `cofactorDet` shallow: the size-4 and size-5 forms recurse only
/// as far as this.
fn det3(comptime M: type, value: *const M) M.Element {
    @setFloatMode(.strict);
    const S = M.Element;
    const m0 = entrySub(
        S,
        entryMul(S, entryAt(M, value, 1, 1), entryAt(M, value, 2, 2)),
        entryMul(S, entryAt(M, value, 1, 2), entryAt(M, value, 2, 1)),
    );
    const m1 = entrySub(
        S,
        entryMul(S, entryAt(M, value, 1, 0), entryAt(M, value, 2, 2)),
        entryMul(S, entryAt(M, value, 1, 2), entryAt(M, value, 2, 0)),
    );
    const m2 = entrySub(
        S,
        entryMul(S, entryAt(M, value, 1, 0), entryAt(M, value, 2, 1)),
        entryMul(S, entryAt(M, value, 1, 1), entryAt(M, value, 2, 0)),
    );
    const t0 = entryMul(S, entryAt(M, value, 0, 0), m0);
    const t1 = entryMul(S, entryAt(M, value, 0, 1), m1);
    const t2 = entryMul(S, entryAt(M, value, 0, 2), m2);
    // The first-row sign vector is (+1, -1, +1).
    return entryAdd(S, t0, entryAdd(S, entryNeg(S, t1), t2));
}

/// The determinant of a square `M`, by cofactor expansion along the first row.
///
/// Cofactors rather than Gaussian elimination, and the reason is in the module
/// comment: elimination's result depends on a pivot order, this does not. The
/// recursion is one function per size and the minors are values, so the number
/// of instantiations grows with the SIZES rather than with the minors.
fn cofactorDet(comptime M: type, value: *const M) M.Element {
    if (comptime M.rows != M.columns) {
        @compileError("zmath.matrix: a cofactor expansion needs a square matrix");
    }
    @setFloatMode(.strict);
    const S = M.Element;
    if (comptime M.rows == 1) return entryAt(M, value, 0, 0);
    if (comptime M.rows == 2) return det2(M, value);
    if (comptime M.rows == 3) return det3(M, value);
    const Minor = Mat(M.rows - 1, M.columns - 1, S);
    var sum: S = undefined;
    inline for (0..M.rows) |column| {
        const minor = minorOf(M, value, 0, column);
        const term = entryMul(S, entryAt(M, value, 0, column), cofactorDet(Minor, &minor));
        // The cofactor sign alternates along the row, starting positive.
        const signed = if (column % 2 == 1) entryNeg(S, term) else term;
        sum = if (column == 0) signed else entryAdd(S, sum, signed);
    }
    return sum;
}

// ---------------------------------------------------------------------------
// Inversion
// ---------------------------------------------------------------------------

/// The smallest pivot, relative to the largest entry magnitude of the matrix,
/// that a binary32 inverse will divide by.
///
/// This is the round-to-nearest unit roundoff of the format, `2^-24`, times the
/// `2^(n-1)` element-growth bound of a Gauss-Jordan reduction. For the largest
/// size this module inverts, `n = 4`, that factor is 8, so the constant is
/// `2^-21`. A pivot at or below it has no significant bits left once the
/// reduction's own rounding is accounted for, and dividing by it would amplify
/// rounding noise instead of recovering an inverse.
///
/// The test is relative to the LARGEST ENTRY OF THE WHOLE MATRIX, read before
/// any elimination, and not to the largest entry of the pivot column. A column
/// magnitude can itself collapse to zero by the time that column's turn comes,
/// which would make a column-relative test vacuous exactly when it is needed
/// most. The global scale cannot collapse, because nothing has been eliminated
/// yet. A column-relative test would also need a per-column constant; this one
/// constant covers every step.
///
/// It is public so a caller can see the threshold this module applies, and
/// reason about whether it is right for a badly scaled matrix. It is not
/// adjustable per call: the function has no epsilon argument, and an implicit
/// epsilon in a signature is how a numeric threshold becomes invisible.
pub const inverse_pivot_threshold_f: f32 = 0x1p-21;

/// The binary64 counterpart of `inverse_pivot_threshold_f`: the unit roundoff
/// `2^-53` times the same growth bound of 8, which is `2^-50`. It is a
/// published constant and not a per-call option, exactly as its binary32
/// sibling is: no `inverseMat*d` entry point takes an epsilon.
pub const inverse_pivot_threshold_d: f64 = 0x1p-50;

fn pivotThresholdOf(comptime Scalar: type) Component(Scalar) {
    return switch (Scalar) {
        f32 => inverse_pivot_threshold_f,
        f64 => inverse_pivot_threshold_d,
        else => @compileError("zmath.matrix: only a real matrix has a pivot magnitude"),
    };
}

/// The inverse of a square real `M` by Gauss-Jordan elimination with partial
/// pivoting, or `null`.
///
/// `null` means one of three things: an entry was not finite, a pivot fell at
/// or below `inverse_pivot_threshold_f` or `inverse_pivot_threshold_d` times the
/// largest entry magnitude of the matrix, or the reduction left a non-finite
/// entry. The finiteness check runs over every entry BEFORE any arithmetic,
/// because a step on a non-finite entry cannot be undone by a check afterwards.
///
/// The reduction below writes only the columns a later step or the result can
/// still read, and that is not a value-dependent claim. Column `j` of the left
/// block is read by the pivot search of step `j` and by nothing else: step `k`
/// reads column `k` for the pivot search and for the multiplier, and writes
/// every column from `k + 1` onward. So after step `j` no cell of column `j` is
/// read again, whatever the entries are, and the result loop reads only the
/// right block. An update skipped for a settled column would have written a
/// cell nothing reads, including a NaN an overflow would otherwise have left
/// behind.
fn gaussJordan(comptime M: type, value: M) ?M {
    if (comptime M.rows != M.columns or M.rows < 2) {
        @compileError("zmath.matrix: an inverse needs a square matrix of at least two rows");
    }
    @setFloatMode(.strict);
    const n = comptime M.rows;
    const S = M.Element;
    const C = Component(S);
    const threshold = pivotThresholdOf(S);

    var scale: C = 0;
    inline for (0..n) |row| {
        inline for (0..n) |column| {
            const entry = entryAt(M, &value, row, column);
            if (!std.math.isFinite(entry)) return null;
            const magnitude = @abs(entry);
            if (magnitude > scale) scale = magnitude;
        }
    }

    // [ value | identity ], held row-major so the reduction reads left to right.
    // The extent is comptime, so nothing is allocated.
    var work: [n][2 * n]C = undefined;
    inline for (0..n) |row| {
        inline for (0..n) |column| work[row][column] = entryAt(M, &value, row, column);
        inline for (0..n) |column| work[row][n + column] = if (row == column) 1 else 0;
    }

    inline for (0..n) |k| {
        // Partial pivoting: the largest magnitude in column `k` among the rows
        // that have not yet supplied a pivot row.
        var pivot_row: usize = k;
        var pivot_magnitude: C = @abs(work[k][k]);
        for (k + 1..n) |row| {
            const magnitude = @abs(work[row][k]);
            if (magnitude > pivot_magnitude) {
                pivot_magnitude = magnitude;
                pivot_row = row;
            }
        }
        if (!(pivot_magnitude > threshold * scale)) return null;
        if (pivot_row != k) {
            // The interchange moves the live part of the two rows. Columns
            // before `k` are settled in EVERY row, so moving them changes
            // nothing that is read again, and they are not read again.
            inline for (k..2 * n) |j| {
                const held = work[k][j];
                work[k][j] = work[pivot_row][j];
                work[pivot_row][j] = held;
            }
        }

        const pivot = work[k][k];
        // Columns `0..k` of the left block stop here, for the reason the
        // module comment on the function gives. Column `k` has just been read
        // as `factor` by every row that will use it, and columns `0..k-1` were
        // read by their own step and by nothing since, so normalizing them
        // again would divide into cells that are never read.
        inline for (k + 1..2 * n) |j| work[k][j] = work[k][j] / pivot;
        for (0..n) |row| {
            if (row == k) continue;
            const factor = work[row][k];
            if (factor == 0) continue;
            // Eliminating the rows above as well as below is what makes this
            // Gauss-Jordan and leaves the identity in the left block. The
            // columns written here are the ones the later steps and the result
            // read; the rest are settled.
            inline for (k + 1..2 * n) |j| {
                work[row][j] = entrySub(S, work[row][j], entryMul(S, factor, work[k][j]));
            }
        }
    }

    var out: M = undefined;
    inline for (0..n) |row| {
        inline for (0..n) |column| {
            const entry = work[row][n + column];
            // A near-singular matrix can pass the pivot test and still drive an
            // entry of the inverse past the format's range. Reporting that as
            // "no inverse" is the honest answer; returning an infinity is not.
            if (!std.math.isFinite(entry)) return null;
            out.data[column][row] = entry;
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Products
// ---------------------------------------------------------------------------

/// The product `left * right`, a `Mat(left.rows, right.columns, Element)`.
///
/// `(A * B) * v` applies `B` first and `A` second, so composition reads
/// right-to-left. The operands fix the row count, the inner dimension, and the
/// entry type, which leaves the output width as the only parameter a return
/// type cannot supply; that is the shape `API.md#matrices` specifies. Each entry
/// sums its products in increasing inner index, starting from the index-zero
/// product, and strictly, so it is neither reassociated nor contracted into an
/// FMA.
///
/// Everything about that shape is a compile-time obligation, and each violation
/// is a compile error naming it: an operand that is not a `Mat`, a right
/// operand whose row count is not the left operand's column count, an output
/// width that is not the right operand's column count, and two operands whose
/// entry types differ. None of those paths computes or mutates anything.
pub fn mul(
    comptime n_out_columns: usize,
    left: anytype,
    right: anytype,
) Product(@TypeOf(left), @TypeOf(right), n_out_columns) {
    @setFloatMode(.strict);
    const L = @TypeOf(left);
    const R = @TypeOf(right);
    const S = L.Element;
    const n = comptime L.rows;
    const inner = comptime L.columns;
    const Out = Product(L, R, n_out_columns);
    var out: Out = undefined;
    inline for (0..n_out_columns) |column| {
        inline for (0..n) |row| {
            var sum = entryMul(S, entryAt(L, &left, row, 0), entryAt(R, &right, 0, column));
            for (1..inner) |index| {
                sum = entryAdd(
                    S,
                    sum,
                    entryMul(S, entryAt(L, &left, row, index), entryAt(R, &right, index, column)),
                );
            }
            out.data[column][row] = sum;
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Named operations
//
// `API.md#matrices` puts these on the module rather than on the type, and that
// is where they are. The size-paired siblings the contract names only three of
// are here too, because `zmath.linalg` consumes the whole family and a caller
// should not have to reach for `Mat` and recompute a determinant by hand.
// ---------------------------------------------------------------------------

/// The determinant of a 2 by 2 `f32` matrix.
pub fn determinantMat2f(value: Mat2f) f32 {
    return value.determinant();
}

/// The determinant of a 2 by 2 `f64` matrix.
pub fn determinantMat2d(value: Mat2d) f64 {
    return value.determinant();
}

/// The determinant of a 3 by 3 `f32` matrix.
pub fn determinantMat3f(value: Mat3f) f32 {
    return value.determinant();
}

/// The determinant of a 3 by 3 `f64` matrix.
pub fn determinantMat3d(value: Mat3d) f64 {
    return value.determinant();
}

/// The determinant of a 4 by 4 `f32` matrix.
pub fn determinantMat4f(value: Mat4f) f32 {
    return value.determinant();
}

/// The determinant of a 4 by 4 `f64` matrix.
pub fn determinantMat4d(value: Mat4d) f64 {
    return value.determinant();
}

/// The determinant of a 5 by 5 `f32` matrix.
pub fn determinantMat5f(value: Mat5f) f32 {
    return value.determinant();
}

/// The determinant of a 5 by 5 `f64` matrix.
pub fn determinantMat5d(value: Mat5d) f64 {
    return value.determinant();
}

/// The inverse of a 2 by 2 `f32` matrix, or `null`.
///
/// `null` means one of three things: an entry was not finite, a pivot fell at
/// or below `inverse_pivot_threshold_f` times the largest entry magnitude of the
/// matrix, or the reduction left a non-finite entry. The threshold is a
/// published constant and not a per-call option; this entry point takes no
/// epsilon. The result is a fresh value returned by value: `value` is not
/// modified, and a `null` leaves no partial result anywhere, because the
/// reduction holds `[ value | identity ]` in its own comptime-bounded work
/// array.
pub fn inverseMat2f(value: Mat2f) ?Mat2f {
    return gaussJordan(Mat2f, value);
}

/// The inverse of a 2 by 2 `f64` matrix, or `null`.
///
/// `null` means one of three things: an entry was not finite, a pivot fell at
/// or below `inverse_pivot_threshold_d` times the largest entry magnitude of the
/// matrix, or the reduction left a non-finite entry. The result is a fresh value
/// returned by value, `value` is not modified, and a `null` leaves no partial
/// result. See `inverseMat2f` for the same contract in binary32.
pub fn inverseMat2d(value: Mat2d) ?Mat2d {
    return gaussJordan(Mat2d, value);
}

/// The inverse of a 3 by 3 `f32` matrix, or `null`.
///
/// `null` means one of three things: an entry was not finite, a pivot fell at
/// or below `inverse_pivot_threshold_f` times the largest entry magnitude of the
/// matrix, or the reduction left a non-finite entry. The result is a fresh value
/// returned by value, `value` is not modified, and a `null` leaves no partial
/// result.
pub fn inverseMat3f(value: Mat3f) ?Mat3f {
    return gaussJordan(Mat3f, value);
}

/// The inverse of a 3 by 3 `f64` matrix, or `null`.
///
/// `null` means one of three things: an entry was not finite, a pivot fell at
/// or below `inverse_pivot_threshold_d` times the largest entry magnitude of the
/// matrix, or the reduction left a non-finite entry. The result is a fresh value
/// returned by value, `value` is not modified, and a `null` leaves no partial
/// result.
pub fn inverseMat3d(value: Mat3d) ?Mat3d {
    return gaussJordan(Mat3d, value);
}

/// The inverse of a 4 by 4 `f32` matrix, or `null`.
///
/// `null` means one of three things: an entry was not finite, a pivot fell at
/// or below `inverse_pivot_threshold_f` times the largest entry magnitude of the
/// matrix, or the reduction left a non-finite entry. The result is a fresh value
/// returned by value, `value` is not modified, and a `null` leaves no partial
/// result.
///
/// The same conditioning caveat as `inverseMat4d` with a coarser threshold:
/// binary32 has a unit roundoff of `2^-24`, so the usable precision of a
/// well-conditioned result is roughly seven decimal digits, and nothing here
/// repairs an ill-conditioned input.
pub fn inverseMat4f(value: Mat4f) ?Mat4f {
    return gaussJordan(Mat4f, value);
}

/// The inverse of a 4 by 4 `f64` matrix, or `null`.
///
/// `null` means one of three things: an entry was not finite, a pivot fell at
/// or below `inverse_pivot_threshold_d` times the largest entry magnitude of the
/// matrix, or the reduction left a non-finite entry. The threshold is a
/// published constant and not a per-call option; this entry point takes no
/// epsilon. The result is a fresh value returned by value: `value` is not
/// modified, and a `null` leaves no partial result anywhere, because the
/// reduction holds `[ value | identity ]` in its own comptime-bounded work
/// array.
///
/// The result is exact only for a well-conditioned input: an ill-conditioned
/// matrix gives an inverse whose error is governed by its condition number,
/// not by this module.
pub fn inverseMat4d(value: Mat4d) ?Mat4d {
    return gaussJordan(Mat4d, value);
}

/// Applies a `Mat4f` to a `Vec4f` column vector.
pub fn mulMat4fVec4f(matrix: Mat4f, v: Vec4f) Vec4f {
    return matrix.mulVec(v);
}

/// Applies a `Mat4d` to a `Vec4d` column vector.
pub fn mulMat4dVec4d(matrix: Mat4d, v: Vec4d) Vec4d {
    return matrix.mulVec(v);
}

// ---------------------------------------------------------------------------
// Named sizes
// ---------------------------------------------------------------------------

/// A 2 by 2 matrix in binary32.
pub const Mat2f = Mat(2, 2, f32);
/// A 3 by 3 matrix in binary32.
pub const Mat3f = Mat(3, 3, f32);
/// A 4 by 4 matrix in binary32: a homogeneous 3D transform or a linear 4D
/// transform.
pub const Mat4f = Mat(4, 4, f32);
/// A 5 by 5 matrix in binary32: a homogeneous 4D affine or projective
/// transform.
pub const Mat5f = Mat(5, 5, f32);
/// A 2 by 2 matrix in binary64.
pub const Mat2d = Mat(2, 2, f64);
/// A 3 by 3 matrix in binary64.
pub const Mat3d = Mat(3, 3, f64);
/// A 4 by 4 matrix in binary64: a homogeneous 3D transform or a linear 4D
/// transform.
pub const Mat4d = Mat(4, 4, f64);
/// A 5 by 5 matrix in binary64: a homogeneous 4D affine or projective
/// transform.
pub const Mat5d = Mat(5, 5, f64);
