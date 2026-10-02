/* Private binary16 ABI helpers for the zmath C kernels.

   The binary16 kernels use a uint16_t ABI instead of a native C half type,
   so that mips, powerpc and wasm32 can build them. The kernels still run
   their accepted binary32/binary64 polynomial graphs: operands are decoded
   exactly into binary32 or binary64, and results are narrowed once, directly
   from binary64, with round-to-nearest ties-to-even. Routing a result through
   binary32 first would double-round it and change the accepted answers.

   Targets whose compiler models IEEE-754 binary16 as _Float16 use native
   conversions except on 32-bit ARM, where Zig's compiler-runtime half
   conversion ABI differs from Clang's lowering. Those targets use the portable
   integer decode and direct round-to-nearest-even pack below.

   Scratch storage is statically bounded. The integer fallback calls no
   library; native conversions use the selected compiler's usual lowering. */

#ifndef ZMATH_COREMATH_HALF_H
#define ZMATH_COREMATH_HALF_H

#include <stdint.h>

typedef union { float f; uint32_t u; } zmath_half_f32u32_u;
typedef union { double f; uint64_t u; } zmath_half_f64u64_u;

/* binary16 encoding landmarks */
#define ZMATH_HALF_SIGN_BIT 0x8000u
#define ZMATH_HALF_QUIET_BIT 0x0200u
#define ZMATH_HALF_MANTISSA_MASK 0x03ffu
#define ZMATH_HALF_INFINITY_BITS 0x7c00u
#define ZMATH_HALF_CANONICAL_NAN_BITS 0x7e00u

/* Feature guard: the compiler must model a 16-bit radix-2 format with an
   11-bit significand and the binary16 exponent range. This is a per-target
   capability probe, never a requirement. The ARM exclusion avoids mismatched
   core-register versus VFP arguments to the half conversion runtime. */
#if !defined(__arm__) && !defined(__thumb__) && \
    defined(__FLT_RADIX__) && defined(__FLT16_MANT_DIG__) && \
    defined(__FLT16_MAX_EXP__) && defined(__FLT16_MIN_EXP__) && \
    __FLT_RADIX__ == 2 && __FLT16_MANT_DIG__ == 11 && \
    __FLT16_MAX_EXP__ == 16 && __FLT16_MIN_EXP__ == -13
#define ZMATH_HALF_NATIVE_FLOAT16 1
typedef union { _Float16 f; uint16_t u; } zmath_half_native_u;
#endif

/* Result of NaN + NaN on a binary16 operand: the payload and the sign are
   kept and the quiet bit is set. This replaces the former `x + x` on a
   native half type, so it never depends on native half support. */
static inline uint16_t zmath_half_nan_bits(uint16_t bits) {
  return (uint16_t)(bits | ZMATH_HALF_QUIET_BIT);
}

/* Exact binary16 to binary32 widening. Every binary16 value is a binary32
   value, so no rounding happens here. */
static inline float zmath_half_to_float_bits(uint16_t bits) {
#ifdef ZMATH_HALF_NATIVE_FLOAT16
  zmath_half_native_u native;
  native.u = bits;
  return (float)native.f;
#else
  uint32_t sign = (uint32_t)(bits & ZMATH_HALF_SIGN_BIT) << 16;
  uint32_t exponent = ((uint32_t)bits >> 10) & 0x1fu;
  uint32_t fraction = (uint32_t)bits & ZMATH_HALF_MANTISSA_MASK;
  zmath_half_f32u32_u v;
  if (__builtin_expect(exponent == 0, 0)) {
    if (__builtin_expect(fraction == 0, 0)) {
      v.u = sign; /* signed zero */
    } else {
      /* A binary16 subnormal is 2^-24 * fraction; renormalize it so that the
         leading 1 lands in the binary32 hidden bit. */
      int shift = __builtin_clz(fraction) - 21;
      fraction <<= shift;
      v.u = sign | ((uint32_t)(113 - shift) << 23) |
            ((fraction & ZMATH_HALF_MANTISSA_MASK) << 13);
    }
  } else if (__builtin_expect(exponent == 0x1fu, 0)) {
    v.u = sign | 0x7f800000u | (fraction << 13); /* infinity or NaN */
  } else {
    v.u = sign | ((exponent + 112u) << 23) | (fraction << 13);
  }
  return v.f;
#endif
}

/* Exact binary16 to binary64 widening, for the kernels whose evaluation graph
   is binary64 (pow). */
static inline double zmath_half_to_double_bits(uint16_t bits) {
#ifdef ZMATH_HALF_NATIVE_FLOAT16
  zmath_half_native_u native;
  native.u = bits;
  return (double)native.f;
#else
  uint64_t sign = (uint64_t)(bits & ZMATH_HALF_SIGN_BIT) << 48;
  uint32_t exponent = ((uint32_t)bits >> 10) & 0x1fu;
  uint32_t fraction = (uint32_t)bits & ZMATH_HALF_MANTISSA_MASK;
  zmath_half_f64u64_u v;
  if (__builtin_expect(exponent == 0, 0)) {
    if (__builtin_expect(fraction == 0, 0)) {
      v.u = sign; /* signed zero */
    } else {
      /* A binary16 subnormal is 2^-24 * fraction with fraction < 2^10, so it
         is a binary64 normal whose leading bit sits at 52 - top. */
      int top = 31 - __builtin_clz(fraction);
      v.u = sign | ((uint64_t)(top - 24 + 1023) << 52) |
            (((uint64_t)fraction << (52 - top)) & 0x000fffffffffffffull);
    }
  } else if (__builtin_expect(exponent == 0x1fu, 0)) {
    v.u = sign | 0x7ff0000000000000ull | ((uint64_t)fraction << 42);
  } else {
    v.u = sign | ((uint64_t)(exponent + 1008u) << 52) | ((uint64_t)fraction << 42);
  }
  return v.f;
#endif
}

/* Direct round-to-nearest-even binary64 to binary16 narrowing.

   The binary64 significand keeps 53 bits and binary16 keeps 11, so a normal
   result drops exactly 42 bits; the guard bit is bit 41 and the sticky value
   is the OR of the 41 bits below it. Subnormal results additionally drop
   enough bits to land on the 2^-24 subnormal grid. This performs one
   rounding and is never routed through binary32.

   NaN inputs are packed in integer arithmetic on every target: LLVM models
   the exact binary16 <-> binary32/binary64 conversions as bit moves and can
   leave a signaling NaN signaling, and a NaN payload must never depend on
   host conversion rules. */
static inline uint16_t zmath_half_from_double_bits(double x) {
  zmath_half_f64u64_u v = {.f = x};
  uint64_t sign = (v.u >> 48) & ZMATH_HALF_SIGN_BIT;
  uint64_t exponent = (v.u >> 52) & 0x7ffu;
  uint64_t fraction = v.u & 0x000fffffffffffffull;

  if (__builtin_expect(exponent == 0x7ffu, 0)) {
    if (__builtin_expect(fraction == 0, 0))
      return (uint16_t)(sign | ZMATH_HALF_INFINITY_BITS);
    /* Narrowing a quiet NaN keeps the leading payload bits and stays quiet. */
    return (uint16_t)(sign | ZMATH_HALF_INFINITY_BITS |
                      ((fraction >> 42) & ZMATH_HALF_MANTISSA_MASK) |
                      ZMATH_HALF_QUIET_BIT);
  }

#ifdef ZMATH_HALF_NATIVE_FLOAT16
  zmath_half_native_u native;
  native.f = (_Float16)x;
  return native.u;
#else
  /* binary64 zeros and binary64 subnormals are far below the smallest
     binary16 subnormal, so both round to a signed zero. */
  if (__builtin_expect(exponent == 0, 0)) return (uint16_t)sign;

  int unbiased = (int)exponent - 1023;
  /* 2^16 is the first binary16 overflow; values below 2^-25 are under half
     the smallest binary16 subnormal, so they round to a signed zero. */
  if (__builtin_expect(unbiased > 15, 0))
    return (uint16_t)(sign | ZMATH_HALF_INFINITY_BITS);
  if (__builtin_expect(unbiased < -25, 0)) return (uint16_t)sign;

  uint64_t significand = (1ull << 52) | fraction;
  int shift = unbiased < -14 ? 28 - unbiased : 42;
  uint64_t rounded = significand >> shift;
  uint64_t guard = (significand >> (shift - 1)) & 1u;
  uint64_t sticky =
      (significand & ((((uint64_t)1) << (shift - 1)) - 1u)) != 0;
  if (guard && (sticky || (rounded & 1u))) rounded++;

  if (__builtin_expect(unbiased < -14, 0)) {
    /* A rounding carry out of the subnormal grid reaches 2^-14 exactly. */
    if (rounded >= 1024) return (uint16_t)(sign | (1u << 10));
    return (uint16_t)(sign | rounded);
  }
  if (__builtin_expect(rounded == 2048, 0)) {
    /* The rounding carry leaves the binary16 exponent range. */
    if (unbiased == 15) return (uint16_t)(sign | ZMATH_HALF_INFINITY_BITS);
    return (uint16_t)(sign | ((uint64_t)(unbiased + 16) << 10));
  }
  return (uint16_t)(sign | ((uint64_t)(unbiased + 15) << 10) |
                    (rounded & ZMATH_HALF_MANTISSA_MASK));
#endif
}

#endif