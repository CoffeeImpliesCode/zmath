/* Private library-free integer rounding for the zmath C kernels.

   Keep the compiler-native operation where LLVM lowers it without a C23 libm
   call; elsewhere round the binary64 encoding with integer arithmetic, so no
   target needs the compiler runtime's round, roundeven or floor. A target whose
   backend has no scalar rounding instruction (RISC-V lowers both round and
   floor to a call into libm) therefore links without libm and still rounds
   nearest-even on the ties. */

#ifndef ZMATH_COREMATH_ROUNDEVEN_H
#define ZMATH_COREMATH_ROUNDEVEN_H

#include <stdint.h>

/* binary64 encoding landmarks */
#define ZMATH_DBL_SIGN_BIT UINT64_C(0x8000000000000000)
#define ZMATH_DBL_MANTISSA_MASK UINT64_C(0x000fffffffffffff)
#define ZMATH_DBL_HIDDEN_BIT UINT64_C(0x0010000000000000)
#define ZMATH_DBL_HALF_BITS UINT64_C(0x3fe0000000000000)
#define ZMATH_DBL_ONE_BITS UINT64_C(0x3ff0000000000000)
#define ZMATH_DBL_TWO_POW_52_BITS UINT64_C(0x4330000000000000)

/* Reassemble the integer (significand >> shift) * 2^shift as a binary64 with
   the sign of the argument, where 1 <= |x| < 2^52 and shift is 1075 minus the
   biased exponent. A significand that reaches 2^(53-shift) is a power of two
   and carries into the exponent. */
static inline double zmath_dbl_from_whole(uint64_t sign, uint64_t exponent,
                                          uint64_t whole, uint32_t shift) {
  union { double f; uint64_t u; } out = { .u = 0 };
  if (whole == (UINT64_C(0x0020000000000000) >> shift)) {
    whole >>= 1;
    exponent++;
  }
  out.u = sign | (exponent << 52) |
          ((whole << shift) & ZMATH_DBL_MANTISSA_MASK);
  return out.f;
}

/* Nearest-even integer of a finite binary64 kernel term, with the sign of the
   argument: a tie goes to the even neighbour and a value below one keeps the
   signed zero. An infinity, a NaN and an already integral value are returned
   unchanged, as the compiler's own rounding operations do. */
static inline double zmath_roundeven_finite(double x) {
#if defined(__SSE4_1__) || \
    (defined(__aarch64__) && defined(__ARM_FP)) || \
    (defined(__ARM_ARCH) && __ARM_ARCH >= 8 && defined(__ARM_FP) && (__ARM_FP & 8)) || \
    defined(__wasm__) || \
    (defined(__loongarch_sx) && defined(__loongarch_frlen) && __loongarch_frlen >= 64)
  return __builtin_roundeven(x);
#else
  union { double f; uint64_t u; } in = { .f = x };
  const uint64_t sign = in.u & ZMATH_DBL_SIGN_BIT;
  const uint64_t magnitude = in.u & ~ZMATH_DBL_SIGN_BIT;
  uint64_t significand, exponent, whole, rest, half_way;
  uint32_t shift;

  /* |x| <= 1/2 is a signed zero and 1/2 < |x| < 1 is the signed unit. */
  if (magnitude <= ZMATH_DBL_HALF_BITS) {
    in.u = sign;
    return in.f;
  }
  if (magnitude < ZMATH_DBL_ONE_BITS) {
    in.u = sign | ZMATH_DBL_ONE_BITS;
    return in.f;
  }
  if (magnitude >= ZMATH_DBL_TWO_POW_52_BITS)
    return x;

  exponent = magnitude >> 52;
  /* Fractional bits below the unit: 52 minus the unbiased exponent, so 1 <=
     shift <= 52 for 1 <= |x| < 2^52. */
  shift = (uint32_t)(1075 - exponent);
  significand = ZMATH_DBL_HIDDEN_BIT | (magnitude & ZMATH_DBL_MANTISSA_MASK);
  whole = significand >> shift;
  rest = significand & ((UINT64_C(1) << shift) - 1);
  half_way = UINT64_C(1) << (shift - 1);
  /* Above the half the remainder rounds away, and a tie keeps the even
     integer: the parity is the lowest retained significand bit, never a count
     of trailing zero bits of an encoding. */
  if (rest > half_way || (rest == half_way && (whole & 1) != 0))
    whole++;
  return zmath_dbl_from_whole(sign, exponent, whole, shift);
#endif
}

/* Largest integral binary64 not greater than any argument, an infinity and a
   NaN included. The rounding is exact: a signed zero stays a signed zero, a
   value below one floors to -1 or to 0, and an integral value, an infinity
   and a NaN are returned unchanged. */
static inline double zmath_floor(double x) {
  union { double f; uint64_t u; } in = { .f = x };
  const uint64_t sign = in.u & ZMATH_DBL_SIGN_BIT;
  const uint64_t magnitude = in.u & ~ZMATH_DBL_SIGN_BIT;
  uint64_t significand, exponent, whole, rest;
  uint32_t shift;

  if (magnitude == 0 || magnitude >= ZMATH_DBL_TWO_POW_52_BITS)
    return x;
  if (magnitude < ZMATH_DBL_ONE_BITS) {
    /* -1 < x < 0 floors to -1; 0 < x < 1 floors to the positive zero. */
    in.u = sign ? sign | ZMATH_DBL_ONE_BITS : 0;
    return in.f;
  }

  exponent = magnitude >> 52;
  shift = (uint32_t)(1075 - exponent);
  significand = ZMATH_DBL_HIDDEN_BIT | (magnitude & ZMATH_DBL_MANTISSA_MASK);
  whole = significand >> shift;
  rest = significand & ((UINT64_C(1) << shift) - 1);
  /* Truncation is toward zero; a negative value with a discarded fraction
     takes the lower neighbour. */
  if (rest != 0 && sign != 0)
    whole++;
  return zmath_dbl_from_whole(sign, exponent, whole, shift);
}

#endif
