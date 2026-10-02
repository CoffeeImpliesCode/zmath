/* Correctly-rounded true gamma function for binary16 value.

Copyright (c) 2023-2025 Alexei Sibidanov and Paul Zimmermann

This file is part of the CORE-MATH project
(https://core-math.gitlabpages.inria.fr/).

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/* This code is based on the binary32 code tgammaf.c:
   the kernel receives binary16 bit patterns and decodes them to binary32
   (exactly),
   we then use the same code than for tgammaf,
   and at output we round once, directly to binary16.
   The changes with respect to tgammaf.c are marked with a comment
   "specific tgammaf16". */

#include <stdint.h>

#include "coremath_half.h"

/* Warning: clang also defines __GNUC__. */
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

/* Round x to nearest integer, breaking ties to even.
   The finite inputs used by this kernel are small enough for an exact
   conversion to int64_t. */
static inline double
roundeven_finite(double x)
{
  int64_t truncated = (int64_t)x;
  double truncated_value = (double)truncated;
  double fraction = x - truncated_value;
  if (fraction > 0.5 ||
      (fraction == 0.5 && (truncated & 1) != 0))
    return truncated_value + 1.0;
  if (fraction < -0.5 ||
      (fraction == -0.5 && (truncated & 1) != 0))
    return truncated_value - 1.0;
  return truncated_value;
}

static inline float
floor_finite(float x)
{
  int32_t truncated = (int32_t)x;
  if ((float)truncated > x)
    --truncated;
  return (float)truncated;
}

typedef union {
  float f;
  uint32_t u;
} b32u32_u;

/* specific tgammaf16: the input is the binary16 bit pattern of x */
uint16_t
zmath_core_tgammah_bits(uint16_t x_bits)
{
  float x = zmath_half_to_float_bits(x_bits);
  b32u32_u t = {.f = x};
  uint32_t ax = t.u << 1;
  if (__builtin_expect(ax >= (0xffu << 24), 0)) { /* x=NaN or +/-Inf */
    if (ax == (0xffu << 24)) {                   /* x=+/-Inf */
      if (t.u >> 31)                             /* x=-Inf */
        return ZMATH_HALF_CANONICAL_NAN_BITS;      /* x=-Inf */
      return ZMATH_HALF_INFINITY_BITS;              /* x=+Inf */
    }
    return zmath_half_nan_bits(x_bits); /* Quiet an sNaN and raise invalid operation. */
  }
  double z = x;
  /* specific tgammaf16: changed threshold to deal with all overflow cases. */
  if (__builtin_expect(ax <= 0x6f000000u, 0)) { /* |x| <= 0x1p-16 */
    volatile double d =
        (0x1.fa658c23b1578p-1 - 0x1.d0a118f324b63p-1 * z) * z -
        0x1.2788cfc6fb619p-1;
    double f = 1.0 / z + d;
    return zmath_half_from_double_bits(f);
  }
  float fx = floor_finite(x);
  if (__builtin_expect(x >= 0x1.274p+3f, 0))
    return zmath_half_from_double_bits((double)(0x1p127f * 0x1p127f));
  int32_t k = fx; /* specific tgammaf16: no overflow possible */
  if (__builtin_expect(fx == x, 0)) { /* x is integer */
    if (x == 0.0f)
      return ZMATH_HALF_INFINITY_BITS;
    if (x < 0.0f)
      return ZMATH_HALF_CANONICAL_NAN_BITS; /* Raise invalid operation. */
    double t0 = 1, x0 = 1;
    for (int32_t i = 1; i < k; i++, x0 += 1.0)
      t0 *= x0;
    return zmath_half_from_double_bits(t0);
  }
  /* specific tgammaf16: changed threshold of -42 to -16 */
  if (__builtin_expect(x < -16.0f, 0)) { /* negative non-integer */
    /* For x < -16, x non-integer, |gamma(x)| < 2^-25. */
    static const float sgn[2] = {0x1p-127f, -0x1p-127f};
    return zmath_half_from_double_bits((double)(0x1p-127f * sgn[k & 1]));
  }
  static const double c[] = {
      0x1.c9a76be577123p+0,  0x1.8f2754ddcf90dp+0,
      0x1.0d1191949419bp+0,  0x1.e1f42cf0ae4a1p-2,
      0x1.82b358a3ab638p-3,  0x1.e1f2b30cd907bp-5,
      0x1.240f6d4071bd8p-6,  0x1.1522c9f3cd012p-8,
      0x1.1fd0051a0525bp-10, 0x1.9808a8b96c37ep-13,
      0x1.b3f78e01152b5p-15, 0x1.49c85a7e1fd04p-18,
      0x1.471ca49184475p-19, -0x1.368f0b7ed9e36p-23,
      0x1.882222f9049efp-23, -0x1.a69ed2042842cp-25,
  };

  double m = z - 0x1.7p+1;
  double i = roundeven_finite(m);
  double step = __builtin_copysign(1.0, i);
  double d = m - i, d2 = d * d, d4 = d2 * d2, d8 = d4 * d4;
  double f = (c[0] + d * c[1]) + d2 * (c[2] + d * c[3]) +
             d4 * ((c[4] + d * c[5]) + d2 * (c[6] + d * c[7])) +
             d8 * ((c[8] + d * c[9]) + d2 * (c[10] + d * c[11]) +
                   d4 * ((c[12] + d * c[13]) +
                         d2 * (c[14] + d * c[15])));
  int jm = __builtin_fabs(i);
  double w = 1;
  if (jm) {
    z -= 0.5 + step * 0.5;
    w = z;
    for (int j = jm - 1; j; j--) {
      z -= step;
      w *= z;
    }
  }
  if (i <= -0.5)
    w = 1 / w;
  f *= w;
  return zmath_half_from_double_bits(f);
}
