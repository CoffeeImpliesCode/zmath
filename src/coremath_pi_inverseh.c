/* Direct correctly-rounded inverse-pi kernels for binary16. */
#include <stdint.h>

/* Pinned CORE-MATH source: revision 3bafed75e3c22399531857d63a3fa0e24bfb42fa, src/binary16/asinpi/asinpif16.c. */
/* Correctly-rounded half-revolution arc-sine for binary16 value.

Copyright (c) 2025 Paul Zimmermann

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


// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif


typedef union {float f; uint32_t u;} zmath_pi_asinpih_b32u32_u;

// the following polynomials were generated using Sollya (cf asin.sollya)

/* degree-7 minimax polynomial for asin(x) over [0,0.25], with relative error
   bounded by 2^-26.989, with coefficients of odd degree only, and degree-1
   coefficient forced to 1
*/
static const float zmath_pi_asinpih_p0[] = {1.0f, 0x1.5555fp-3f, 0x1.32b088p-4f, 0x1.8baf46p-5f};

/* degree-7 minimax polynomial for asin(x) over [0.25,0.5], with relative error
   bounded by 2^-19.202, with coefficients of odd degree only, and degree-1
   coefficient forced to 1. This polynomial was optimized to reduce the
   number of exceptions.
*/
static const float zmath_pi_asinpih_p1[] = {0x1.000002p+0f, 0x1.55a7ep-3f, 0x1.258e54p-4f, 0x1.08e44p-4f};

#define ZMATH_PI_ASINPIH_HALF_PI 0x1.921fb4p+0f
#define ZMATH_PI_ASINPIH_HALF 0.5f
#define ZMATH_PI_ASINPIH_INV_PI 0x1.45f306p-2f

_Float16 zmath_core_asinpih (_Float16 x)
{
  zmath_pi_asinpih_b32u32_u v = {.f = x};
  uint32_t u = v.u;
  uint32_t au = u & 0x7fffffffu;

  if (au >= 0x3f800000u) { // NaN, Inf, or |x| >= 1
    if (au == 0x3f800000u) // |x| = 1
      return (u == 0x3f800000u) ? ZMATH_PI_ASINPIH_HALF : -ZMATH_PI_ASINPIH_HALF;
    if ((au >> 23) == 0x3ff && ((au & 0x7fffff) != 0)) // qNaN or sNaN
      return x + x;
    return 0.0f / 0.0f; // will signal invalid and return sNaN
  }

  float t = v.f, tt, c1, c3, c5, y;

  /* for |x| < 0.5 we don't reduce t into -t, since the polynomials
     zmath_pi_asinpih_p0 and zmath_pi_asinpih_p1 are valid also on the negative side */

  int reduce = au >= 0x3f000000u; // |x| >= 0.5
  if (reduce) {
    if (au == 0x3f10a000u) // |x|=0x1.214p-1
      return (au == u) ? 0x1.876016p-3f : -0x1.876016p-3f;
    if (au == 0x3f208000u) // |x|=0x1.41p-1
      return (au == u) ? 0x1.b9c002p-3f : -0x1.b9c002p-3f;
    if (au == 0x3f3ca000u) // |x|=0x1.794p-1
      return (au == u) ? 0x1.0dfffap-2f : -0x1.0dfffap-2f;
    if (u >> 31) // x < 0
      t = -t;
    // argument reduction: asin(x) = pi/2 - 2*asin(sqrt((1-x)/2))
    t = __builtin_sqrtf ((1.0f - t) * 0.5f);
    zmath_pi_asinpih_b32u32_u w = {.f = t};
    au = w.u & 0x7fffffffu;
  }

  // now 0 <= t <= 0.5
  tt = t * t;

  if (au < 0x3e800000u) { // 0 <= t < 0.25
    /* For |x| <= 0x1.92p-13, we get underflow, except when x=0x1.92p-13
       with RNDU and x=-0x1.92p-13 with RNDD. */
    if (au <= 0x39490000u) {
      // ZMATH_PI_ASINPIH_INV_PI is smaller than 1/pi
      return __builtin_fmaf (t, ZMATH_PI_ASINPIH_INV_PI, 0x1p-25f * t);
    }
    if (au == 0x3e002000) // |x| = 0x1.004p-3
      return (au == u) ? 0x1.472002p-5f : -0x1.472002p-5f;
    c5 = __builtin_fmaf (zmath_pi_asinpih_p0[3], tt, zmath_pi_asinpih_p0[2]); // degree 5
    c3 = __builtin_fmaf (c5, tt, zmath_pi_asinpih_p0[1]);    // degree 3
    y = __builtin_fmaf (c3, tt * t, t);
  }
  else { // 0.25 <= t <= 0.5
    if (au == 0x3eb7a000u) // |x|=0x1.6f4p-2
      return (au == u) ? 0x1.de400ep-4f : -0x1.de400ep-4f;
    if (au == 0x3ec4c000u) // |x|=0x1.898p-2
      return (au == u) ? 0x1.012006p-3f : -0x1.012006p-3f;
    c5 = __builtin_fmaf (zmath_pi_asinpih_p1[3], tt, zmath_pi_asinpih_p1[2]);
    c1 = __builtin_fmaf (zmath_pi_asinpih_p1[1], tt, zmath_pi_asinpih_p1[0]);
    y = __builtin_fmaf (c5, tt * tt, c1);
    y = t * y;
  }

  if (reduce) // argument reconstruction
  {
    y = ZMATH_PI_ASINPIH_HALF_PI - 2.0f * y;
    if (u >> 31) // x < 0
      y = -y;
  }

  /* now y approximates asin(x), we divide by pi */

  return y * ZMATH_PI_ASINPIH_INV_PI;
}

/* Pinned CORE-MATH source: revision 3bafed75e3c22399531857d63a3fa0e24bfb42fa, src/binary16/acospi/acospif16.c. */
/* Correctly-rounded half-revolution arc-cosine for binary16 value.

Copyright (c) 2025 Paul Zimmermann

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


// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif


typedef union {float f; uint32_t u;} zmath_pi_acospih_b32u32_u;

// the following polynomials were generated using Sollya (cf acos.sollya)

/* degree-4 minimax polynomial for acos(x) over [0,0.25], with relative error
   bounded by 2^-22.943, manually optimized to reduce the number of exceptions
*/
static const float zmath_pi_acospih_p0[] = {0x1.921fb4p0f, -0x1.fffb44p-1f, -0x1.25e6b8p-10f,
                           -0x1.3cc114p-3f, -0x1.a85b22p-5f};

/* degree-4 minimax polynomial for acos(x) over [0.25,0.5], with relative error
   bounded by 2^-20.789, manually optimized to reduce the number of exceptions
*/
static const float zmath_pi_acospih_p1[] = {0x1.91b67p0f, -0x1.f515cap-1f, -0x1.bd043ap-4f,
                           0x1.7e2d5ap-4f, -0x1.190806p-2f};

/* degree-4 minimax polynomial for acos(x)/sqrt(1-x) over [0.5,1],
   with relative error bounded by 2^-23.583 */
static const float zmath_pi_acospih_p2[] = {0x1.91fa1cp0f, -0x1.ae5c5ep-3f, 0x1.31640cp-4f,
                           -0x1.98038p-6f, 0x1.251b5p-8f};

_Float16 zmath_core_acospih (_Float16 x)
{
  zmath_pi_acospih_b32u32_u v = {.f = x};
  uint32_t u = v.u;
  uint32_t au = u & 0x7fffffffu;

  if (au >= 0x3f800000u) { // NaN, Inf, or |x| >= 1
    if (au == 0x3f800000u)
      return (u == 0x3f800000u) ? 0.0f16 : 1.0f;
    if ((au >> 23) == 0x3ff && ((au & 0x7fffff) != 0)) // qNaN or sNaN
      return x;
    return 0.0f / 0.0f; // will signal invalid and return sNaN
  }

  float t = v.f, tt = t * t, c1, c3, y;

  if (u >> 31) // x < 0
    t = -t;

#define ZMATH_PI_ACOSPIH_INV_PI 0x1.45f306p-2f

  if (au < 0x3e800000u) { // |x| < 0.25
    if (au == 0) return 0.5f; // would be inexact with the code below
    if (au == 0x3c55a000u) // |x| = 0x1.ab4p-7
      return (au == u) ? 0x1.fbbffcp-2f : 0x1.022002p-1f;
    if (au == 0x3d8a2000u) // |x| = 0x1.144p-4
      return (au == u) ? 0x1.ea0002p-2f : 0x1.0afffep-1f;
    if (u == 0x3d21c000u) // x = 0x1.438p-5
      return 0x1.f31ffep-2f;
    if (u == 0x3db3a000u) // x = 0x1.674p-4
      return 0x1.e36002p-2f;
    c1 = __builtin_fmaf (zmath_pi_acospih_p0[2], t, zmath_pi_acospih_p0[1]);
    c3 = __builtin_fmaf (zmath_pi_acospih_p0[4], t, zmath_pi_acospih_p0[3]);
    y = __builtin_fmaf (c3, tt, c1);
    y = __builtin_fmaf (y, t, zmath_pi_acospih_p0[0]);
    y = y * ZMATH_PI_ACOSPIH_INV_PI;
    if (u >> 31) // x < 0
      y = 1.0f - y;
  }
  else if (au < 0x3f000000u) { // 0.25 <= |x| < 0.5
    c1 = __builtin_fmaf (zmath_pi_acospih_p1[2], t, zmath_pi_acospih_p1[1]);
    c3 = __builtin_fmaf (zmath_pi_acospih_p1[4], t, zmath_pi_acospih_p1[3]);
    y = __builtin_fmaf (c3, tt, c1);
    y = __builtin_fmaf (y, t, zmath_pi_acospih_p1[0]);
    y = y * ZMATH_PI_ACOSPIH_INV_PI;
    if (u >> 31) // x < 0
      y = 1.0f - y;
  }
  else { // 0.5 <= |x| <= 1
    if (au == 0x3f3ca000u) // |x| = 0x1.794p-1
      return (au == u) ? 0x1.e4000ap-3f : 0x1.86fffep-1f;
    c1 = __builtin_fmaf (zmath_pi_acospih_p2[2], t, zmath_pi_acospih_p2[1]);
    c3 = __builtin_fmaf (zmath_pi_acospih_p2[4], t, zmath_pi_acospih_p2[3]);
    y = __builtin_fmaf (c3, tt, c1);
    y = __builtin_fmaf (y, t, zmath_pi_acospih_p2[0]);
    y = __builtin_sqrtf (1.0f - t) * y;
    y = y * ZMATH_PI_ACOSPIH_INV_PI;
    if (u >> 31) // x < 0
      y = 1.0f - y;
  }

  return y;
}

/* Pinned CORE-MATH source: revision 3bafed75e3c22399531857d63a3fa0e24bfb42fa, src/binary16/atanpi/atanpif16.c. */
/* Correctly-rounded half-revolution arc-tangent for binary16 value.

Copyright (c) 2025 Paul Zimmermann

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


// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif


typedef union {float f; uint32_t u;} zmath_pi_atanpih_b32u32_u;

// the following polynomials were generated using Sollya (cf atan.sollya)

/* Degree-7 minimax polynomial for atan(x) over [0,0.25], with relative error
   bounded by 2^-25.419, with coefficients of odd degree only, and degree-1
   coefficient forced to 1.
*/
static const float zmath_pi_atanpih_p0[] = {0x1.fffffcp-1, -0x1.55546cp-2, 0x1.98d0ep-3,
                           -0x1.0c7c54p-3};

/* degree-4 minimax polynomial for atan(x) over [0.25,0.5], with relative error
   bounded by 2^-22.573, and coefficients later optimized */
static const float zmath_pi_atanpih_p1[] = {0x1.411612p-14, 0x1.ff076ep-1, 0x1.1ee64cp-6,
                           -0x1.a6fc96p-2, 0x1.81dd3ep-3};

/* degree-4 minimax polynomial for atan(x) over [0.5,0.75], with relative error
   bounded by 2^-21.757, and coefficients later optimized */
static const float zmath_pi_atanpih_p2[] = {-0x1.95964cp-8, 0x1.0bad74p+0, -0x1.e652e2p-4,
                           -0x1.e70ab4p-3, 0x1.a5eb88p-4};

/* degree-4 minimax polynomial for atan(x) over [0.75,1], with relative error
   bounded by 2^-23.027, and coefficients later optimized */
static const float zmath_pi_atanpih_p3[] = {-0x1.f75b26p-6, 0x1.2d3d1cp+0, -0x1.882376p-2,
                           0x1.d18c84p-13, 0x1.6aa26p-6};

_Float16 zmath_core_atanpih (_Float16 x)
{
  zmath_pi_atanpih_b32u32_u v = {.f = x};
  uint32_t u = v.u;
  uint32_t au = u & 0x7fffffffu;

#define ZMATH_PI_ATANPIH_HALF_PI 0x1.921fb6p+0f // approximation of pi/2
#define ZMATH_PI_ATANPIH_INV_PI  0x1.45f306p-2f // approximation of 1/pi

  if (au >= 0x7f800000u) { // NaN or Inf
    if (au == 0x7f800000u) // +/-Inf
      return (u == 0x7f800000u) ? 0.5f : -0.5f;
    return x + x; // will signal invalid for sNaN and return qNaN
  }
  
  float t = v.f;

  // we first approximate atan(x), then divide by pi

  // for x < 0 we use atan(-x) = -atan(x)
  int neg = u >> 31;
  static const float Neg[] = {1.0f, -1.0f};
  float s = Neg[neg];
  t = t * s;

  // now t >= 0

  // for x > 1 we use atan(x) = pi/2 - atan(1/x)
  int reduce = au > 0x3f800000u;
  if (reduce) {
    // deal with exceptional cases
    if (au == 0x3fd24000u) return s * 0x1.4dbffep-2f; // |x| = 0x1.a48p+0
    if (au == 0x42ba4000u) return s * 0x1.fc8002p-2f; // |x| = 0x1.748p+6
    if (au == 0x412ca000u) return s * 0x1.e1dffep-2f;  // |x| = 0x1.594p+3
    if (au == 0x41fb4000u) return s * 0x1.f5a002p-2f; // |x| = 0x1.f68p+4
    t = 1.0f / t;
    v.f = t;
    au = v.u;
  }

  // now 0 <= t <= 1

  const float *p;
  float y, c1, c2, c3, c5, tt = t * t;
  if (au < 0x3e800000u) { // t < 0.25
    if (!reduce && au < 0x3a800000u) { // t < 0x1p-10
      // use direct degree-1 polynomial for atanpi (cf atanpi.sollya)
      // deal with exceptional cases
      if (au == 0x3a052000u) return s * 0x1.52ffep-13f;  // |x| = 0x1.0a4p-11
      if (au == 0x3a318000u) return s * 0x1.c3fffep-13f; // |x| = 0x1.63p-11
      return s * t * __builtin_fmaf (-0x1.cb4bcep-14f, t, 0x1.45f30ap-2f);
    }
    // deal with exceptional cases 
    if (au == 0x3af4c000u) return s * 0x1.37a002p-11f; // |x| = 0x1.e98p-10
    if (au == 0x3be78000u) return s * 0x1.26c004p-9f;  // |x| = 0x1.cfp-8
    if (au == 0x3b92a000u) return s * 0x1.756002p-10f; // |x| = 0x1.254p-8
    if (au == 0x3e1ea000u) return s * 0x1.90c002p-5f;  // |x| = 0x1.3d4p-3
    if (au == 0x3e4b2000u) return s * 0x1.fea002p-5f;  // |x| = 0x1.964p-3
    p = zmath_pi_atanpih_p0;
    c5 = __builtin_fmaf (p[3], tt, p[2]);
    c1 = __builtin_fmaf (p[1], tt, p[0]);
    c1 = __builtin_fmaf (c5, tt * tt, c1);
    y = t * c1;
  }
  else {
    // use zmath_pi_atanpih_p1 for t <= 0.5, zmath_pi_atanpih_p2 for 0.5 < t <= 0.75, zmath_pi_atanpih_p3 for 0.75 < t <= 1
    if (au <= 0x3f000000u)
      p = zmath_pi_atanpih_p1;
    else if (au <= 0x3f400000u)
      p = zmath_pi_atanpih_p2;
    else
    {
      p = zmath_pi_atanpih_p3;
      // deal with exceptional cases
      if (au == 0x3f800000u) return s * 0x1p-2f;        // |x| = 1
      if (au == 0x3f712000u) return s * 0x1.ec7feap-3f; // |x| = 0x1.e24p-1
    }
    c3 = __builtin_fmaf (p[4], t, p[3]);
    c2 = __builtin_fmaf (c3, t, p[2]);
    y = __builtin_fmaf (p[1], t, p[0]);
    y = __builtin_fmaf (c2, tt, y);
  }

  if (reduce) y = ZMATH_PI_ATANPIH_HALF_PI - y; // argument reconstruction

  if (neg) y = -y;

  // divide by pi
  return y * ZMATH_PI_ATANPIH_INV_PI;
}

/* Pinned CORE-MATH source: revision 3bafed75e3c22399531857d63a3fa0e24bfb42fa, src/binary16/atan2pi/atan2pif16.c. */
/* Correctly-rounded half revolution arctangent function of two binary16 values.

Copyright (c) 2022-2025 Alexei Sibidanov and Paul Zimmermann

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

/* This code is based on the binary32 code atan2pif.c:
   at input we convert the inputs to _Float16 (exactly),
   we then use the same code than for atan2pif,
   and at output we round to _Float16.
   The changes with respect to atan2pif.c are marked with a comment
   "specific atan2pif16". */


// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif


typedef union {float f; uint32_t u;} zmath_pi_atan2pih_b32u32_u;
typedef union {double f; uint64_t u;} zmath_pi_atan2pih_b64u64_u;

static inline double zmath_pi_atan2pih_muldd(double xh, double xl, double ch, double cl, double *l){
  double ahlh = ch*xl, alhh = cl*xh, ahhh = ch*xh, ahhl = __builtin_fma(ch, xh, -ahhh);
  ahhl += alhh + ahlh;
  ch = ahhh + ahhl;
  *l = (ahhh - ch) + ahhl;
  return ch;
}

static double zmath_pi_atan2pih_polydd(double xh, double xl, int n, const double c[][2], double *l){
  int i = n-1;
  double ch = c[i][0], cl = c[i][1];
  while(--i>=0){
    ch = zmath_pi_atan2pih_muldd(xh,xl,ch,cl,&cl);
    double th = ch + c[i][0], tl = (c[i][0] - th) + ch;
    ch = th;
    cl += tl + c[i][1];
  }
  *l = cl;
  return ch;
}

// specific atan2pif16: inputs name changed from y,x to yf16,xf16
_Float16 zmath_core_atan2pih(_Float16 yf16, _Float16 xf16){
  float y = yf16, x = xf16; // specific atan2pif16
  static const double cn[] =
    {0x1.45f306dc9c883p-2, 0x1.988d83a142adap-1, 0x1.747bebf492057p-1, 0x1.2cc5645094ff3p-2,
     0x1.a0521c711ab66p-5, 0x1.881b8058b9a0dp-9, 0x1.b16ff514a0afp-16};
  static const double cd[] =
    {0x1p+0, 0x1.6b8b143a3f6dap+1, 0x1.8421201d18ed5p+1, 0x1.8221d086914ebp+0, 0x1.670657e3a07bap-2,
     0x1.0f4951fd1e72dp-5, 0x1.b3874b8798286p-11};
  static const double m[] = {0, 1};
  static const double off[] = {0.0f, 0.5f, 1.0f, 0.5f, -0.0f, -0.5f, -1.0f, -0.5f};
  static const float sgnf[] = {1,-1};
  static const double sgn[] = {1,-1};
  zmath_pi_atan2pih_b32u32_u tx = {.f = x}, ty = {.f = y};
  uint32_t ux = tx.u, uy = ty.u, ax = ux&(~0u>>1), ay = uy&(~0u>>1);
  if(__builtin_expect(ay >= (0xff<<23)||ax >= (0xff<<23), 0)){
    if(ay > (0xff<<23)) return x + y; // nan
    if(ax > (0xff<<23)) return x + y; // nan
    uint32_t yinf = ay==(0xff<<23), xinf = ax==(0xff<<23);
    if(yinf&xinf){
      if(ux>>31)
	return 0.75f*sgnf[uy>>31];
      else
	return 0.25f*sgnf[uy>>31];
    }
    if(xinf){
      if(ux>>31)
	return sgnf[uy>>31];
      else
	return 0.0f*sgnf[uy>>31];
    }
    if(yinf){
      return 0.5f*sgnf[uy>>31];
    }
  }
  if(__builtin_expect(ay==0, 0)){
    if(__builtin_expect(!(ay|ax),0)){
      uint32_t i = (uy>>31)*4 + (ux>>31)*2;
      return off[i];
    }
    if(!(ux>>31)) return 0.0f*sgnf[uy>>31];
  }
  if(__builtin_expect(ax==ay, 0)){
    static const float s[] = {0.25,0.75,-0.25,-0.75};
    uint32_t i = (uy>>31)*2 + (ux>>31);
    return s[i];
  }
  uint32_t gt = ay>ax, i = (uy>>31)*4 + (ux>>31)*2 + gt;
  
  double zx = x, zy = y;
  double z = (m[gt]*zx + m[1-gt]*zy)/(m[gt]*zy + m[1-gt]*zx);
  double r = cn[0], z2 = z*z;
  z *= sgn[gt];
  // avoid spurious underflow in the polynomial evaluation excluding extremely small arguments
  if(__builtin_expect(z2>0x1p-54, 1)){
    double z4 = z2*z2, z8 = z4*z4;
    double cn0 =     r + z2*cn[1];
    double cn2 = cn[2] + z2*cn[3];
    double cn4 = cn[4] + z2*cn[5];
    double cn6 = cn[6];
    cn0 += z4*cn2;
    cn4 += z4*cn6;
    cn0 += z8*cn4;
    double cd0 = cd[0] + z2*cd[1];
    double cd2 = cd[2] + z2*cd[3];
    double cd4 = cd[4] + z2*cd[5];
    double cd6 = cd[6];
    cd0 += z4*cd2;
    cd4 += z4*cd6;
    cd0 += z8*cd4;
    r = cn0/cd0;
  }
  r = z*r + off[i];
  zmath_pi_atan2pih_b64u64_u res = {.f = r};
  if(__builtin_expect((res.u<<1) > 0x6d40000000000000 && ((res.u + 8)&0xfffffff) <= 16, 0)){
    // |res| > 0x1p-149
    if(ax==ay) {
      static const double off2[] = {0.25, 0.75, -0.25, -0.75};
      r = off2[(uy>>31)*2 + (ux>>31)];
    } else {
      double zh,zl;
      if(!gt){
	zh = zy/zx;
	zl = __builtin_fma(zh,-zx,zy)/zx;
      } else {
	zh = zx/zy;
	zl = __builtin_fma(zh,-zy,zx)/zy;
      }
      double z2l, z2h = zmath_pi_atan2pih_muldd(zh,zl,zh,zl,&z2l);
      static const double c[][2] =
	{{0x1.45f306dc9c883p-2, -0x1.6b01ec5513324p-56}, {-0x1.b2995e7b7b604p-4, 0x1.e402b0c13eedcp-58},
	 {0x1.04c26be3b06cfp-4, -0x1.571d178a53efp-60}, {-0x1.7483758e69c03p-5, 0x1.819a6ed7aaf38p-63},
	 {0x1.21bb9452523ffp-5, -0x1.234d866fb9807p-60}, {-0x1.da1bace3cc54ep-6, -0x1.c84f6ada49294p-64},
	 {0x1.912b1c23345ddp-6, -0x1.534890fbc165p-60}, {-0x1.5bade52f5f52ap-6, 0x1.f783bafc832f6p-60},
	 {0x1.32c69d084c5cp-6, 0x1.042d155953025p-60}, {-0x1.127bcfb3e8c7dp-6, -0x1.85aae199a7b6bp-60},
	 {0x1.f0af43b11a731p-7, 0x1.8f0356356663p-61}, {-0x1.c57e86801029ep-7, 0x1.dcdf3e3b38eb4p-61},
	 {0x1.a136408617ea1p-7, 0x1.a71affb36c6c4p-63}, {-0x1.824ac7814ba37p-7, 0x1.8928b295c0898p-61},
	 {0x1.6794e32ea5471p-7, 0x1.0b4334fb41e63p-61}, {-0x1.501d57f643d97p-7, 0x1.516785bf1376ep-61},
	 {0x1.3adf02ff2400ap-7, -0x1.b0e30bb8c8076p-62}, {-0x1.267702f94faap-7, -0x1.7a4d3a1850cc6p-62},
	 {0x1.10dce97099686p-7, 0x1.fcc208eee2571p-61}, {-0x1.eee49cdad8002p-8, -0x1.9109b3f1bab82p-64},
	 {0x1.af93bc191a929p-8, 0x1.069fd3b47d7bp-62}, {-0x1.6240751b54675p-8, -0x1.72dc8cfd03b6fp-62},
	 {0x1.0b61e84080884p-8, 0x1.825824c80941bp-63}, {-0x1.6a72a8a74e3a5p-9, 0x1.8786a82fd117ep-63},
	 {0x1.aede3217d939dp-10, -0x1.93b626982e1fep-68}, {-0x1.b66568f09ebeep-11, -0x1.704a39121d0a5p-66},
	 {0x1.73af3977fa973p-12, -0x1.aa050e2244ea3p-68}, {-0x1.fc69d85ed28c9p-14, 0x1.867f17b764cap-68},
	 {0x1.0c883a9270162p-15, -0x1.6842833896dd9p-70}, {-0x1.9a0b27b6dfe15p-18, 0x1.427fc2f4e1327p-73},
	 {0x1.91e15e7ab5bdcp-21, -0x1.730dbc6279d0dp-77}, {-0x1.7b1119c1ff867p-25, 0x1.145f9980759c4p-79}};
      double pl, ph = zmath_pi_atan2pih_polydd(z2h, z2l, 32, c, &pl);
      zh *= sgn[gt];
      zl *= sgn[gt];
      ph = zmath_pi_atan2pih_muldd(zh,zl,ph,pl,&pl);
      double sh = ph + off[i], sl = ((off[i] - sh) + ph) + pl;
      float rf = sh;
      double th = rf, dh = sh - th, tm = dh + sl;
      r = th + tm;
      zmath_pi_atan2pih_b64u64_u d = {.f = r - th};
      if( !(d.u<<12) ) {
	double ad = __builtin_fabs(d.f), am = __builtin_fabs(tm);
	if(ad > am) r -= d.f*0x1p-10;
	if(ad < am) r += d.f*0x1p-10;
      }
    }
  }
  _Float16 rf = r; // specific atan2pif16
  return rf;
}
