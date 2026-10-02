/* Direct correctly-rounded inverse-pi kernels for binary32. */
#include <stdint.h>

#include "coremath_fma.h"

/* Pinned CORE-MATH source: revision 3bafed75e3c22399531857d63a3fa0e24bfb42fa, src/binary32/asinpi/asinpif.c. */
/* Correctly-rounded half-revolution arc-sine function for binary32 value.

Copyright (c) 2022 Alexei Sibidanov.

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
typedef union {float f; uint32_t u;} zmath_pi_asinpif_b32u32_u;

float zmath_core_asinpif(float x){
  float ax = __builtin_fabsf(x);
  double az = ax, z = x;
  zmath_pi_asinpif_b32u32_u t = {.f = x};
  int32_t e = (t.u>>23)&0xff;
  if(__builtin_expect(e>=127, 0)){ // |x| >= 1 or nan
    if(ax == 1.0f) return __builtin_copysignf(0.5f, x); // |x| = 1
    if(e==0xff && (t.u<<9)) return x+x; // x = nan
return __builtin_nanf("1"); // |x| > 1
  }
  int32_t s = 146 - e, i = 0;
  // s<32 corresponds to |x| >= 2^-12
  if(__builtin_expect(s<32, 1)) i = ((t.u&(~0u>>9))|1<<23)>>s;
  static const double ch[][8] = {
    {0x1.45f306dc9c882p-2, 0x1.b2995e7b7dc2fp-5, 0x1.8723a1cf50c7ep-6, 0x1.d1a4591d16a29p-7,
    0x1.3ce3aa68ddaeep-7, 0x1.d3182ab0cc1bfp-8, 0x1.62b379a8b88e3p-8, 0x1.6811411fcfec2p-8},
    {0x1.ffffffffd3cdap-2, -0x1.17cc1b3355fddp-4, 0x1.d067a1e8d5a99p-6, -0x1.08e16fb09314ap-6,
    0x1.5eed43d42dcb2p-7, -0x1.f58baca7acc71p-8, 0x1.5dab64e2dcf15p-8, -0x1.59270e30797acp-9},
    {0x1.fffffff7c4617p-2, -0x1.17cc149ded3a2p-4, 0x1.d0654d4cb2c1ap-6, -0x1.08c3ba713d33ap-6,
    0x1.5d2053481079cp-7, -0x1.e485ebc545e7ep-8, 0x1.303baca167dddp-8, -0x1.dee8d16d06b38p-10},
    {0x1.ffffffa749848p-2, -0x1.17cbe7155935p-4, 0x1.d05a312269adfp-6, -0x1.0862b3ee617d7p-6,
    0x1.5920708db2a73p-7, -0x1.cb0463b3862c3p-8, 0x1.02b82478f95d7p-8, -0x1.52a7b8579e729p-10},
    {0x1.fffffe1f92bb5p-2, -0x1.17cb3e74c64e3p-4, 0x1.d03af67311cbfp-6, -0x1.079441cbfc7ap-6,
    0x1.52b4287805a61p-7, -0x1.ac3286d604a98p-8, 0x1.b2f1210d9701bp-9, -0x1.e740ddc25afd6p-11},
    {0x1.fffff92beb6e2p-2, -0x1.17c986fe9518bp-4, 0x1.cff98167c9a5ep-6, -0x1.0638b591eae52p-6,
    0x1.4a0803828959ep-7, -0x1.8adeca229f11dp-8, 0x1.6b9a7ba05dfcep-9, -0x1.640521a43b2dp-11},
    {0x1.ffffeccee5bfcp-2, -0x1.17c5f1753f5eap-4, 0x1.cf874e4fe258fp-6, -0x1.043e6cf77b256p-6,
    0x1.3f7db42227d92p-7, -0x1.691a6fa2a2882p-8, 0x1.2f6543162bc61p-9, -0x1.07d5da05822b6p-11},
    {0x1.ffffd2f64431dp-2, -0x1.17bf8208c10c1p-4, 0x1.ced7487cdb124p-6, -0x1.01a0d30932905p-6,
    0x1.3388f99b254dap-7, -0x1.4844e245c65bdp-8, 0x1.fa777150197c6p-10, -0x1.8c1ecf16a05c8p-12},
    {0x1.ffffa36d1712ep-2, -0x1.17b523971bd4ep-4, 0x1.cddee26de2deep-6, -0x1.fccb00abaaabcp-7,
    0x1.269afc3622342p-7, -0x1.2933152686752p-8, 0x1.a76d4956cc9a3p-10, -0x1.2ce7d6dc651cep-12},
    {0x1.ffff5402ab3a1p-2, -0x1.17a5ba85da77ap-4, 0x1.cc96894e05c02p-6, -0x1.f532143cb832ep-7,
    0x1.19180b660ff09p-7, -0x1.0c57417a78b3cp-8, 0x1.62e26cbd7bb1ep-10, -0x1.ce28d33fe1df3p-13},
    {0x1.fffed8d639751p-2, -0x1.1790349f3ae76p-4, 0x1.caf9a4fd1b398p-6, -0x1.ec986b111342ep-7,
    0x1.0b53c3ad4baa4p-7, -0x1.e3c2282eeace4p-9, 0x1.2a55369f55bbep-10, -0x1.667fe48c396e8p-13},
    {0x1.fffe24b714161p-2, -0x1.177394fbcb719p-4, 0x1.c90652d920ebdp-6, -0x1.e3239197bddf1p-7,
    0x1.fb2188525b025p-8, -0x1.b3aadd451afc7p-9, 0x1.f74020f31fdabp-11, -0x1.18b0cb246768dp-13},
    {0x1.fffd298bec9e2p-2, -0x1.174efbfd34648p-4, 0x1.c6bcfe48ea92bp-6, -0x1.d8f9f2a16157cp-7,
    0x1.e0044f56c8864p-8, -0x1.883e2347fe76cp-9, 0x1.a9f0e3c1b7af5p-11, -0x1.bb5acc0e60825p-14},
    {0x1.fffbd8b784c4dp-2, -0x1.1721abdd3722ep-4, 0x1.c41fee756d4bp-6, -0x1.ce40bccf8065fp-7,
    0x1.c59b684b70ef9p-8, -0x1.6133d027996b3p-9, 0x1.69cad01106397p-11, -0x1.60f8e45494156p-14},
    {0x1.fffa23749cf88p-2, -0x1.16eb0a8285c06p-4, 0x1.c132d762e1b0dp-6, -0x1.c31a959398f4ep-7,
    0x1.ac1c5b46bc8ap-8, -0x1.3e34f1abe51dcp-9, 0x1.346738737c0b9p-11, -0x1.1b227a3f5c75p-14},
    {0x1.fff7fb25bb407p-2, -0x1.16aaa14d7564p-4, 0x1.bdfa75fca5ff2p-6, -0x1.b7a6e260d079cp-7,
    0x1.93ab06911033cp-8, -0x1.1ee5560967fd5p-9, 0x1.07d31060838bfp-11, -0x1.c96f33a283115p-15}};
  const double *c = ch[i];
  double z2 = z*z, z4 = z2*z2;
  if(__builtin_expect(!i, 0)){ // |x| < 2^-4
    double c0 = c[0] + z2*c[1];
    double c2 = c[2] + z2*c[3];
    double c4 = c[4] + z2*c[5];
    double c6 = c[6] + z2*c[7];
    c0 += c2*z4;
    c4 += c6*z4;
    c0 += c4*(z4*z4);
    return z*c0;
  } else { // |x| >= 2^-4
    double f = __builtin_sqrt(1-az);
    double c0 = c[0] + az*c[1];
    double c2 = c[2] + az*c[3];
    double c4 = c[4] + az*c[5];
    double c6 = c[6] + az*c[7];
    c0 += c2*z2;
    c4 += c6*z2;
    c0 += c4*z4;
    double r = zmath_fma(-c0, __builtin_copysign(f, x), __builtin_copysign(0.5, x));
    return r;
  }
}

/* Pinned CORE-MATH source: revision 3bafed75e3c22399531857d63a3fa0e24bfb42fa, src/binary32/acospi/acospif.c. */
/* Correctly-rounded half-revolution arc-cosine function for binary32 value.

Copyright (c) 2022 Alexei Sibidanov.

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
typedef union {float f; uint32_t u;} zmath_pi_acospif_b32u32_u;

float zmath_core_acospif(float x){
  float ax = __builtin_fabsf(x);
  double az = ax, z = x;
  zmath_pi_acospif_b32u32_u t = {.f = x};
  int e = (t.u>>23)&0xff;
  if(__builtin_expect(e>=127, 0)){
    if(x == 1.0f) return 0.0f;
    if(x ==-1.0f) return 1.0f;
    if(e==0xff && (t.u<<9)) return x + x; // nan
return __builtin_nanf("1");
  }
  int s = 146 - e, i = 0;
  if(__builtin_expect(s<32, 1)) i = ((t.u&(~0u>>9))|1<<23)>>s;
  static const double ch[][8] = {
    {0x1.45f306dc9c882p-2, 0x1.b2995e7b7dc2fp-5, 0x1.8723a1cf50c7ep-6, 0x1.d1a4591d16a29p-7,
     0x1.3ce3aa68ddaeep-7, 0x1.d3182ab0cc1bfp-8, 0x1.62b379a8b88e3p-8, 0x1.6811411fcfec2p-8},
    {0x1.ffffffffd3cd9p-2, -0x1.17cc1b3355fd5p-4, 0x1.d067a1e8d5a99p-6, -0x1.08e16fb09314ap-6,
    0x1.5eed43d42dcb2p-7, -0x1.f58baca7acc71p-8, 0x1.5dab64e2dcf15p-8, -0x1.59270e30797acp-9},
    {0x1.fffffff7c4622p-2, -0x1.17cc149ded3a2p-4, 0x1.d0654d4cb2c1ap-6, -0x1.08c3ba713d33ap-6,
    0x1.5d2053481079cp-7, -0x1.e485ebc545e7ep-8, 0x1.303baca167dddp-8, -0x1.dee8d16d06b38p-10},
    {0x1.ffffffa749848p-2, -0x1.17cbe7155935p-4, 0x1.d05a312269adfp-6, -0x1.0862b3ee617d7p-6,
    0x1.5920708db2a73p-7, -0x1.cb0463b3862c3p-8, 0x1.02b82478f95d7p-8, -0x1.52a7b8579e729p-10},
    {0x1.fffffe1f92bb5p-2, -0x1.17cb3e74c64e3p-4, 0x1.d03af67311cbfp-6, -0x1.079441cbfc7ap-6,
    0x1.52b4287805a61p-7, -0x1.ac3286d604a98p-8, 0x1.b2f1210d9701bp-9, -0x1.e740ddc25afd6p-11},
    {0x1.fffff92beb6e2p-2, -0x1.17c986fe9518bp-4, 0x1.cff98167c9a5ep-6, -0x1.0638b591eae52p-6,
    0x1.4a0803828959ep-7, -0x1.8adeca229f11dp-8, 0x1.6b9a7ba05dfcep-9, -0x1.640521a43b2dp-11},
    {0x1.ffffeccee5bfcp-2, -0x1.17c5f1753f5eap-4, 0x1.cf874e4fe258fp-6, -0x1.043e6cf77b256p-6,
    0x1.3f7db42227d92p-7, -0x1.691a6fa2a2882p-8, 0x1.2f6543162bc61p-9, -0x1.07d5da05822b6p-11},
    {0x1.ffffd2f64431dp-2, -0x1.17bf8208c10c1p-4, 0x1.ced7487cdb124p-6, -0x1.01a0d30932905p-6,
    0x1.3388f99b254dap-7, -0x1.4844e245c65bdp-8, 0x1.fa777150197c6p-10, -0x1.8c1ecf16a05c8p-12},
    {0x1.ffffa36d1712ep-2, -0x1.17b523971bd4ep-4, 0x1.cddee26de2deep-6, -0x1.fccb00abaaabcp-7,
    0x1.269afc3622342p-7, -0x1.2933152686752p-8, 0x1.a76d4956cc9a3p-10, -0x1.2ce7d6dc651cep-12},
    {0x1.ffff5402ab3a1p-2, -0x1.17a5ba85da77ap-4, 0x1.cc96894e05c02p-6, -0x1.f532143cb832ep-7,
    0x1.19180b660ff09p-7, -0x1.0c57417a78b3cp-8, 0x1.62e26cbd7bb1ep-10, -0x1.ce28d33fe1df3p-13},
    {0x1.fffed8d639751p-2, -0x1.1790349f3ae76p-4, 0x1.caf9a4fd1b398p-6, -0x1.ec986b111342ep-7,
    0x1.0b53c3ad4baa4p-7, -0x1.e3c2282eeace4p-9, 0x1.2a55369f55bbep-10, -0x1.667fe48c396e8p-13},
    {0x1.fffe24b714161p-2, -0x1.177394fbcb719p-4, 0x1.c90652d920ebdp-6, -0x1.e3239197bddf1p-7,
    0x1.fb2188525b025p-8, -0x1.b3aadd451afc7p-9, 0x1.f74020f31fdabp-11, -0x1.18b0cb246768dp-13},
    {0x1.fffd298bec9e2p-2, -0x1.174efbfd34648p-4, 0x1.c6bcfe48ea92bp-6, -0x1.d8f9f2a16157cp-7,
    0x1.e0044f56c8864p-8, -0x1.883e2347fe76cp-9, 0x1.a9f0e3c1b7af5p-11, -0x1.bb5acc0e60825p-14},
    {0x1.fffbd8b784c4dp-2, -0x1.1721abdd3722ep-4, 0x1.c41fee756d4bp-6, -0x1.ce40bccf8065fp-7,
    0x1.c59b684b70ef9p-8, -0x1.6133d027996b3p-9, 0x1.69cad01106397p-11, -0x1.60f8e45494156p-14},
    {0x1.fffa23749cf88p-2, -0x1.16eb0a8285c06p-4, 0x1.c132d762e1b0dp-6, -0x1.c31a959398f4ep-7,
    0x1.ac1c5b46bc8ap-8, -0x1.3e34f1abe51dcp-9, 0x1.346738737c0b9p-11, -0x1.1b227a3f5c75p-14},
    {0x1.fff7fb25bb407p-2, -0x1.16aaa14d7564p-4, 0x1.bdfa75fca5ff2p-6, -0x1.b7a6e260d079cp-7,
    0x1.93ab06911033cp-8, -0x1.1ee5560967fd5p-9, 0x1.07d31060838bfp-11, -0x1.c96f33a283115p-15},
  };
  const double *c = ch[i];
  double z2 = z*z, z4 = z2*z2;
  if(__builtin_expect(!i, 0)){
    double c0 = c[0] + z2*c[1];
    double c2 = c[2] + z2*c[3];
    double c4 = c[4] + z2*c[5];
    double c6 = c[6] + z2*c[7];
    c0 += c2*z4;
    c4 += c6*z4;
    /* For |x| <= 0x1.0fd288p-127, c0 += c4*(z4*z4) would raise a spurious
       underflow exception, we use an FMA instead, where c4 * z4 does not
       underflow. */
#ifndef __FP_FAST_FMA
    if (__builtin_expect (ax > 0x1.0fd288p-127f, 0))
      c0 = (c4 * z4) * z4 + c0;
    else
#endif
      c0 = zmath_fma (c4 * z4, z4, c0);
    return 0.5 - z*c0;
  } else {
    double f = __builtin_sqrt(1-az);
    double c0 = c[0] + az*c[1];
    double c2 = c[2] + az*c[3];
    double c4 = c[4] + az*c[5];
    double c6 = c[6] + az*c[7];
    c0 += c2*z2;
    c4 += c6*z2;
    c0 += c4*z4;
    static const double o[] = {0,1};
    double r = o[t.u>>31] + c0*__builtin_copysign(f, x);
    return r;
  }
}

/* Pinned CORE-MATH source: revision 3bafed75e3c22399531857d63a3fa0e24bfb42fa, src/binary32/atanpi/atanpif.c. */
/* Correctly-rounded half-revolution arc-tangent of binary32 value.

Copyright (c) 2022 Alexei Sibidanov.

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
typedef union {float f; uint32_t u;} zmath_pi_atanpif_b32u32_u;

float zmath_core_atanpif(float x){
  zmath_pi_atanpif_b32u32_u t = {.f = x};
  int32_t e = (t.u>>23)&0xff, gt = e>=127;
  if(__builtin_expect(e>127+24, 0)) { // |x| >= 2^25
    float f = __builtin_copysignf(0.5f, x);
    if(__builtin_expect(e==0xff, 0)) {
      if(t.u<<9) return x + x; // nan
      return f; // inf
    }
    // Warning: 0x1.45f306p-2f / x underflows for |x| >= 0x1.45f306p+124
    if (__builtin_fabsf (x) >= 0x1.45f306p+124f)
      return f - __builtin_copysignf(0x1p-26f, x);
    else
      return f - 0x1.45f306p-2f / x;
  }
  double z = x;
  if (__builtin_expect(e<127-13, 0)){ // |x| < 2^-13
    double sx = z*0x1.45f306dc9c883p-2;
    if (__builtin_expect(e<127-25, 0)) { // |x| < 2^-25
      return sx;
    }
    return sx - (0x1.5555555555555p-2*sx)*(z*z);
  }
  uint32_t ax = t.u&(~0u>>1);
  if(__builtin_expect(ax == 0x3fa267ddu, 0)) return __builtin_copysignf(0x1.267004p-2f,x) - __builtin_copysignf(0x1p-55f,x);
  if(__builtin_expect(ax == 0x3f693531u, 0)) return __builtin_copysignf(0x1.e1a662p-3f,x) + __builtin_copysignf(0x1p-28f,x);
  if(__builtin_expect(ax == 0x3f800000u, 0)) return __builtin_copysignf(0x1p-2f,x);
  if(gt) z = 1/z;
  double z2 = z*z, z4 = z2*z2, z8 = z4*z4;
  static const double cn[] =
    {0x1.45f306dc9c882p-2, 0x1.733b561bc23d5p-1, 0x1.28d9805bdfbf2p-1,
     0x1.8c3ba966ae287p-3, 0x1.94a7f81ee634bp-6, 0x1.a6bbf6127a6dfp-11};
  static const double cd[] =
    {0x1p+0, 0x1.4e3b3ecc2518fp+1, 0x1.3ef4a360ff063p+1, 0x1.0f1dc55bad551p+0,
     0x1.8da0fecc018a4p-3, 0x1.8fa87803776bfp-7, 0x1.dadf2ca0acb43p-14};
  double cn0 = cn[0] + z2*cn[1];
  double cn2 = cn[2] + z2*cn[3];
  double cn4 = cn[4] + z2*cn[5];
  cn0 += z4*cn2;
  cn0 += z8*cn4;
  cn0 *= z;
  double cd0 = cd[0] + z2*cd[1];
  double cd2 = cd[2] + z2*cd[3];
  double cd4 = cd[4] + z2*cd[5];
  double cd6 = cd[6];
  cd0 += z4*cd2;
  cd4 += z4*cd6;
  cd0 += z8*cd4;
  double r = cn0/cd0;
  if (gt) r = __builtin_copysign(0.5, z) - r;
  return r;
}

/* Pinned CORE-MATH source: revision 3bafed75e3c22399531857d63a3fa0e24bfb42fa, src/binary32/atan2pi/atan2pif.c. */
/* Correctly-rounded half revolution arctangent function of two binary32 values.

Copyright (c) 2022-2025 Alexei Sibidanov.

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
typedef union {float f; uint32_t u;} zmath_pi_atan2pif_b32u32_u;
typedef union {double f; uint64_t u;} zmath_pi_atan2pif_b64u64_u;

static inline double zmath_pi_atan2pif_muldd(double xh, double xl, double ch, double cl, double *l){
  double ahlh = ch*xl, alhh = cl*xh, ahhh = ch*xh, ahhl = zmath_fma(ch, xh, -ahhh);
  ahhl += alhh + ahlh;
  ch = ahhh + ahhl;
  *l = (ahhh - ch) + ahhl;
  return ch;
}

static double zmath_pi_atan2pif_polydd(double xh, double xl, int n, const double c[][2], double *l){
  int i = n-1;
  double ch = c[i][0], cl = c[i][1];
  while(--i>=0){
    ch = zmath_pi_atan2pif_muldd(xh,xl,ch,cl,&cl);
    double th = ch + c[i][0], tl = (c[i][0] - th) + ch;
    ch = th;
    cl += tl + c[i][1];
  }
  *l = cl;
  return ch;
}

float zmath_core_atan2pif(float y, float x){
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
  zmath_pi_atan2pif_b32u32_u tx = {.f = x}, ty = {.f = y};
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
  } else if (gt || (ux >> 31)) {
    // With z2 <= 2^-54, atan(|z|)/pi is below the f32 midpoint
    // distance from 0.5 or 1; positive-x tiny angles remain nonzero.
    return off[i];
  }
  r = z*r + off[i];
  zmath_pi_atan2pif_b64u64_u res = {.f = r};
  if(__builtin_expect((res.u<<1) > 0x6d40000000000000 && ((res.u + 8)&0xfffffff) <= 16, 0)){
    // |res| > 0x1p-149
    if(ax==ay) {
      static const double off2[] = {0.25, 0.75, -0.25, -0.75};
      r = off2[(uy>>31)*2 + (ux>>31)];
    } else {
      double zh,zl;
      if(!gt){
	zh = zy/zx;
	zl = zmath_fma(zh,-zx,zy)/zx;
      } else {
	zh = zx/zy;
	zl = zmath_fma(zh,-zy,zx)/zy;
      }
      double z2l, z2h = zmath_pi_atan2pif_muldd(zh,zl,zh,zl,&z2l);
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
      double pl, ph = zmath_pi_atan2pif_polydd(z2h, z2l, 32, c, &pl);
      zh *= sgn[gt];
      zl *= sgn[gt];
      ph = zmath_pi_atan2pif_muldd(zh,zl,ph,pl,&pl);
      double sh = ph + off[i], sl = ((off[i] - sh) + ph) + pl;
      float rf = sh;
      double th = rf, dh = sh - th, tm = dh + sl;
      r = th + tm;
      zmath_pi_atan2pif_b64u64_u d = {.f = r - th};
      if( !(d.u<<12) ) {
	double ad = __builtin_fabs(d.f), am = __builtin_fabs(tm);
	if(ad > am) r -= d.f*0x1p-10;
	if(ad < am) r += d.f*0x1p-10;
      }
    }
  }
  float rf = r;
  return rf;
}
