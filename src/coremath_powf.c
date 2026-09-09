/* Correctly-rounded power function for binary32 values.

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

#include <stdint.h>

typedef union {double f; uint64_t u;} b64u64_u;

/* x_bits must encode a positive, finite, nonzero binary32 value. */
static inline uint32_t normalizef_fraction(uint32_t x_bits, int *exponent) {
  uint32_t exponent_field = x_bits >> 23;
  uint32_t fraction = x_bits & 0x7fffff;
  if (__builtin_expect(exponent_field != 0, 1)) {
    *exponent = (int)exponent_field - 127;
  } else {
    int shift = __builtin_clz(fraction) - 8;
    fraction = (fraction << shift) & 0x7fffff;
    *exponent = -126 - shift;
  }
  return fraction;
}

/* __builtin_roundeven was introduced in gcc 10:
   https://gcc.gnu.org/gcc-10/changes.html,
   and in clang 17 */
#if ((defined(__GNUC__) && __GNUC__ >= 10) || (defined(__clang__) && __clang_major__ >= 17)) && !defined(_MSC_VER) && (defined(__aarch64__) || defined(__x86_64__) || defined(__i386__))
# define roundeven_finite(x) __builtin_roundeven (x)
#else
/* round x to nearest integer, breaking ties to even */
static double
roundeven_finite (double x)
{
  double ix;
# if (defined(__GNUC__) || defined(__clang__)) && (defined(__AVX__) || defined(__SSE4_1__) || (__ARM_ARCH >= 8))
#  if defined __AVX__
   __asm__("vroundsd $0x8,%1,%1,%0":"=x"(ix):"x"(x));
#  elif __ARM_ARCH >= 8
   __asm__ ("frintn %d0, %d1":"=w"(ix):"w"(x));
#  else /* __SSE4_1__ */
   __asm__("roundsd $0x8,%1,%0":"=x"(ix):"x"(x));
#  endif
# else
  ix = __builtin_round (x); /* nearest, away from 0 */
  if (__builtin_fabs (ix - x) == 0.5)
  {
    /* if ix is odd, we should return ix-1 if x>0, and ix+1 if x<0 */
    union { double f; uint64_t n; } u, v;
    u.f = ix;
    v.f = ix - __builtin_copysign (1.0, x);
    if (__builtin_ctz (v.n) > __builtin_ctz (u.n))
      ix = v.f;
  }
# endif
  return ix;
}
#endif


static inline double muldd(double xh, double xl, double ch, double cl, double *l){
  double ahlh = ch*xl, alhh = cl*xh, ahhh = ch*xh, ahhl = __builtin_fma(ch, xh, -ahhh);
  ahhl += alhh + ahlh;
  ch = ahhh + ahhl;
  *l = (ahhh - ch) + ahhl;
  return ch;
}

static inline double mulddd(double xh, double xl, double ch, double *l){
  double ahlh = ch*xl, ahhh = ch*xh, ahhl = __builtin_fma(ch, xh, -ahhh);
  ahhl += ahlh;
  ch = ahhh + ahhl;
  *l = (ahhh - ch) + ahhl;
  return ch;
}

static __attribute__((noinline)) double polydd(double xh, double xl, int n, const double c[][2], double *l){
  int i = n-1;
  double ch = c[i][0], cl = c[i][1];
  while(--i>=0){
    ch = muldd(xh,xl,ch,cl,&cl);
    double th = ch + c[i][0], tl = (c[i][0] - th) + ch;
    ch = th;
    cl += tl + c[i][1];
  }
  *l = cl;
  return ch;
}

static __attribute__((noinline)) float as_powf_accurate2(uint64_t, int, float);


float cr_powfFiniteBits(uint32_t x_bits, float y0){
  static const double ix[] = {
    0x1p+0, 0x1.f07c1f07cp-1, 0x1.e1e1e1e1ep-1, 0x1.d41d41d42p-1,
    0x1.c71c71c72p-1, 0x1.bacf914c2p-1, 0x1.af286bca2p-1, 0x1.a41a41a42p-1,
    0x1.99999999ap-1, 0x1.8f9c18f9cp-1, 0x1.861861862p-1, 0x1.7d05f417dp-1,
    0x1.745d1745dp-1, 0x1.6c16c16c1p-1, 0x1.642c8590bp-1, 0x1.5c9882b93p-1,
    0x1.555555555p-1, 0x1.4e5e0a72fp-1, 0x1.47ae147aep-1, 0x1.414141414p-1,
    0x1.3b13b13b1p-1, 0x1.3521cfb2bp-1, 0x1.2f684bda1p-1, 0x1.29e4129e4p-1,
    0x1.249249249p-1, 0x1.1f7047dc1p-1, 0x1.1a7b9611ap-1, 0x1.15b1e5f75p-1,
    0x1.111111111p-1, 0x1.0c9714fbdp-1, 0x1.084210842p-1, 0x1.041041041p-1, 0x1p-1
  };
  
  /* Scale both logarithm table halves and the polynomial by exactly 16
     to remove the dependent run-time y *= 16 multiplication. */
  static const double lix16[][2] = {
    {0x0p+0, 0x0p+0}, {-0x1.6cp-1, 0x1.4b229b87f3f89p-11},
    {-0x1.66p+0, -0x1.fb7d654235799p-11}, {-0x1.08p+1, -0x1.8b119b2c9c87bp-8},
    {-0x1.5cp+1, -0x1.a39fa6533294dp-15}, {-0x1.acp+1, -0x1.ebc5b663dd4b8p-8},
    {-0x1.fcp+1, 0x1.f4a37fe0fa46fp-10}, {-0x1.24p+2, -0x1.01eac33103e6bp-8},
    {-0x1.4ap+2, 0x1.61ed0d15725ep-8}, {-0x1.6ep+2, -0x1.10e6ceb499ba9p-9},
    {-0x1.92p+2, 0x1.115db8ada837dp-8}, {-0x1.b4p+2, -0x1.fafdce266d7aep-8},
    {-0x1.d6p+2, -0x1.d4f80cd19906fp-8}, {-0x1.f8p+2, 0x1.5ea5ccd0a7396p-8},
    {0x1.e8p+2, -0x1.0500d67fe62ebp-9}, {0x1.c8p+2, 0x1.9dc2d41aa4626p-10},
    {0x1.a8p+2, 0x1.ff2e2ff321344p-7}, {0x1.8ap+2, 0x1.130157f4c3a3ep-7},
    {0x1.6cp+2, 0x1.61ed0cad929ccp-7}, {0x1.5p+2, -0x1.2089a632d7949p-7},
    {0x1.32p+2, 0x1.7fdc6dfb2d21ap-7}, {0x1.16p+2, 0x1.380a6c36088f3p-7},
    {0x1.f6p+1, -0x1.3ab7dc7ba81acp-14}, {0x1.cp+1, -0x1.cc2c0061ef1a2p-10},
    {0x1.8ap+1, 0x1.130157c97bbep-8}, {0x1.56p+1, 0x1.ee14ff34c4128p-10},
    {0x1.22p+1, 0x1.b5b854c4fde69p-8}, {0x1.ep+0, 0x1.635d1df7cb0b5p-9},
    {0x1.7ep+0, -0x1.3f6d2636c101ep-9}, {0x1.1cp+0, -0x1.33567f1b193a4p-10},
    {0x1.78p-1, -0x1.8d66c5313a71dp-10}, {0x1.74p-2, 0x1.f7430ee200ep-13},
    {0x0p+0, 0x0p+0}
  };
  double y = y0;
  int e;
  uint32_t fraction = normalizef_fraction(x_bits, &e);
  int j = (fraction + (1u << (23 - 6))) >> (23 - 5), k = j > 13;
  e += k;
  b64u64_u xd = {.u = (0x3ffull << 52) | ((uint64_t)fraction << 29)};
  double z = __builtin_fma(xd.f, ix[j], -1.0);
  double z2 = z*z, z4 = z2*z2;
  double c6 = __builtin_fma(z, -0x1.717d824a520f7p+1, 0x1.a6406efd4b877p+1);
  double c4 = __builtin_fma(z, -0x1.ec709bdf453ecp+1, 0x1.2776c441b72ep+2);
  double c2 = __builtin_fma(z, -0x1.71547652bc4a9p+2, 0x1.ec709dc3a2d0bp+2);
  double c0 = __builtin_fma(z, -0x1.71547652b82fep+3, 0x1.71547652b82fep+4);
  c0 = __builtin_fma(z2, c2, c0);
  c4 = __builtin_fma(z2, c6, c4);
  c0 = __builtin_fma(z4, c4, c0);
  double l = __builtin_fma(z, c0, -lix16[j][1]);
  double zt = (e * 16 - lix16[j][0])*y;
  z = l*y + zt;
  if(__builtin_expect(z>2048, 0)){
#ifdef CORE_MATH_SUPPORT_ERRNO
    errno = ERANGE; // overflow
#endif
    return 0x1p127f*0x1p127f;
  }
  if(__builtin_expect(z<-2400, 0)){
#ifdef CORE_MATH_SUPPORT_ERRNO
    errno = ERANGE; // underflow
#endif
    /* Warning: the expression 0x1p-126f*0x1p-126f should be evaluated
       at run time, with the current rounding mode, and not constant folded
       at compile time. See the comment about -frounding-math in README. */
    return 0x1p-126f*0x1p-126f;
  }
  if(__builtin_fabs(z)<0x1p-26) return 1.0 + z;
  double ia = __builtin_floor(z), h = __builtin_fma(l, y, zt - ia);
  static const uint64_t tb[] = {
    UINT64_C(0x3ff0000000000000), UINT64_C(0x3fefb5586cf9890f), UINT64_C(0x3fef72b83c7d517b), UINT64_C(0x3fef387a6e756238),
    UINT64_C(0x3fef06fe0a31b715), UINT64_C(0x3feedea64c123422), UINT64_C(0x3feebfdad5362a27), UINT64_C(0x3feeab07dd485429),
    UINT64_C(0x3feea09e667f3bcd), UINT64_C(0x3feea11473eb0187), UINT64_C(0x3feeace5422aa0db), UINT64_C(0x3feec49182a3f090),
    UINT64_C(0x3feee89f995ad3ad), UINT64_C(0x3fef199bdd85529c), UINT64_C(0x3fef5818dcfba487), UINT64_C(0x3fefa4afa2a490da)
  };
  /* ia is an integer in [-2400,2048]. The exact bias addition retains
     ia modulo 2^16. Each table entry subtracts j*2^48 from its encoding;
     adding ia*2^48 therefore scales it by 2^floor(ia/16), exactly. */
  b64u64_u integral_bits = {.f = ia + 0x1.8p52};
  b64u64_u scaled = {.u = tb[integral_bits.u & 0xf] + (integral_bits.u << 48)};
  double s = scaled.f;
  double h2 = h*h;
  c0 = __builtin_fma(h, 0x1.ebfbdff84555ap-11, 0x1.62e42fefa398bp-5);
  c2 = __builtin_fma(h, 0x1.3b2ad1b1716a2p-23, 0x1.c6b08d4ad86d3p-17);
  c4 = __builtin_fma(h, 0x1.4a1d7f457ac56p-37, 0x1.5d7472718ce9dp-30);
  c0 = __builtin_fma(h2, __builtin_fma(h2, c4, c2), c0);
  double w = s*h;
  b64u64_u rr = {.f = __builtin_fma(w, c0, s)};
  // with off=467, fails for x,y=0x1.fd12b4p-1,-0x1.d0b058p+13 and RNDZ
  uint64_t off = 468;
  if(((rr.u+off)&0xfffffff) <= 2*off)
    return as_powf_accurate2 (xd.u, e - k, y0);
  float res = rr.f;
  return res;
}



static __attribute__((noinline)) float as_powf_accurate2(uint64_t x_bits, int e, float y0){
  static const double o[] = {1, 2};
  static const double ch[][2] =
    {{0x1.71547652b82fep+1, 0x1.777d0ffda2b89p-55}, {0x1.ec709dc3a03fdp-1, 0x1.d27f04ff73b3ap-55},
     {0x1.2776c50ef9bfep-1, 0x1.e4b514251d0ecp-55}, {0x1.a61762a7aded9p-2, 0x1.de632dc7f6998p-57},
     {0x1.484b13d7c02aep-2, 0x1.a320ec342ddb3p-56}, {0x1.0c9a84993fd48p-2, -0x1.e6425ce9a74a4p-57},
     {0x1.c68f568d8beafp-3, -0x1.03a175487feabp-57}, {0x1.89f3b14657dfbp-3, 0x1.f04a3acf0bcf7p-57},
     {0x1.5b9ad2f2d12ap-3, -0x1.68fdff6815a6fp-58}, {0x1.3702165b88acbp-3, 0x1.45b052ace6c8ep-60},
     {0x1.1998f60f2f005p-3, -0x1.79a94f62fb524p-57}, {0x1.f9bc428e30809p-4, -0x1.51f063387e47p-59},
     {0x1.1ac0ab871296ap-3, 0x1.2ba6a2e1a625bp-57}
    };
  static const double ce[][2] =
    {{0x1p+0, 0x1.f7d70599926c4p-98}, {0x1.62e42fefa39efp-1, 0x1.abc9e3b39856bp-56},
     {0x1.ebfbdff82c58fp-3, -0x1.5e43a540c283dp-57}, {0x1.c6b08d704a0cp-5, -0x1.d3316277451e6p-59},
     {0x1.3b2ab6fba4e77p-7, 0x1.4e66003ba7f85p-62}, {0x1.5d87fe78a6731p-10, 0x1.07183d46a9697p-66},
     {0x1.430912f86c787p-13, 0x1.bc81afca4c93p-67}, {0x1.ffcbfc588b0c7p-17, -0x1.e63f6f0116f4cp-71},
     {0x1.62c0223a5c826p-20, -0x1.30542d98ea4a5p-74}, {0x1.b5253d395e7c6p-24, -0x1.9285a132ce05ep-80},
     {0x1.e4cf5158b7b01p-28, -0x1.9ac1facae1b88p-83}, {0x1.e8cac7351a7a8p-32, -0x1.4fb82adebd76bp-91},
     {0x1.c3bd65182746dp-36, 0x1.84ad0689d30ep-91}, {0x1.8161931d765c3p-40, -0x1.254c6535279cep-95},
     {0x1.314943a26c9e2p-44, -0x1.f4f2fdc14fb82p-98}, {0x1.c36e53b459602p-49, 0x1.f0d06a5a63c41p-103},
     {0x1.397637b3876a4p-53, -0x1.5632c551ae458p-107}, {0x1.98fbfefdddb51p-58, -0x1.fd134923d52b4p-115}};
  double y = y0;
  b64u64_u t = {.u = x_bits};
  int k = (t.u & (~(uint64_t)0 >> 12)) > 0x6a09e667f3bcdull;
  e += k;
  double x = t.f;
  double xm = x-o[k], xp = x+o[k], zh = xm/xp, zl = __builtin_fma(zh,-xp,xm)/xp;
  double z2l, z2h = muldd(zh, zl, zh, zl, &z2l);
  z2h = polydd(z2h,z2l, 13, ch, &z2l);
  zh = muldd(zh,zl, z2h,z2l, &zl);
  zh = mulddd(zh,zl, y, &zl);
  double ey = e*y, eh = ey + zh, el = ((ey - eh) + zh) + zl, ee = roundeven_finite(eh);
  eh -= ee;
  eh = polydd(eh, el, 18, ce, &el);
  b64u64_u r = {.u = ((uint64_t)0x3ff+(int64_t)ee)<<52};
  b64u64_u ll = {.f = el}, lh = {.f = eh};
  if(((ll.u>>(6*4-1))&((1<<29)-1)) == ((1<<29)-1)){
    if(eh<1){
      if(el>=0x1p-54){
	el -= 0x1p-53;
	eh += 0x1p-53;
      } else if(el<=-0x1p-54){
	el += 0x1p-53;
	eh -= 0x1p-53;
      }
    } else {
      if(el>=0x1p-53){
	el -= 0x1p-52;
	eh += 0x1p-52;
      }else if(el<=-0x1p-53){
	el += 0x1p-52;
	eh -= 0x1p-52;
      }
    }
  } else if(((ll.u>>(6*4-1))&((1<<29)-1)) == 0 ){
    if(el>0){
      if(eh<1){
	if(el>=0x1p-53){
	  el -= 0x1p-53;
	  eh += 0x1p-53;
	}
      } else {
	if(el>=0x1p-52){
	  el -= 0x1p-52;
	  eh += 0x1p-52;
	}
      }
    } else {
      if(eh<1){
	if(el<=-0x1p-53){
	  el += 0x1p-53;
	  eh -= 0x1p-53;
	}
      } else {
	if(el<=-0x1p-52){
	  el += 0x1p-52;
	  eh -= 0x1p-52;
	}
      }
    }
  }
  ll.f = el;
  lh.f = eh;
  if((lh.u&0xfffffff) == 0){
    if(__builtin_fabs(ll.f)>0x1p-91){
      if(el<0){
	lh.u--;
	eh = lh.f;
      } else {
	lh.u++;
	eh = lh.f;
      }
    }
  }
  
  eh *= r.f;
  el *= r.f;
  float res = eh;
  return res;
}
