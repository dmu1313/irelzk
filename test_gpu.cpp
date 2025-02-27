
#include <iostream>
#include <immintrin.h>

#include "consts.h"
// #include "poly.h"
#include "ntt.h"
#include "gpu_ntt.cuh"

using namespace std;

typedef struct {
  int32_t coeffs[N];
} ppoly __attribute__((aligned(32)));

void poly_reduce(ppoly *r) {
  int i;
  __m256i f,t;
  const __m256i mask = _mm256_set1_epi32((1 << 30) - 1);

  for(i=0;i<N/8;i++) {
    f = _mm256_load_si256((__m256i *)&r->coeffs[8*i]);
    t = _mm256_srai_epi32(f,30);
    f = _mm256_and_si256(f,mask);
    f = _mm256_sub_epi32(f,t);
    t = _mm256_slli_epi32(t,18);
    f = _mm256_add_epi32(f,t);
    _mm256_store_si256((__m256i *)&r->coeffs[8*i],f);
  }
}

void poly_scale_montgomery(ppoly *r, const ppoly *a, int32_t s)
{
  int i;
  __m256i f0,f1,g0,g1;
  const __m256i q = _mm256_load_si256((__m256i *)&qdata[_8XQ]);
  const __m256i lo = _mm256_set1_epi32(s*QINV);
  const __m256i hi = _mm256_set1_epi32(s);

  for(i=0;i<N;i+=8) {
    f0 = _mm256_load_si256((__m256i *)&a->coeffs[i]);
    f1 = (__m256i)_mm256_movehdup_ps((__m256)f0);
    g0 = _mm256_mul_epi32(f0,lo);
    g1 = _mm256_mul_epi32(f1,lo);
    f0 = _mm256_mul_epi32(f0,hi);
    f1 = _mm256_mul_epi32(f1,hi);
    g0 = _mm256_mul_epi32(g0,q);
    g1 = _mm256_mul_epi32(g1,q);
    f0 = _mm256_sub_epi32(f0,g0);
    f1 = _mm256_sub_epi32(f1,g1);
    f0 = (__m256i)_mm256_movehdup_ps((__m256)f0);
    f0 = _mm256_blend_epi32(f0,f1,0xAA);
    _mm256_store_si256((__m256i *)&r->coeffs[i],f0);
  }
}

void poly_ntt(ppoly *r) {
  ntt_avx(r->coeffs,qdata);
  poly_reduce(r);
}

void poly_invntt(ppoly *r) {
  poly_scale_montgomery(r,r,33554432);
  invntt_avx(r->coeffs,qdata);
}

int main(int argc, char *argv[]) {
    ppoly p;
    for (int i = 0; i < 128; i++) {
        p.coeffs[i] = i;
    }

    testFunc();
    
    for (int i = 0; i < 128; i++) {
        cout << p.coeffs[i] << " ";
    }
    printf("\n");
    printf("\n");

    poly_ntt(&p);
    for (int i = 0; i < 128; i++) {
        cout << p.coeffs[i] << " ";
    }
    printf("\n");
    printf("\n");

    poly_invntt(&p);
    for (int i = 0; i < 128; i++) {
        cout << p.coeffs[i] << " ";
    }
    printf("\n");
    printf("\n");

    return 0;







    int32_t input[128];
    for (int i = 0; i < 128; i++) {
        input[i] = i;
    }

    testFunc();
    
    for (int i = 0; i < 128; i++) {
        cout << input[i] << " ";
    }
    printf("\n");
    printf("\n");

    ntt_avx(input, qdata);
    for (int i = 0; i < 128; i++) {
        cout << input[i] << " ";
    }
    printf("\n");
    printf("\n");

    invntt_avx(input, qdata);

    for (int i = 0; i < 128; i++) {
        cout << input[i] << " ";
    }
    printf("\n");
    printf("\n");
    
    return 0;
}
