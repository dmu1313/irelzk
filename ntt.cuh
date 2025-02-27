
#ifndef NTT_H
#define NTT_H

#include <vector>
#include <cstdint>

#define Q 1073479681

__global__ void ntt_kernel(
  int* input,
  unsigned int n,
  bool invert,
  int mod,
  int root,
  int root_inverse,
  int root_pw,
  int* twiddle_factors,
  int len,
  int level);

__host__ void ntt(
  int * input,
  unsigned int n,
  int* d_twiddle_factors,
  bool invert,
  int mod,
  int root,
  int root_inverse,
  int root_pw);

#endif
