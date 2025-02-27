
#include <algorithm>

#include "ntt.cuh"

// #define DEBUG

#ifdef DEBUG
#define DEBUG_PRINT(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...) do {} while (0)
#endif

__device__ int extended_euclidean_device(int a, int b, int& x, int& y) {
  x = 1;
  y = 0;
  int new_x = 0;
  int new_y = 1;
  int new_gcd = b;
  int gcd = a;
  while (new_gcd) {
    int quotient = gcd / new_gcd;

    int temp = gcd;
    gcd = new_gcd;
    new_gcd = temp - quotient * new_gcd;

    temp = x;
    x = new_x;
    new_x = temp - quotient * new_x;

    temp = y;
    y = new_y;
    new_y = temp - quotient * new_y;
  }
  return gcd;
}

__device__ int mod_inverse_device(const int a, int modulus) {
  int x;
  int y;
  int g = extended_euclidean_device(a, modulus, x, y);
  if (g != 1) {
    return -1;
  } else {
    // Handles negative numbers with the additional logic.
    return (x % modulus + modulus) % modulus;
  }
}

__device__ inline bool in_range(int lower_bound, int upper_bound, int index) {
  return index >= lower_bound && index <= upper_bound;
}

__global__ void bit_reverse_positions_kernel(int* input, unsigned int n) {
  unsigned int id = threadIdx.x + blockDim.x * blockIdx.x;
  if (id >= n) {
    return;
  }

  // get the power of 2 that n is
  unsigned int n_power_of_2 = 0;
  unsigned int temp = n;
  while (temp) {
    n_power_of_2++;
    temp >>= 1;
  }
  n_power_of_2--;

  unsigned int reversed_id = 0;
  unsigned int select_bit = 1;
  for (int i = 0; i < n_power_of_2; i++) {
    unsigned int bit = id & select_bit;

    reversed_id |= (bit << ((n_power_of_2 - 1) - i)) >> i;
    select_bit <<= 1;
  }
  if (id < reversed_id) {
    int temp = input[id];
    input[id] = input[reversed_id];
    input[reversed_id] = temp;
  }
}

__global__ void ntt_kernel(
  int* input,
  unsigned int n,
  bool invert,
  int mod,
  int root,
  int root_inverse,
  int root_pw,
  int* twiddle_factors
) {
  extern __shared__ unsigned char sharedMem[];
  unsigned int id = threadIdx.x + blockDim.x * blockIdx.x;
  if (id >= n) {
    return;
  }
  
  for (int len = 2, level = 0; len <= n; len <<= 1, level++) {

    int butterfly_index = threadIdx.x;
    int shamt = level;
    shamt = shamt < 0 ? 0 : shamt;
    int offset_multiplier = threadIdx.x >> shamt;
    int offset = offset_multiplier * len/2;
    butterfly_index += offset;

    // if (threadIdx.x < n) {
    //   printf("butterfly_index=%d, threadIdx.x=%d\n", butterfly_index, threadIdx.x);
    // }

    int num_rounds_per_butterfly_group = (len>>1) / blockDim.x;
    int butterfly_index_increase = num_rounds_per_butterfly_group > 1 ? blockDim.x : (blockDim.x * 2);
    int num_rounds_for_current_butterfly_group = 0;
    // printf("num_rounds_per_butterfly_group=%d, len>>1=%d, blockDim.x=%d\n", num_rounds_per_butterfly_group, len>>1, blockDim.x);
    int num_rounds = (n/2 + blockDim.x - 1) / blockDim.x;
    for (int i = 0; i < num_rounds; i++) {
      if (butterfly_index >= n) {
        break;
      }
      // int twiddle_layer_offset = (butterfly_index + 1) / 2;
      int twiddle_layer_offset = blockDim.x * i + threadIdx.x; // 
      int twiddle_index = (level*n/2)+twiddle_layer_offset;
      int w = twiddle_factors[twiddle_index];
      // printf(
      //   "twiddle_index=%d, butterfly_index=%d, level=%d, twiddle_layer_offset=%d, num_rounds=%d, i=%d, threadIdx.x=%d\n",
      //   twiddle_index, butterfly_index, level, twiddle_layer_offset, num_rounds, i, threadIdx.x);

      int u = input[butterfly_index];
      int v = (int)(1LL * input[butterfly_index + len / 2] * w % mod);

      input[butterfly_index] = u+v < mod ? u + v : u + v - mod;
      input[butterfly_index + len / 2] = u-v >= 0 ? u-v : u-v+mod;

      butterfly_index += butterfly_index_increase;
      num_rounds_for_current_butterfly_group++;
      if (num_rounds_for_current_butterfly_group == num_rounds_per_butterfly_group &&
          num_rounds_for_current_butterfly_group > 1) {
        butterfly_index += len>>1;
        num_rounds_for_current_butterfly_group = 0;
      }
    }
    __syncthreads();
  }
}

__global__ void ntt_inverse_final(int* input, unsigned int n, int mod) {
  unsigned int id = threadIdx.x + blockDim.x * blockIdx.x;
  if (id < n) {
    int n_inverse = mod_inverse_device((int)n, mod);
    input[id] = (int)(1LL * input[id] * n_inverse % mod);
  }
}

__host__ void ntt(
  int * input,
  unsigned int n,
  int* d_twiddle_factors,
  bool invert,
  int mod,
  int root,
  int root_inverse,
  int root_pw
) {
  int threads_per_block = 1024;
  int num_blocks = (n + threads_per_block - 1) / threads_per_block;
  bit_reverse_positions_kernel<<<num_blocks, threads_per_block>>>(input, n);
  ntt_kernel<<<1, threads_per_block>>>(
    input, n, invert, mod, root, root_inverse, root_pw, d_twiddle_factors);
  if (invert) {
    ntt_inverse_final<<<num_blocks, threads_per_block>>>(input, n, mod);
  }
}
