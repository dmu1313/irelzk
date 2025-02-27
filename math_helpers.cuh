#ifndef MATH_HELPERS_H
#define MATH_HELPERS_H

#include <cuda_runtime.h>
#include <vector>
#include <cstdint>

using namespace std;

__host__ __device__ int mod_inverse(const int a, int modulus);

__host__ __device__ int extended_euclidean(int a, int b, int& x, int& y);

__host__ __device__ bool isPowerOfTwo(unsigned int n);

__host__ unsigned int get_power_of_2(unsigned int n);

#endif
