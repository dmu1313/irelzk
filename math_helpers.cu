
#include "math_helpers.cuh"

// Serial algorithms borrowed from various sources:
// https://cp-algorithms.com/algebra/fft.html
// https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/
// https://github.com/wandering007/algorithms/blob/master/extended_gcd.cpp
// https://cp-algorithms.com/algebra/extended-euclid-algorithm.html
// https://cp-algorithms.com/algebra/module-inverse.html

__host__ __device__ int mod_inverse(const int a, int modulus) {
  int x;
  int y;
  int g = extended_euclidean(a, modulus, x, y);
  if (g != 1) {
    return -1;
  } else {
    // Handles negative numbers with the additional logic.
    return (x % modulus + modulus) % modulus;
  }
}

__host__ __device__ int extended_euclidean(int a, int b, int& x, int& y) {
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

__host__ __device__ bool isPowerOfTwo(unsigned int n) {
  return (n & (n - 1)) == 0;
}

// given an integer, get the power of 2 that it is
__host__ unsigned int get_power_of_2(unsigned int n) {
  unsigned int n_power_of_2 = 0;
  while (n) {
    n_power_of_2++;
    n >>= 1;
  }
  n_power_of_2--;
  return n_power_of_2;
}
