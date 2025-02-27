#ifndef NTT_COMMON_H
#define NTT_COMMON_H

#include <iostream>
#include <random>

#include "consts.h"
#include "math_helpers.cuh"

using namespace std;

std::vector<int> generate_random_ntt_input(unsigned int num_elements) {
  // Set up random number generators
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<> dist(0, MODULUS - 1);
  
  std::vector<int> input;
  for (unsigned int i = 0; i < num_elements; i++) {
    input.push_back(dist(mt));
  }

  if (!isPowerOfTwo(num_elements)) {
    unsigned int newN = pow(2, ceil(log2(num_elements)));
    for (unsigned int i = num_elements; i < newN; i++) {
      input.push_back(0);
    }
  }
  return input;
}

// output will be of size nlogn since we need to store n twiddle factors for each
// of the logn levels of the FFT.
void cpu_compute_twiddle_factors(int* output, unsigned int n, bool invert, int root, int root_inverse, int root_pw, int modulus) {
  int orig_root = invert ? root_inverse : root;
  int counter = 0;
  for (int len = 2; len <= n; len <<= 1) {
    int root_of_unity = orig_root;
    for (int i = len; i < root_pw; i <<= 1) {
      root_of_unity = (int)(1LL * root_of_unity * root_of_unity % modulus);
    }

    for (int i = 0; i < n; i += len) {
      int w = 1;
      for (int j = 0; j < len / 2; j++) {
        output[counter] = w;
        counter++;
        w = (int)(1LL * w * root_of_unity % modulus);
      }
    }
  }
}

// Returns the number of twiddle factors needed for 1 direction.
unsigned int compute_forward_and_inverse_twiddle_factors(int** forward_twiddle_factors, int** inverse_twiddle_factors, unsigned int n) {
  unsigned int power = get_power_of_2(n);
  unsigned int num_twiddle_factors = n * power / 2;

  *forward_twiddle_factors = new int[num_twiddle_factors];
  *inverse_twiddle_factors = new int[num_twiddle_factors];
  cpu_compute_twiddle_factors(*forward_twiddle_factors, n, false, ROOT, ROOT_INVERSE, ROOT_PW, MODULUS);
  cpu_compute_twiddle_factors(*inverse_twiddle_factors, n, true, ROOT, ROOT_INVERSE, ROOT_PW, MODULUS);
  return num_twiddle_factors;
}

#endif
