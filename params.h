#ifndef PARAMS_H
#define PARAMS_H

//#define ADDITION_PROOF
//#define MULTIPLICATION_PROOF
//#define MULTIPLICATION_PROOF_2

#define NAMESPACE(s) s

#if defined(ADDITION_PROOF)
#define M 6
#elif defined(MULTIPLICATION_PROOF)
#define M 10
#elif defined(MULTIPLICATION_PROOF_2)
#define M 15
#endif

// Make thread blocks as big as possible to maximize GPU utilization.
#define N 128
#define Q 1073479681
#define GAMMA1 (1 << 18) // 262144
#define GAMMA2 (Q-1)/(1 << 13) // 131040
#define D 14
#define BETA 32
#define R 4
#define K 10
#define L 10

#define SYMBYTES 32

// 1073479681 = 4095 * 2^18 + 1
#define MODULUS 1073479681
#define ROOT 
#define ROOT_INVERSE 
#define ROOT_PW (1 << 18)

// const int MODULUS = 7340033;
// const int ROOT = 2187;//3;//5;
// const int ROOT_INVERSE = 4665133;//2446678;//4404020;
// const int ROOT_PW = 1 << 20;

#endif
