#include <stdint.h>
#include <string.h>
#include "fips202.h"
#include "randombytes.h"
#include "params.h"
#include "addition.h"
#include "opening.h"
#include "product.h"
#include "linear.h"
#include "comm.h"
#include "polyvec.h"
#include "poly.h"

void addition_proof(proof *p, comm *t, commrnd *r, const uint8_t rho[SYMBYTES],
                    const uint64_t a[2], const uint64_t b[2])
{
  int i, rej;
  uint16_t nonce = 0;
  uint32_t f, g, x = 0;
  uint8_t symbuf[2*SYMBYTES+SHAKE128_RATE];
  uint8_t *seed = symbuf;
  uint8_t *thash = symbuf+SYMBYTES;
  uint8_t *chash = symbuf+2*SYMBYTES;
  commkey ck;
  polyvecm msg;
  polyveck t0low, w1[R];
  poly v, vprime[R], tmp;
  keccak_state kecst;
  kecst.mdctx = NULL;

  // Get random bytes as a seed for future use
  randombytes(seed,SYMBYTES);

  // Convert a and b to Two's Complement
  // and store in msg.vec[0] and msg.vec[1]
  // Store a + b in msg.vec[2] (Let c = a + b)
  // Store the carry bit in msg.vec[3]
  for(i=0;i<64;i++) {
    f = (a[0] >> i) & 1;
    g = (b[0] >> i) & 1;
    x = x + f + g;
    msg.vec[0].coeffs[i] = f;
    msg.vec[1].coeffs[i] = g;
    msg.vec[2].coeffs[i] = x & 1;
    x >>= 1;
    msg.vec[3].coeffs[i] = x;
  }
  for(i=0;i<64;i++) {
    f = (a[1] >> i) & 1;
    g = (b[1] >> i) & 1;
    x = x + f + g;
    msg.vec[0].coeffs[64+i] = f;
    msg.vec[1].coeffs[64+i] = g;
    msg.vec[2].coeffs[64+i] = x & 1;
    x >>= 1;
    msg.vec[3].coeffs[64+i] = x;
  }

  memset(&msg.vec[M-2],0,sizeof(poly)); // Zero out the 2nd to last polynomial
  poly_uniform(&msg.vec[M-1],seed,nonce++); // Get uniformly random polynomial coefficients.
  // Not sure why they zero out R=4 of the coefficients of the last polynomial.
  for(i=0;i<R;i++)
    msg.vec[M-1].coeffs[i] = 0;
  poly_ntt(&msg.vec[M-1]); // Convert last polynomial to NTT domain.

  expand_commkey(&ck,rho); // Init commitment key.
  commit(t,r,&msg,&ck); // Commit to the message containing the numbers being summed and their result.
  polyveck_invntt(&t->t0); // Inverse NTT
  polyveck_power2round(&t->t0,&t0low,&t->t0); // Perform rounding on `t`
  polyveck_ntt(&t->t0); // Convert to NTT
  polyveck_ntt(&t0low); // Convert to NTT

  opening_init(&t0low,r,&ck); // Sets up variables necessary in the opening proof.

  // Do lots of hashing of `t`.
  shake128_init(&kecst);
  shake128_absorb(&kecst,rho,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&t->t0,sizeof(polyveck));
  shake128_absorb(&kecst,(uint8_t *)&t->tm,(M-2)*sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)&t->tm.vec[M-1],sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeeze(thash,SYMBYTES,&kecst);
  do {
    opening_first(w1,seed,nonce); // Fiat shamir first part
    nonce += R*(K+L+M);

    shake128_init(&kecst);
    shake128_absorb(&kecst,thash,SYMBYTES);
    shake128_absorb(&kecst,(uint8_t *)w1,R*sizeof(polyveck));
    shake128_finalize(&kecst);
    shake128_squeezeblocks(chash,1,&kecst);


    product(&v,&msg,chash);
    poly_add(&tmp,&t->tm.vec[M-2],&msg.vec[M-2]);
    linear(vprime,&p->h,&msg,chash+SYMBYTES);

    shake128_init(&kecst);
    shake128_absorb(&kecst,chash,2*SYMBYTES);
    shake128_absorb(&kecst,(uint8_t *)&tmp,sizeof(poly));
    shake128_absorb(&kecst,(uint8_t *)&v,sizeof(poly));
    shake128_absorb(&kecst,(uint8_t *)&p->h,sizeof(poly));
    shake128_absorb(&kecst,(uint8_t *)vprime,R*sizeof(poly));
    shake128_finalize(&kecst);
    shake128_squeezeblocks(chash,1,&kecst);

    challenge_prehash(p->c,chash);
    rej = opening_last(p->z,p->c,w1); // Fiat shamir second part with responses `z`
  } while(rej);

  t->tm.vec[M-2] = tmp;
}

int addition_proof_verify(const proof *p, const comm *t, const uint8_t rho[SYMBYTES]) {
  uint8_t symbuf[SYMBYTES+SHAKE128_RATE];
  uint8_t *thash = symbuf;
  uint8_t *chash = symbuf+SYMBYTES;
  commkey ck;
  poly v, vprime[R];
  polyveck w1[R];
  keccak_state kecst;
  kecst.mdctx = NULL;

  shake128_init(&kecst);
  shake128_absorb(&kecst,rho,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&t->t0,sizeof(polyveck));
  shake128_absorb(&kecst,(uint8_t *)&t->tm,R*sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)&t->tm.vec[M-1],sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeeze(thash,SYMBYTES,&kecst);

  expand_commkey(&ck,rho);
  if(opening_verify_first(w1,p->c,p->z,t,&ck))
    return 1;

  shake128_init(&kecst);
  shake128_absorb(&kecst,thash,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)w1,R*sizeof(polyveck));
  shake128_finalize(&kecst);
  shake128_squeezeblocks(chash,1,&kecst);

  if(product_verify(&v,chash,p->c,p->z,t,&ck))
    return 1;
  if(linear_verify(vprime,chash+SYMBYTES,&p->h,p->c,p->z,t,&ck))
    return 1;

  shake128_init(&kecst);
  shake128_absorb(&kecst,chash,2*SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&t->tm.vec[M-2],sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)&v,sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)&p->h,sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)vprime,R*sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeezeblocks(chash,1,&kecst);

  if(opening_verify_last(p->c,chash))
    return 1;

  return 0;
}
