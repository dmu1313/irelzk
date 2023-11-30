#include "debug.h"

void print_digest(unsigned char *digest, unsigned int digest_length) {
  fprintf_digest(stdout, digest, digest_length);
}

void print_digest_ln(unsigned char *digest, unsigned int digest_length) {
  fprintf_digest_ln(stdout, digest, digest_length);
}

void fprintf_digest(FILE *f, unsigned char *digest, unsigned int digest_length) {
  for(int i = 0; i < digest_length; i++)
    fprintf(f, "%02x", digest[i]);
}

void fprintf_digest_ln(FILE *f, unsigned char *digest, unsigned int digest_length) {
  fprintf_digest(f, digest, digest_length);
  fprintf(f, "\n");
}
