
#include <stdio.h>

#define DEBUG

#ifdef DEBUG
#define DEBUG_PRINT(fmt, args...) printf(fmt, ##args); fflush(stdout)
#define DEBUG_PRINT_DIGEST(digest, digest_length) print_digest(digest, digest_length); fflush(stdout);
#define DEBUG_PRINT_DIGEST_LN(digest, digest_length) print_digest_ln(digest, digest_length); fflush(stdout);
#else
#define DEBUG_PRINT(fmt, args...) /* do nothing */
#define DEBUG_PRINT_DIGEST(digest, digest_length) /* do nothing */
#define DEBUG_PRINT_DIGEST_LN(digest, digest_length) /* do nothing */
#endif


#ifndef DEBUG_H
#define DEBUG_H

void print_digest(unsigned char *digest, unsigned int digest_length);
void print_digest_ln(unsigned char *digest, unsigned int digest_length);
void fprintf_digest(FILE *f, unsigned char *digest, unsigned int digest_length);
void fprintf_digest_ln(FILE *f, unsigned char *digest, unsigned int digest_length);

#endif
