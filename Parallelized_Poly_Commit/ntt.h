#ifndef NTT_H
#define NTT_H

#include <stdint.h>

typedef uint64_t goldilocks_t;

// NTT function declaration
void ntt(goldilocks_t* a, uint32_t m, uint32_t n, bool is_inverse);

// Helper function declarations
goldilocks_t goldilocks_root(size_t size);
goldilocks_t goldilocks_add(goldilocks_t a, goldilocks_t b);
goldilocks_t goldilocks_sub(goldilocks_t a, goldilocks_t b);
goldilocks_t goldilocks_mul(goldilocks_t a, goldilocks_t b);
goldilocks_t goldilocks_inv(goldilocks_t a);

#endif // NTT_H 
