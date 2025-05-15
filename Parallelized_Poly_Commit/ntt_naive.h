#ifndef NTT_NAIVE_H
#define NTT_NAIVE_H

#include <stdint.h>

typedef uint64_t goldilocks_t;

// NTT function declaration for naive implementation
void ntt_naive(goldilocks_t* a, uint32_t m, uint32_t n, bool is_inverse);

// Helper function declarations
goldilocks_t goldilocks_root_naive(size_t size);
goldilocks_t goldilocks_add_naive(goldilocks_t a, goldilocks_t b);
goldilocks_t goldilocks_sub_naive(goldilocks_t a, goldilocks_t b);
goldilocks_t goldilocks_mul_naive(goldilocks_t a, goldilocks_t b);
goldilocks_t goldilocks_inv_naive(goldilocks_t a);

#endif // NTT_NAIVE_H 