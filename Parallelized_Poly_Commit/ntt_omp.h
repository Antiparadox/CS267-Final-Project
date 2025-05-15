#ifndef NTT_OMP_H
#define NTT_OMP_H

#include <stdint.h>

typedef uint64_t goldilocks_t;

// NTT function declaration for OpenMP implementation
void ntt_omp(goldilocks_t* a, uint32_t m, uint32_t n, bool is_inverse);

// Helper function declarations
goldilocks_t goldilocks_root_omp(size_t size);
goldilocks_t goldilocks_add_omp(goldilocks_t a, goldilocks_t b);
goldilocks_t goldilocks_sub_omp(goldilocks_t a, goldilocks_t b);
goldilocks_t goldilocks_mul_omp(goldilocks_t a, goldilocks_t b);
goldilocks_t goldilocks_inv_omp(goldilocks_t a);

#endif // NTT_OMP_H 