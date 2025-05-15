#ifndef NTT_MPI_H
#define NTT_MPI_H

#include <cstdint>

// Define the Goldilocks field type
typedef uint64_t goldilocks_t;

// Goldilocks field constants
constexpr goldilocks_t GOLDILOCKS_PRIME = 0xFFFFFFFF00000001;
constexpr goldilocks_t ROOT_OF_UNITY = 7;

// Function declarations for MPI NTT implementation
void ntt_mpi(goldilocks_t* data, uint32_t m, uint32_t n, bool is_inverse, int rank, int size);

// Helper function declarations
goldilocks_t goldilocks_add_mpi(goldilocks_t a, goldilocks_t b);
goldilocks_t goldilocks_sub_mpi(goldilocks_t a, goldilocks_t b);
goldilocks_t goldilocks_mul_mpi(goldilocks_t a, goldilocks_t b);
goldilocks_t goldilocks_inv_mpi(goldilocks_t a);

#endif // NTT_MPI_H 
