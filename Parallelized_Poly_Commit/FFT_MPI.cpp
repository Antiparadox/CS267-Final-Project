#include <mpi.h>
#include <vector>
#include <cstdint>
#include <cmath>
#include <iostream>
#include "ntt_mpi.h"

// Helper function to compute modular exponentiation
goldilocks_t mod_pow(goldilocks_t base, goldilocks_t exp) {
    goldilocks_t result = 1;
    while (exp > 0) {
        if (exp & 1) {
            result = goldilocks_mul_mpi(result, base);
        }
        base = goldilocks_mul_mpi(base, base);
        exp >>= 1;
    }
    return result;
}

// Goldilocks field arithmetic operations
goldilocks_t goldilocks_add_mpi(goldilocks_t a, goldilocks_t b) {
    goldilocks_t sum = a + b;
    if (sum >= GOLDILOCKS_PRIME) {
        sum -= GOLDILOCKS_PRIME;
    }
    return sum;
}

goldilocks_t goldilocks_sub_mpi(goldilocks_t a, goldilocks_t b) {
    if (a < b) {
        a += GOLDILOCKS_PRIME;
    }
    return a - b;
}

goldilocks_t goldilocks_mul_mpi(goldilocks_t a, goldilocks_t b) {
    __uint128_t prod = (__uint128_t)a * b;
    return prod % GOLDILOCKS_PRIME;
}

goldilocks_t goldilocks_inv_mpi(goldilocks_t a) {
    return mod_pow(a, GOLDILOCKS_PRIME - 2);
}

// Helper function to compute bit reversal
uint32_t bit_reverse(uint32_t x, uint32_t bits) {
    uint32_t y = 0;
    for (uint32_t i = 0; i < bits; i++) {
        y = (y << 1) | (x & 1);
        x >>= 1;
    }
    return y;
}

// Helper function to compute modular inverse
goldilocks_t mod_inv(goldilocks_t a) {
    return mod_pow(a, GOLDILOCKS_PRIME - 2);
}

// Local 1D FFT computation
void local_fft(goldilocks_t* data, uint32_t size, bool is_inverse) {
    // Bit-reverse permutation
    for (uint32_t i = 0; i < size; i++) {
        uint32_t j = bit_reverse(i, static_cast<uint32_t>(log2(size)));
        if (j > i) {
            std::swap(data[i], data[j]);
        }
    }

    // Cooley-Tukey FFT
    for (uint32_t len = 2; len <= size; len <<= 1) {
        goldilocks_t wlen = is_inverse ? 
            mod_pow(ROOT_OF_UNITY, GOLDILOCKS_PRIME - 1 - (GOLDILOCKS_PRIME - 1) / len) :
            mod_pow(ROOT_OF_UNITY, (GOLDILOCKS_PRIME - 1) / len);

        for (uint32_t i = 0; i < size; i += len) {
            goldilocks_t w = 1;
            for (uint32_t j = 0; j < len/2; j++) {
                goldilocks_t u = data[i + j];
                goldilocks_t t = goldilocks_mul_mpi(w, data[i + j + len/2]);
                data[i + j] = goldilocks_add_mpi(u, t);
                data[i + j + len/2] = goldilocks_sub_mpi(u, t);
                w = goldilocks_mul_mpi(w, wlen);
            }
        }
    }

    // Scale by 1/N for inverse FFT
    if (is_inverse) {
        goldilocks_t n_inv = mod_inv(size);
        for (uint32_t i = 0; i < size; i++) {
            data[i] = goldilocks_mul_mpi(data[i], n_inv);
        }
    }
}

// Transpose-based parallel FFT implementation
void ntt_mpi(goldilocks_t* data, uint32_t m, uint32_t n, bool is_inverse, int rank, int size) {
    uint32_t total_size = 1 << (m + n);
    uint32_t block_size = total_size / size;
    uint32_t total_order = m + n;
    
    // Ensure the number of processes is a power of 2
    if ((size & (size - 1)) != 0) {
        if (rank == 0) {
            std::cerr << "Error: Number of processes must be a power of 2" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Ensure block_size is a power of 2
    if ((block_size & (block_size - 1)) != 0) {
        if (rank == 0) {
            std::cerr << "Error: Block size must be a power of 2" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Find primitive root of unity
    goldilocks_t omega = 7;  // 7 is a generator of the multiplicative group
    goldilocks_t order = (GOLDILOCKS_PRIME - 1) >> total_order;
    omega = mod_pow(omega, order);
    if (is_inverse) {
        omega = mod_pow(omega, GOLDILOCKS_PRIME - 2);
    }

    // Phase 1: Bit-reverse permutation
    for (uint32_t i = 0; i < block_size; i++) {
        uint32_t j = 0;
        for (uint32_t k = 0; k < static_cast<uint32_t>(log2(block_size)); k++) {
            j = (j << 1) | ((i >> k) & 1);
        }
        if (j > i) {
            std::swap(data[i], data[j]);
        }
    }

    // Phase 2: Local FFT on rows
    for (uint32_t s = 1; s <= static_cast<uint32_t>(log2(block_size)); s++) {
        uint32_t m_s = 1 << s;
        goldilocks_t w_m = mod_pow(omega, 1 << (total_order - s));
        
        for (uint32_t i = 0; i < block_size; i += m_s) {
            goldilocks_t w = 1;
            for (uint32_t j = 0; j < m_s/2; j++) {
                goldilocks_t u = data[i + j];
                goldilocks_t t = goldilocks_mul_mpi(w, data[i + j + m_s/2]);
                data[i + j] = goldilocks_add_mpi(u, t);
                data[i + j + m_s/2] = goldilocks_sub_mpi(u, t);
                w = goldilocks_mul_mpi(w, w_m);
            }
        }
    }

    // Phase 3: Transpose
    std::vector<goldilocks_t> sendbuf(block_size);
    std::vector<goldilocks_t> recvbuf(block_size);

    // Prepare data for transpose
    for (uint32_t i = 0; i < block_size; i++) {
        sendbuf[i] = data[i];
    }

    // Calculate and print total communication size
    int send_count = block_size / size;
    size_t total_comm_size = (size_t)send_count * size * sizeof(goldilocks_t);
    
    // Measure communication time
    double comm_start = MPI_Wtime();
    
    // Perform all-to-all communication
    MPI_Alltoall(sendbuf.data(), send_count, MPI_UINT64_T,
                 recvbuf.data(), send_count, MPI_UINT64_T,
                 MPI_COMM_WORLD);
    
    double comm_end = MPI_Wtime();
    double comm_time = comm_end - comm_start;
    
    // Calculate communication bandwidth
    double bandwidth = total_comm_size * size / (comm_time * 1e6); // MB/s
    
    if (rank == 0) {
        printf("\nTranspose Communication Statistics:\n");
        printf("  Communication size: %zu bytes per process, %zu bytes total across all processes\n", 
               total_comm_size, total_comm_size * size);
        printf("  Communication time: %.6f seconds\n", comm_time);
        printf("  Effective bandwidth: %.2f MB/s\n", bandwidth);
        printf("  Average latency per message: %.6f seconds\n", comm_time / (size - 1));
    }

    // Phase 4: Second local FFT
    for (uint32_t s = 1; s <= static_cast<uint32_t>(log2(block_size)); s++) {
        uint32_t m_s = 1 << s;
        goldilocks_t w_m = mod_pow(omega, 1 << (total_order - s));
        
        for (uint32_t i = 0; i < block_size; i += m_s) {
            goldilocks_t w = 1;
            for (uint32_t j = 0; j < m_s/2; j++) {
                goldilocks_t u = recvbuf[i + j];
                goldilocks_t t = goldilocks_mul_mpi(w, recvbuf[i + j + m_s/2]);
                recvbuf[i + j] = goldilocks_add_mpi(u, t);
                recvbuf[i + j + m_s/2] = goldilocks_sub_mpi(u, t);
                w = goldilocks_mul_mpi(w, w_m);
            }
        }
    }

    // Phase 5: Final bit-reverse permutation
    for (uint32_t i = 0; i < block_size; i++) {
        uint32_t j = 0;
        for (uint32_t k = 0; k < static_cast<uint32_t>(log2(block_size)); k++) {
            j = (j << 1) | ((i >> k) & 1);
        }
        if (j > i) {
            std::swap(recvbuf[i], recvbuf[j]);
        }
    }

    // Copy result back to original array
    for (uint32_t i = 0; i < block_size; i++) {
        data[i] = recvbuf[i];
    }

    // Scale by 1/N for inverse FFT
    if (is_inverse) {
        goldilocks_t n_inv = mod_inv(total_size);
        for (uint32_t i = 0; i < block_size; i++) {
            data[i] = goldilocks_mul_mpi(data[i], n_inv);
        }
    }
} 
