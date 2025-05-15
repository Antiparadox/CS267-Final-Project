#ifndef GOLDILOCKS_H
#define GOLDILOCKS_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "ntt_naive.h"

typedef uint64_t goldilocks_t;

// Helper function for min
size_t min(size_t a, size_t b) {
    return a < b ? a : b;
}

// Forward declarations
goldilocks_t goldilocks_root_naive(size_t size);

static const goldilocks_t GOLDILOCKS_P = 0xffffffff00000001ULL;

// a + b mod p
inline goldilocks_t goldilocks_add_naive(goldilocks_t a, goldilocks_t b) {
    goldilocks_t s = a + b;
    if (s < a || s >= GOLDILOCKS_P) s -= GOLDILOCKS_P;
    return s;
}

// a − b mod p
inline goldilocks_t goldilocks_sub_naive(goldilocks_t a, goldilocks_t b) {
    return (a >= b) ? (a - b) : (GOLDILOCKS_P - (b - a));
}

// reduce a 128‑bit value x mod p
static inline goldilocks_t goldilocks_reduce128(__uint128_t x) {
    // Since p = 2^64 - 2^32 + 1
    // We can write x = a*2^64 + b
    uint64_t a = (uint64_t)(x >> 64);
    uint64_t b = (uint64_t)x;
    
    // First reduce a*2^64 mod p
    // 2^64 ≡ 2^32 - 1 (mod p)
    // So a*2^64 ≡ a*2^32 - a (mod p)
    __uint128_t t = (__uint128_t)a << 32;  // a*2^32
    t = t - a + b;  // a*2^32 - a + b
    
    // Now t might be up to ~2^96, so we need to reduce it further
    // We can use the same trick again for the high bits
    uint64_t hi = (uint64_t)(t >> 64);
    uint64_t lo = (uint64_t)t;
    
    if (hi > 0) {
        // hi*2^64 ≡ hi*2^32 - hi (mod p)
        __uint128_t t2 = (__uint128_t)hi << 32;
        t2 = t2 - hi + lo;
        lo = (uint64_t)t2;
        
        // One more reduction might be needed
        if (t2 >> 64) {
            lo += (uint64_t)(t2 >> 64) << 32;
            lo -= (uint64_t)(t2 >> 64);
        }
    }
    
    // Final reduction
    while (lo >= GOLDILOCKS_P) {
        lo -= GOLDILOCKS_P;
    }
    
    return lo;
}

// a * b mod p
inline goldilocks_t goldilocks_mul_naive(goldilocks_t a, goldilocks_t b) {
    if (a == 0 || b == 0) return 0;
    return goldilocks_reduce128((__uint128_t)a * b);
}

// −a mod p
static inline goldilocks_t goldilocks_neg(goldilocks_t a) {
    return a ? (GOLDILOCKS_P - a) : 0;
}

// a^e mod p (binary exponentiation)
static inline goldilocks_t goldilocks_pow(goldilocks_t a, uint64_t e) {
    if (a == 0) {
        if (e == 0) return 1;  // 0^0 is undefined, but we return 1 by convention
        return 0;  // 0^e = 0 for e > 0
    }
    
    goldilocks_t res = 1;
    goldilocks_t base = a;
    while (e) {
        if (e & 1) res = goldilocks_mul_naive(res, base);
        base = goldilocks_mul_naive(base, base);
        e >>= 1;
    }
    return res;
}

// a^(−1) mod p using extended Euclidean algorithm
// optimized inverse modulo p = 2^64 - 2^32 + 1
inline goldilocks_t goldilocks_inv_naive(goldilocks_t a) {
    if (a == 0) return 0;  // zero has no inverse; choose 0 by convention

    // initialize
    uint64_t u = GOLDILOCKS_P;
    uint64_t v = a;
    uint64_t t0 = 0;
    uint64_t t1 = 1;

    // pull out all factors of two from v
    int c = __builtin_ctzll(v);
    v >>= c;
    // track total "twos" removed; we also add 96 for the fact that
    // swapping u<->v effectively multiplies by -1 ≡ 2^96
    uint64_t twos = (uint64_t)c + 96;

    // binary‑GCD style loop
    while (u != v) {
        if (u > v) {
            u -= v;
            t0 += t1;
            c = __builtin_ctzll(u);
            u >>= c;
            t1 <<= c;
            twos += c;
        } else {
            v -= u;
            t1 += t0;
            c = __builtin_ctzll(v);
            v >>= c;
            t0 <<= c;
            twos += c;
        }
    }

    // now u == v == 1, and t0·a ≡ 2^twos (mod p)
    // so a^{-1} ≡ t0 · 2^{-twos} ≡ t0 · 2^{191·twos mod 192}
    int k = (int)((191 * twos) % 192);

    // multiply by 2^k via k successive doublings mod p
    goldilocks_t res = t0;
    for (int i = 0; i < k; i++) {
        res = goldilocks_add_naive(res, res);
    }
    return res;
}

// Find a primitive root of unity of order 2^k
static inline goldilocks_t find_primitive_root(uint32_t k) {
    // The multiplicative group has order p-1 = 2^32 * (2^32 - 1)
    // We need a root of order 2^k where k < 32
    goldilocks_t g = 7;  // 7 is a generator of the multiplicative group
    goldilocks_t order = (GOLDILOCKS_P - 1) >> k;  // Divide by 2^k to get element of order 2^k
    return goldilocks_pow(g, order);
}

// In-place NTT implementation using Cooley-Tukey algorithm
// Input: vector of size 2^m, output: NTT of the vector over a subgroup of order 2^(m+n)
// is_inverse: if true, computes inverse NTT
void ntt_naive(goldilocks_t* a, uint32_t m, uint32_t n, bool is_inverse) {
    uint32_t poly_size = 1 << m;
    uint32_t eval_size = 1 << (m + n);
    uint32_t total_order = m + n;
    
    // Find primitive root of unity
    goldilocks_t omega = find_primitive_root(total_order);
    if (is_inverse) {
        omega = goldilocks_inv_naive(omega);
    }
    
    // Bit-reverse permutation of the polynomial coefficients
    for (uint32_t i = 0; i < poly_size; i++) {
        uint32_t j = 0;
        for (uint32_t k = 0; k < m; k++) {
            j = (j << 1) | ((i >> k) & 1);
        }
        if (j > i) {
            goldilocks_t temp = a[i];
            a[i] = a[j];
            a[j] = temp;
        }
    }
    
    // Cooley-Tukey FFT for the polynomial
    for (uint32_t s = 1; s <= m; s++) {
        uint32_t m_s = 1 << s;
        goldilocks_t w_m = goldilocks_pow(omega, 1 << (total_order - s));
        
        for (uint32_t k = 0; k < poly_size; k += m_s) {
            goldilocks_t w = 1;
            for (uint32_t j = 0; j < m_s/2; j++) {
                goldilocks_t t = goldilocks_mul_naive(w, a[k + j + m_s/2]);
                goldilocks_t u = a[k + j];
                a[k + j] = goldilocks_add_naive(u, t);
                a[k + j + m_s/2] = goldilocks_sub_naive(u, t);
                w = goldilocks_mul_naive(w, w_m);
            }
        }
    }
    
    // Extend the evaluations to 2^(m+n) points
    for (uint32_t s = m + 1; s <= total_order; s++) {
        uint32_t m_s = 1 << s;
        uint32_t prev_m_s = 1 << (s - 1);
        goldilocks_t w_m = goldilocks_pow(omega, 1 << (total_order - s));
        
        for (uint32_t k = 0; k < eval_size; k += m_s) {
            goldilocks_t w = 1;
            for (uint32_t j = 0; j < prev_m_s; j++) {
                goldilocks_t t = goldilocks_mul_naive(w, a[k + j + prev_m_s]);
                goldilocks_t u = a[k + j];
                a[k + j] = goldilocks_add_naive(u, t);
                a[k + j + prev_m_s] = goldilocks_sub_naive(u, t);
                w = goldilocks_mul_naive(w, w_m);
            }
        }
    }
    
    // For inverse NTT, multiply by 1/N
    if (is_inverse) {
        goldilocks_t n_inv = goldilocks_inv_naive(poly_size);
        for (uint32_t i = 0; i < poly_size; i++) {
            a[i] = goldilocks_mul_naive(a[i], n_inv);
        }
    }
}

// Returns the primitive root of unity of order size
goldilocks_t goldilocks_root_naive(size_t size) {
    // For Goldilocks field, primitive root is 7
    goldilocks_t root = 7;
    
    // Calculate the maximum order (p-1)/2
    size_t max_order = (GOLDILOCKS_P - 1) / 2;
    
    // Calculate the required order
    size_t required_order = max_order;
    while (required_order > size) {
        required_order >>= 1;
        root = goldilocks_mul_naive(root, root);
    }
    
    return root;
}

#endif // GOLDILOCKS_H 
