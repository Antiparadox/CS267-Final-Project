#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <utility>
#include <openssl/sha.h>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <omp.h>
#include "ntt_omp.h"
#include "ntt_naive.h"  // Include naive implementation for comparison

// Function to compute SHA256 hash
std::string sha256(const std::string& str) {
    unsigned char hash[SHA256_DIGEST_LENGTH];
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    SHA256_Update(&sha256, str.c_str(), str.size());
    SHA256_Final(hash, &sha256);
    
    std::stringstream ss;
    for(int i = 0; i < SHA256_DIGEST_LENGTH; i++) {
        ss << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(hash[i]);
    }
    return ss.str();
}

// Function to build Merkle tree (naive implementation)
std::string build_merkle_tree_naive(const std::vector<std::string>& leaves, std::vector<double>& layer_times) {
    if (leaves.empty()) return "";
    if (leaves.size() == 1) return leaves[0];
    
    std::vector<std::string> current = leaves;
    layer_times.clear();  // Clear any previous timing data
    
    while (current.size() > 1) {
        auto layer_start = std::chrono::high_resolution_clock::now();
        
        std::vector<std::string> next;
        for (size_t i = 0; i < current.size(); i += 2) {
            if (i + 1 < current.size()) {
                next.push_back(sha256(current[i] + current[i + 1]));
            } else {
                next.push_back(sha256(current[i] + current[i]));
            }
        }
        current = next;
        
        auto layer_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> layer_time = layer_end - layer_start;
        layer_times.push_back(layer_time.count());
    }
    return current[0];
}

// Function to build Merkle tree (OpenMP implementation)
std::string build_merkle_tree_parallel(const std::vector<std::string>& leaves, std::vector<double>& layer_times) {
    if (leaves.empty()) return "";
    if (leaves.size() == 1) return leaves[0];
    
    std::vector<std::string> current = leaves;
    layer_times.clear();  // Clear any previous timing data
    
    while (current.size() > 1) {
        auto layer_start = std::chrono::high_resolution_clock::now();
        
        std::vector<std::string> next(current.size() / 2 + (current.size() % 2));
        
        #pragma omp parallel for
        for (size_t i = 0; i < current.size(); i += 2) {
            if (i + 1 < current.size()) {
                next[i/2] = sha256(current[i] + current[i + 1]);
            } else {
                next[i/2] = sha256(current[i] + current[i]);
            }
        }
        current = next;
        
        auto layer_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> layer_time = layer_end - layer_start;
        layer_times.push_back(layer_time.count());
    }
    return current[0];
}

// Function to find the next power of 2 and its log2
std::pair<int, uint32_t> getPaddedSizeAndLog2(int n) {
    int power = 1;
    uint32_t log2 = 0;
    while (power < n) {
        power *= 2;
        log2++;
    }
    return std::make_pair(power, log2);
}

// Function to read vector from file
std::vector<goldilocks_t> readVectorFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }
    
    std::vector<goldilocks_t> vec;
    goldilocks_t val;
    while (file >> val) {
        vec.push_back(val);
    }
    return vec;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <n> <test/no_test> <num_threads>" << std::endl;
        return 1;
    }

    int n = std::stoi(argv[1]);
    std::string test_mode = argv[2];
    int num_threads = std::stoi(argv[3]);
    
    if (test_mode != "test" && test_mode != "no_test") {
        std::cerr << "Second argument must be 'test' or 'no_test'" << std::endl;
        return 1;
    }

    // Set number of OpenMP threads
    omp_set_num_threads(num_threads);

    // If test mode, set n to 0
    if (test_mode == "test") {
        n = 0;
    }

    // Start timing initialization
    auto start_init = std::chrono::high_resolution_clock::now();

    // Read input polynomial from file
    std::string input_file = "test_dataset.txt";  // Default to test dataset
    if (const char* env_dataset = std::getenv("DATASET_FILE")) {
        input_file = env_dataset;
    }
    std::vector<goldilocks_t> polynomial = readVectorFromFile(input_file);
    
    // Calculate m based on input size
    int m = 0;
    size_t size = polynomial.size();
    while (size > 1) {
        size >>= 1;
        m++;
    }
    
    std::cout << "Input size: " << polynomial.size() << ", calculated m: " << m << std::endl;
    std::cout << "Using n: " << n << std::endl;
    std::cout << "Number of OpenMP threads: " << num_threads << std::endl;

    // Validate polynomial size
    if (polynomial.size() != (1ULL << m)) {
        std::cerr << "Error: Input polynomial size must be a power of 2" << std::endl;
        return 1;
    }

    // Pad or truncate to the required size
    if (polynomial.size() < (1ULL << m)) {
        polynomial.resize((1ULL << m), 0);
    } else if (polynomial.size() > (1ULL << m)) {
        polynomial.resize((1ULL << m));
    }

    // Resize to evaluation domain size and zero-pad
    polynomial.resize((1ULL << (m + n)), 0);

    // End timing initialization
    auto end_init = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> init_time = end_init - start_init;

    if (test_mode == "test") {
        // Store first 5 coefficients for verification
        std::vector<goldilocks_t> first_five_coeffs;
        for (int i = 0; i < 5; i++) {
            first_five_coeffs.push_back(polynomial[i]);
        }

        // Print first 5 coefficients
        std::cout << "First 5 input polynomial coefficients:" << std::endl;
        for (int i = 0; i < 5; i++) {
            std::cout << first_five_coeffs[i] << " ";
        }
        std::cout << std::endl;

        // Perform NTT using OpenMP
        ntt_omp(polynomial.data(), m, n, false);

        // Print first 5 evaluations
        std::cout << "\nFirst 5 evaluations:" << std::endl;
        for (int i = 0; i < 5; i++) {
            std::cout << polynomial[i] << " ";
        }
        std::cout << std::endl;

        // Perform inverse NTT using OpenMP
        ntt_omp(polynomial.data(), m, n, true);

        // Print first 5 recovered coefficients
        std::cout << "\nFirst 5 recovered coefficients:" << std::endl;
        for (int i = 0; i < 5; i++) {
            std::cout << polynomial[i] << " ";
        }
        std::cout << std::endl;

        // Verify if recovered coefficients match input
        bool ntt_match = true;
        for (int i = 0; i < 5; i++) {
            if (polynomial[i] != first_five_coeffs[i]) {
                ntt_match = false;
                break;
            }
        }
        if (ntt_match) {
            std::cout << "\nNTT Correctness test passed: Recovered coefficients match input coefficients" << std::endl;
        } else {
            std::cout << "\nNTT Correctness test failed: Recovered coefficients do not match input coefficients" << std::endl;
        }

        // Build Merkle tree from evaluations using both implementations
        std::vector<std::string> leaves;
        for (uint32_t i = 0; i < (1ULL << (m + n)); i++) {
            leaves.push_back(std::to_string(polynomial[i]));
        }

        std::cout << "\nComputing Merkle root using naive implementation..." << std::endl;
        std::vector<double> naive_layer_times;
        std::string naive_root = build_merkle_tree_naive(leaves, naive_layer_times);
        std::cout << "Merkle root from naive implementation: " << naive_root << std::endl;
        
        // Print naive implementation layer times
        std::cout << "\nNaive Implementation Layer Times:" << std::endl;
        double total_naive_layer_time = 0.0;
        for (size_t i = 0; i < naive_layer_times.size(); i++) {
            std::cout << "Layer " << i + 1 << " time: " << naive_layer_times[i] << " seconds" << std::endl;
            total_naive_layer_time += naive_layer_times[i];
        }
        std::cout << "Total layer computation time: " << total_naive_layer_time << " seconds" << std::endl;

        std::cout << "\nComputing Merkle root using OpenMP implementation..." << std::endl;
        std::vector<double> omp_layer_times;
        std::string omp_root = build_merkle_tree_parallel(leaves, omp_layer_times);
        std::cout << "Merkle root from OpenMP implementation: " << omp_root << std::endl;
        
        // Print OpenMP implementation layer times
        std::cout << "\nOpenMP Implementation Layer Times:" << std::endl;
        double total_omp_layer_time = 0.0;
        for (size_t i = 0; i < omp_layer_times.size(); i++) {
            std::cout << "Layer " << i + 1 << " time: " << omp_layer_times[i] << " seconds" << std::endl;
            total_omp_layer_time += omp_layer_times[i];
        }
        std::cout << "Total layer computation time: " << total_omp_layer_time << " seconds" << std::endl;

        // Verify if Merkle roots match
        if (naive_root == omp_root) {
            std::cout << "\nMerkle tree Correctness test passed: OpenMP implementation matches naive implementation" << std::endl;
        } else {
            std::cout << "\nMerkle tree Correctness test failed: OpenMP implementation does not match naive implementation" << std::endl;
        }
    } else {
        // Performance mode: measure and output timing
        // Print initialization time
        std::cout << "\nInitialization time: " << init_time.count() << " seconds" << std::endl;

        // Perform NTT
        auto start_ntt = std::chrono::high_resolution_clock::now();
        ntt_omp(polynomial.data(), m, n, false);
        auto end_ntt = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> ntt_time = end_ntt - start_ntt;

        // Perform inverse NTT
        auto start_intt = std::chrono::high_resolution_clock::now();
        ntt_omp(polynomial.data(), m, n, true);
        auto end_intt = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> intt_time = end_intt - start_intt;

        // Build Merkle tree from evaluations
        std::vector<std::string> leaves;
        for (uint32_t i = 0; i < (1ULL << (m + n)); i++) {
            leaves.push_back(std::to_string(polynomial[i]));
        }

        auto start_merkle = std::chrono::high_resolution_clock::now();
        std::vector<double> layer_times;
        std::string root = build_merkle_tree_parallel(leaves, layer_times);
        auto end_merkle = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> merkle_time = end_merkle - start_merkle;

        // Print timing results
        std::cout << "\nTiming Results:" << std::endl;
        std::cout << "Initialization time: " << init_time.count() << " seconds" << std::endl;
        std::cout << "NTT computation time: " << ntt_time.count() << " seconds" << std::endl;
        std::cout << "Inverse NTT computation time: " << intt_time.count() << " seconds" << std::endl;
        
        // Print detailed Merkle tree timing
        std::cout << "\nMerkle Tree Timing Details:" << std::endl;
        double total_layer_time = 0.0;
        for (size_t i = 0; i < layer_times.size(); i++) {
            std::cout << "Layer " << i + 1 << " time: " << layer_times[i] << " seconds" << std::endl;
            total_layer_time += layer_times[i];
        }
        std::cout << "Total Merkle tree computation time: " << merkle_time.count() << " seconds" << std::endl;
        std::cout << "Total layer computation time: " << total_layer_time << " seconds" << std::endl;
        std::cout << "Overhead time: " << (merkle_time.count() - total_layer_time) << " seconds" << std::endl;
        
        std::cout << "Total execution time: " << (init_time.count() + ntt_time.count() + intt_time.count() + merkle_time.count()) << " seconds" << std::endl;
        std::cout << "Merkle root: " << root << std::endl;
    }

    return 0;
} 
