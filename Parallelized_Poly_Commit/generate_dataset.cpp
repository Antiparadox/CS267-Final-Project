#include <iostream>
#include <fstream>
#include <random>
#include <cstdint>
#include <string>

void generate_test_dataset(const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open output file " << filename << std::endl;
        return;
    }

    // Generate numbers from 1 to 2^10
    const size_t num_elements = 1ULL << 10;
    for (size_t i = 1; i <= num_elements; i++) {
        outfile << i << std::endl;
    }

    outfile.close();
    std::cout << "Generated test dataset with " << num_elements << " sequential numbers in " << filename << std::endl;
}

void generate_random_dataset(const std::string& filename, size_t power, const std::string& size_desc) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(1, 1ULL << 15);

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open output file " << filename << std::endl;
        return;
    }

    const size_t num_elements = 1ULL << power;
    for (size_t i = 0; i < num_elements; i++) {
        outfile << dis(gen) << std::endl;
    }

    outfile.close();
    std::cout << "Generated " << size_desc << " dataset with " << num_elements << " random numbers in " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <dataset_type>" << std::endl;
        std::cerr << "Dataset types:" << std::endl;
        std::cerr << "  test   - Sequential numbers 1 to 2^10" << std::endl;
        std::cerr << "  small  - 2^10 random numbers" << std::endl;
        std::cerr << "  medium - 2^15 random numbers" << std::endl;
        std::cerr << "  large  - 2^22 random numbers" << std::endl;
        return 1;
    }

    std::string type = argv[1];
    
    if (type == "test") {
        generate_test_dataset("test_dataset.txt");
    }
    else if (type == "small") {
        generate_random_dataset("small_dataset.txt", 10, "small");
    }
    else if (type == "medium") {
        generate_random_dataset("medium_dataset.txt", 15, "medium");
    }
    else if (type == "large") {
        generate_random_dataset("large_dataset.txt", 22, "large");
    }
    else {
        std::cerr << "Invalid dataset type. Use: test, small, medium, or large" << std::endl;
        return 1;
    }

    return 0;
} 