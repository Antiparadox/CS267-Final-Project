# CS267 Final Project: Parallelized Polynomial Commitment Scheme

This project implements a polynomial commitment scheme with three different implementations:
1. Naive (sequential) implementation
2. OpenMP (shared-memory parallel) implementation
3. MPI (distributed-memory parallel) implementation

## Prerequisites

- C++17 compatible compiler (g++ 7.0 or later)
- OpenMP (usually comes with g++)
- MPI (for distributed implementation)
- OpenSSL development libraries (for SHA256 hashing)

### Installing Dependencies

#### Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install g++ make libopenmpi-dev libssl-dev
```

#### macOS (using Homebrew):
```bash
brew install gcc make open-mpi openssl
```

## Project Structure

```
Parallelized_Poly_Commit/
├── FFT.cpp                 # FFT implementation
├── FFT_MPI.cpp            # MPI version of FFT
├── FFT_omp.cpp            # OpenMP version of FFT
├── Makefile               # Main Makefile
├── ntt.h                  # Common NTT definitions
├── ntt_mpi.h             # MPI-specific NTT definitions
├── ntt_naive.h           # Naive NTT implementation
├── ntt_omp.h             # OpenMP NTT definitions
├── polynomial_commitment_mpi.cpp    # MPI implementation
├── polynomial_commitment_naive.cpp  # Naive implementation
├── polynomial_commitment_omp.cpp    # OpenMP implementation
└── launch_polynomial_commitment.sh  # Interactive launcher script
```

## Compilation

### Compiling All Implementations

```bash
cd Parallelized_Poly_Commit
make all
```

This will create three executables:
- `polynomial_commitment_naive`
- `polynomial_commitment_omp`
- `polynomial_commitment_mpi`

### Compiling Individual Implementations

```bash
# Naive implementation
make naive

# OpenMP implementation
make omp

# MPI implementation
make mpi
```

## Usage

### Interactive Launcher

The easiest way to run the code is using the interactive launcher script:

```bash
./launch_polynomial_commitment.sh
```

This script provides a menu-driven interface to:
1. Choose the implementation (naive, OpenMP, or MPI)
2. Select the mode (test or performance)
3. Configure parameters (dataset size, number of threads/processes)

### Manual Execution

#### Naive Implementation
```bash
./polynomial_commitment_naive <n> <test_mode>
# Example: ./polynomial_commitment_naive 2 test
```

#### OpenMP Implementation
```bash
OMP_NUM_THREADS=<num_threads> ./polynomial_commitment_omp <n> <test_mode> <num_threads>
# Example: OMP_NUM_THREADS=4 ./polynomial_commitment_omp 2 test 4
```

#### MPI Implementation
```bash
# First, allocate compute nodes (if using a cluster)
salloc -N <num_nodes> -C cpu -q interactive -t 00:30:00

# Then run the program
srun -n <num_processes> ./polynomial_commitment_mpi <n> <test_mode>
# Example: srun -n 4 ./polynomial_commitment_mpi 2 test
```

### Parameters

- `n`: Blowup factor (0 for test mode)
- `test_mode`: Either "test" or "no_test"
  - "test": Uses a small dataset and verifies correctness
  - "no_test": Uses the specified dataset size for performance measurement
- `num_threads`: Number of OpenMP threads (for OpenMP implementation)
- `num_processes`: Number of MPI processes (for MPI implementation)

## Test Datasets

The project includes several test datasets:
- `test_dataset.txt`: Small dataset (1024 elements) for testing
- `small_dataset.txt`: 1024 elements
- `medium_dataset.txt`: 32768 elements
- `large_dataset.txt`: 4194304 elements

You can specify a custom dataset using the `DATASET_FILE` environment variable:
```bash
export DATASET_FILE=path/to/your/dataset.txt
```

## Performance Considerations

1. For OpenMP:
   - Optimal thread count is typically equal to the number of CPU cores
   - Use `OMP_NUM_THREADS` to control thread count

2. For MPI:
   - Number of processes should be a power of 2
   - Each process should have enough memory for its portion of the data
   - Consider network bandwidth when choosing number of nodes

## Troubleshooting

1. If you get compilation errors about OpenSSL:
   - Make sure libssl-dev is installed
   - You might need to specify the OpenSSL include path:
     ```bash
     make CXXFLAGS="-I/usr/local/opt/openssl/include" LDFLAGS="-L/usr/local/opt/openssl/lib"
     ```

2. If MPI programs fail to start:
   - Make sure you've allocated compute nodes using `salloc`
   - Check that the number of processes is a power of 2
   - Verify MPI installation with `mpirun --version`

3. If you get "permission denied" errors:
   - Make the launcher script executable:
     ```bash
     chmod +x launch_polynomial_commitment.sh
     ```

## License

This project is part of the CS267 course at UC Berkeley.

## Authors

- Your Name
- Your Partner's Name (if applicable) 