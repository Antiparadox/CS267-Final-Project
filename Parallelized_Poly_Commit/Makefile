CXX = CC
CXXFLAGS = -O3 -Wall -Wextra
OPENMP_FLAGS = -fopenmp
MPI_FLAGS = 

# Source files
NAIVE_SOURCES = polynomial_commitment_naive.cpp FFT.cpp
OMP_SOURCES = polynomial_commitment_omp.cpp FFT_omp.cpp
MPI_SOURCES = polynomial_commitment_mpi.cpp FFT_MPI.cpp

# Object files
NAIVE_OBJS = $(NAIVE_SOURCES:.cpp=.o)
OMP_OBJS = $(OMP_SOURCES:.cpp=.o)
MPI_OBJS = $(MPI_SOURCES:.cpp=.o)

# Executables
NAIVE_EXEC = polynomial_commitment_naive
OMP_EXEC = polynomial_commitment_omp
MPI_EXEC = polynomial_commitment_mpi
GEN_DATASET = generate_dataset

# Libraries
LIBS = -lcrypto

# Default target
all: $(NAIVE_EXEC) $(OMP_EXEC) $(MPI_EXEC) $(GEN_DATASET) setup_runner
	@echo "\nPolynomial Commitment Interactive Console"
	@echo "======================================="
	@echo "Quick Start Tips:"
	@echo "1. Run the interactive console:"
	@echo "   ./launch_polynomial_commitment.sh"
	@echo
	@echo "2. Start with test mode (option 1) to verify correctness"
	@echo "3. For OpenMP, try different thread counts"
	@echo "4. For MPI, first run: salloc -N 2 -C cpu -q interactive -t 00:30:00"
	@echo
	@echo "For detailed usage guide, run: make help"
	@echo "=======================================\n"

# Help target
help:
	@echo "Polynomial Commitment Interactive Console Usage Guide"
	@echo "=================================================="
	@echo
	@echo "1. Building the Project"
	@echo "   make all        - Build all implementations"
	@echo "   make naive      - Build only naive implementation"
	@echo "   make omp        - Build only OpenMP implementation"
	@echo "   make mpi        - Build only MPI implementation"
	@echo "   make clean      - Remove all built files"
	@echo
	@echo "2. Running the Interactive Console"
	@echo "   ./launch_polynomial_commitment.sh"
	@echo
	@echo "3. Implementation Options"
	@echo "   1) Naive (sequential, single-threaded)"
	@echo "   2) OpenMP (shared-memory parallel)"
	@echo "   3) MPI (distributed-memory parallel)"
	@echo
	@echo "4. Running Modes"
	@echo "   a) Test Mode"
	@echo "      - Uses test dataset with n=0"
	@echo "      - Verifies NTT correctness"
	@echo "      - Good for debugging"
	@echo
	@echo "   b) Performance Mode"
	@echo "      - Choose dataset size:"
	@echo "        * Small  (2^10 = 1,024 elements)"
	@echo "        * Medium (2^15 = 32,768 elements)"
	@echo "        * Large  (2^22 = 4,194,304 elements)"
	@echo "      - Configure evaluation domain size (n)"
	@echo "      - Set parallelization parameters"
	@echo
	@echo "5. MPI Requirements"
	@echo "   Before running MPI implementation:"
	@echo "   salloc -N 2 -C cpu -q interactive -t 00:30:00"
	@echo "   (Replace '2' with desired number of nodes)"
	@echo
	@echo "6. Tips"
	@echo "   - Start with test mode to verify correctness"
	@echo "   - Use small dataset for quick testing"
	@echo "   - Increase dataset size gradually"
	@echo "   - For OpenMP, try different thread counts"
	@echo "   - For MPI, ensure proper node allocation"
	@echo
	@echo "For more details, see the script comments and documentation"

# Naive implementation
$(NAIVE_EXEC): $(NAIVE_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# OpenMP implementation
$(OMP_EXEC): $(OMP_OBJS)
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) -o $@ $^ $(LIBS)

# MPI implementation
$(MPI_EXEC): $(MPI_OBJS) FFT.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Generate dataset
$(GEN_DATASET): generate_dataset.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

# Make run script executable
setup_runner:
	chmod +x launch_polynomial_commitment.sh

# Compilation rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

# OpenMP specific compilation
FFT_omp.o: FFT_omp.cpp
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) -c $<

polynomial_commitment_omp.o: polynomial_commitment_omp.cpp
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) -c $<

# MPI specific compilation
FFT_MPI.o: FFT_MPI.cpp
	$(CXX) $(CXXFLAGS) -c $<

polynomial_commitment_mpi.o: polynomial_commitment_mpi.cpp
	$(CXX) $(CXXFLAGS) -c $<

# Individual targets
naive: $(NAIVE_EXEC)
omp: $(OMP_EXEC)
mpi: $(MPI_EXEC)
dataset: $(GEN_DATASET)

# Clean target
clean:
	rm -f $(NAIVE_OBJS) $(OMP_OBJS) $(MPI_OBJS) $(NAIVE_EXEC) $(OMP_EXEC) $(MPI_EXEC) $(GEN_DATASET) *.txt

.PHONY: all naive omp mpi dataset clean setup_runner help 
