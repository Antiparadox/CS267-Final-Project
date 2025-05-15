#!/bin/bash

# Colors for better readability
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Function to print colored text
print_color() {
    color=$1
    text=$2
    echo -e "${color}${text}${NC}"
}

# Function to print warning messages
print_warning() {
    print_color $RED "ERROR: $1"
    print_color $RED "This computation size is not allowed."
    print_color $YELLOW "Please choose a smaller n value or dataset size."
    print_color $YELLOW "Returning to main menu..."
    sleep 2  # Give user time to read the message
    return 1
}

# Function to check if m+n is too large
check_computation_size() {
    local m=$1
    local n=$2
    if ((m + n > 31)); then
        print_warning "Computation size 2^($m+$n) = 2^$((m+n)) exceeds maximum allowed size of 2^31!"
        return 1
    fi
    return 0
}

# Function to print a separator line
print_separator() {
    echo -e "${CYAN}==================================================${NC}"
}

# Function to print a welcome message
print_welcome() {
    clear
    print_separator
    print_color $BOLD "Welcome to the Polynomial Commitment Scheme Implementation!"
    print_separator
    echo
    print_color $GREEN "This tool implements a complete polynomial commitment scheme:"
    print_color $GREEN "1. NTT: Evaluates polynomial over a smooth subgroup of size 2^m"
    print_color $GREEN "   (where 2^m is the number of polynomial coefficients)"
    print_color $GREEN "2. Inverse NTT: Recovers coefficients over a coset of size 2^(m+n)"
    print_color $GREEN "   (where n is the blowup factor for improved codeword distance)"
    print_color $GREEN "3. Merkle Tree: Creates commitment over NTT evaluations"
    echo
    print_color $YELLOW "Implementation Options:"
    print_color $YELLOW "- Naive: Sequential implementation"
    print_color $YELLOW "- OpenMP: Shared-memory parallel version"
    print_color $YELLOW "- MPI: Distributed-memory parallel version"
    echo
    print_color $YELLOW "Quick Tips:"
    print_color $YELLOW "- Start with test mode to verify NTT correctness"
    print_color $YELLOW "- Use small datasets for initial testing"
    print_color $YELLOW "- For MPI, ensure you've allocated nodes with salloc"
    echo
    print_separator
    echo
}

# Print welcome message
print_welcome

# Function to run naive implementation
run_naive() {
    m=$1
    n=$2
    test_mode=$3
    print_color $BLUE "Running naive implementation..."
    print_color $BLUE "Using m=${m}, n=${n}, test_mode=${test_mode}"
    if [ "$test_mode" = "test" ]; then
        export DATASET_FILE="test_dataset.txt"
    fi
    ./polynomial_commitment_naive $n $test_mode
}

# Function to run OpenMP implementation
run_omp() {
    m=$1
    n=$2
    threads=$3
    test_mode=$4
    print_color $BLUE "Running OpenMP implementation with ${threads} threads..."
    print_color $BLUE "Using m=${m}, n=${n}, test_mode=${test_mode}"
    if [ "$test_mode" = "test" ]; then
        export DATASET_FILE="test_dataset.txt"
        threads=8  # Force 8 threads in test mode
    fi
    OMP_NUM_THREADS=$threads ./polynomial_commitment_omp $n $test_mode $threads
}

# Function to check MPI environment
check_mpi_env() {
    if [ -z "$SLURM_JOB_ID" ]; then
        print_color $RED "Error: MPI environment not detected!"
        print_color $YELLOW "To run MPI implementation, you need to:"
        print_color $YELLOW "1. Request compute nodes using:"
        print_color $YELLOW "   salloc -N 2 -C cpu -q interactive -t 00:30:00"
        print_color $YELLOW "2. Then run this script again"
        print_color $YELLOW "Note: Replace '2' with desired number of nodes"
        return 1
    fi
    return 0
}

# Function to run MPI implementation
run_mpi() {
    m=$1
    n=$2
    processes=$3
    test_mode=$4
    
    # Check MPI environment
    if ! check_mpi_env; then
        return 1
    fi
    
    print_color $BLUE "Running MPI implementation with ${processes} processes..."
    print_color $BLUE "Using m=${m}, n=${n}, test_mode=${test_mode}"
    if [ "$test_mode" = "test" ]; then
        export DATASET_FILE="test_dataset.txt"
        processes=8  # Force 8 processes in test mode
    fi
    srun -n $processes ./polynomial_commitment_mpi $n $test_mode
}

# Function to generate test dataset
generate_test_dataset() {
    size=$1
    output_file=$2
    print_color $GREEN "Generating ${output_file} with ${size} elements..."
    python3 -c "
with open('${output_file}', 'w') as f:
    for i in range(1, ${size} + 1):
        f.write(f'{i}\\n')
"
}

# Generate datasets
print_color $GREEN "Generating datasets..."
generate_test_dataset 1024 "test_dataset.txt"
generate_test_dataset 1024 "small_dataset.txt"
generate_test_dataset 32768 "medium_dataset.txt"
generate_test_dataset 4194304 "large_dataset.txt"

# Main loop
while true; do
    # Main menu
    echo
    print_separator
    print_color $BOLD "Select implementation:"
    print_color $CYAN "1) Run naive implementation (sequential, single-threaded)"
    print_color $CYAN "2) Run OpenMP implementation (shared-memory parallel)"
    print_color $CYAN "3) Run MPI implementation (distributed-memory parallel)"
    print_color $CYAN "4) Exit"
    print_separator
    read -p "Enter your choice (1-4): " choice

    case $choice in
        1|2|3) 
            # Common menu for all implementations
            echo
            print_separator
            print_color $BOLD "Select mode:"
            print_color $CYAN "1) Test mode - Verify correctness"
            print_color $CYAN "   - Uses test dataset with n=0 (no blowup)"
            print_color $CYAN "   - Performs NTT over subgroup of size 2^m"
            print_color $CYAN "   - Performs inverse NTT over the same subgroup"
            print_color $CYAN "   - Verifies recovered coefficients match original"
            print_color $CYAN "   - Useful for debugging and validation"
            echo
            print_color $CYAN "2) Performance mode - Run with custom parameters"
            print_color $CYAN "   - Select dataset size (small/medium/large)"
            print_color $CYAN "   - Choose evaluation domain size (n)"
            print_color $CYAN "   - Configure parallelization parameters"
            print_color $CYAN "   - Measure execution time"
            print_separator
            read -p "Enter choice (1-2): " mode_choice

            case $mode_choice in
                1)
                    # Test mode - always use test_dataset.txt and 8 threads/processes
                    case $choice in
                        1) run_naive 10 0 "test" ;;
                        2) run_omp 10 0 8 "test" ;;
                        3) 
                            if ! check_mpi_env; then
                                continue
                            fi
                            run_mpi 10 0 8 "test" 
                            ;;
                    esac
                    ;;
                2)
                    # No test mode
                    echo
                    print_separator
                    print_color $BOLD "Select dataset size:"
                    print_color $CYAN "1) Small dataset (2^10 = 1,024 elements)"
                    print_color $CYAN "   - Good for quick testing and development"
                    print_color $CYAN "   - Suitable for all implementations"
                    echo
                    print_color $CYAN "2) Medium dataset (2^15 = 32,768 elements)"
                    print_color $CYAN "   - Moderate computational load"
                    print_color $CYAN "   - Good for benchmarking parallel implementations"
                    echo
                    print_color $CYAN "3) Large dataset (2^22 = 4,194,304 elements)"
                    print_color $CYAN "   - Heavy computational load"
                    print_color $CYAN "   - Best for measuring parallel performance"
                    print_color $CYAN "   - May require significant memory"
                    print_separator
                    read -p "Enter choice (1-3): " dataset_choice

                    case $dataset_choice in
                        1) 
                            m=10
                            export DATASET_FILE="small_dataset.txt"
                            ;;
                        2) 
                            m=15
                            export DATASET_FILE="medium_dataset.txt"
                            ;;
                        3) 
                            m=22
                            export DATASET_FILE="large_dataset.txt"
                            ;;
                        *) 
                            print_color $RED "Invalid choice"
                            continue
                            ;;
                    esac

                    print_separator
                    print_color $BOLD "Enter blowup factor (n):"
                    print_color $CYAN "This determines the size of the evaluation domain"
                    print_color $CYAN "Total evaluation domain size will be 2^(m+n)"
                    print_color $CYAN "where m=$m (current dataset size)"
                    print_color $CYAN "Larger n provides better zero-knowledge properties"
                    print_color $CYAN "but requires more computation"
                    print_separator
                    read -p "Enter n value: " n
                    
                    # Check if computation size is too large
                    if ! check_computation_size $m $n; then
                        continue  # This will return to the main menu
                    fi

                    case $choice in
                        1)
                            run_naive $m $n "no_test"
                            ;;
                        2)
                            read -p "Enter number of threads: " threads
                            run_omp $m $n $threads "no_test"
                            ;;
                        3)
                            if ! check_mpi_env; then
                                continue
                            fi
                            read -p "Enter number of processes: " processes
                            run_mpi $m $n $processes "no_test"
                            ;;
                    esac
                    ;;
                *)
                    print_color $RED "Invalid choice"
                    ;;
            esac
            ;;
        4)
            print_separator
            print_color $GREEN "Thank you for using the Polynomial Commitment Interactive Console!"
            print_color $GREEN "Goodbye!"
            print_separator
            exit 0
            ;;
        *)
            print_color $RED "Invalid choice"
            ;;
    esac

    # Wait for user to press Enter before showing menu again
    echo
    print_separator
    print_color $GREEN "Execution completed. Press Enter to return to main menu..."
    print_separator
    read
done 
