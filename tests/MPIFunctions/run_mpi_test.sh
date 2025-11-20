#!/bin/bash
# Wrapper script to run MPI test with 2 ranks

# Don't exit on error - we want to try all commands
set +e

echo "=== MPI Test Runner ==="
echo "Current directory: $(pwd)"
echo "Date: $(date)"
echo ""

# Check if test executable exists
if [ ! -f ./test_mpi_parallel ]; then
    echo "ERROR: test_mpi_parallel executable not found"
    echo "Contents of current directory:"
    ls -la
    exit 1
fi

echo "Found test executable: ./test_mpi_parallel"
echo ""

# List of MPI commands to try
declare -a MPI_COMMANDS=(
    "mpirun"
    "mpiexec"
    "/usr/bin/mpirun"
    "/usr/bin/mpiexec"
    "/usr/local/bin/mpirun"
    "/usr/local/bin/mpiexec"
)

declare -a FLAGS_TO_TRY=(
    ""
    "--oversubscribe"
    "--allow-run-as-root"
    "--oversubscribe --allow-run-as-root"
)

echo "=== Environment Information ==="
echo "PATH: $PATH"
echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH:-<not set>}"
echo "USER: $USER"
echo ""

echo "=== Searching for MPI commands ==="
for cmd in "${MPI_COMMANDS[@]}"; do
    if command -v "$cmd" &> /dev/null || [ -x "$cmd" ]; then
        echo "✓ Found: $cmd"
        full_path=$(command -v "$cmd" 2>/dev/null || echo "$cmd")
        echo "  Path: $full_path"

        # Try to get version
        version_output=$("$cmd" --version 2>&1 || echo "Unable to get version")
        echo "  Version info: $(echo "$version_output" | head -n 1)"
    else
        echo "✗ Not found: $cmd"
    fi
done
echo ""

# Try each MPI command with each flag combination
echo "=== Attempting to run MPI test ==="
test_passed=0

for cmd in "${MPI_COMMANDS[@]}"; do
    if command -v "$cmd" &> /dev/null || [ -x "$cmd" ]; then
        for flags in "${FLAGS_TO_TRY[@]}"; do
            echo "----------------------------------------"
            echo "Trying: $cmd $flags -n 2 ./test_mpi_parallel"

            output=$($cmd $flags -n 2 ./test_mpi_parallel 2>&1)
            exit_code=$?

            echo "Output:"
            echo "$output"
            echo "Exit code: $exit_code"

            if [ $exit_code -eq 0 ]; then
                echo "SUCCESS: Test passed with '$cmd $flags'"
                test_passed=1
                break 2
            else
                echo "FAILED: Test failed with exit code $exit_code"
            fi
            echo ""
        done
    fi
done

if [ $test_passed -eq 1 ]; then
    echo "=== TEST PASSED ==="
    exit 0
else
    echo "=== ALL ATTEMPTS FAILED ==="
    echo "None of the MPI command/flag combinations succeeded"
    exit 1
fi
