#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

EXIT_CODE=0

APPLY_FIXES=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --apply) APPLY_FIXES=$1;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

echo "🚀 Starting Full Project Linting..."
echo "Location: $SCRIPT_DIR"
echo "----------------------------------------"

# Function to run sibling scripts
run_sub_script() {
    local script_name=$1
    local label=$2
    local full_path="$SCRIPT_DIR/$script_name"

    echo "🏃 Running $label..."
    
    if [ ! -f "$full_path" ]; then
        echo "⚠️  Error: $script_name not found in $SCRIPT_DIR"
        EXIT_CODE=1
        return
    fi

    # Run the script and capture its exit code
    bash "$full_path" ${APPLY_FIXES:+"$APPLY_FIXES"}
    local status=$?

    if [ $status -eq 0 ]; then
        echo "✅ $label passed."
    else
        echo "❌ $label failed (Exit Code: $status)."
        EXIT_CODE=1
    fi
    echo "----------------------------------------"
}

# Run all three scripts sequentially
# Note: Even if one fails, the next one will still execute
run_sub_script "cpp-formatting.sh" "C++ Formatting"
run_sub_script "python-formatting.sh" "Python Formatting"
run_sub_script "cmake-formatting.sh" "CMake Formatting"

# Final Summary for CI/CD
if [ $EXIT_CODE -eq 0 ]; then
    echo "🎉 ALL CHECKS PASSED"
    exit 0
else
    echo "⛔ SOME CHECKS FAILED"
    echo "Please review the logs above and run with --apply locally."
    exit 1
fi