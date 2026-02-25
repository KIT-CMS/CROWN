#!/bin/bash

# Configuration
APPLY_FIXES=false

# Parse CLI arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --apply) APPLY_FIXES=true ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

echo "🔍 Finding CMake files ..."

# Find tracked CMake files:
# 1. Any file named CMakeLists.txt
# 2. Any file ending in .cmake or .cmake.in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILES=$(git ls-files | grep -E '(CMakeLists\.txt|\.cmake(\.in)?)$' | grep -vEf "$SCRIPT_DIR/format_ignore.txt")

if [[ -z "$FILES" ]]; then
    echo "No CMake files found."
    exit 0
fi

if [ "$APPLY_FIXES" = true ]; then
    echo "🛠  Running cmake-format to fix formatting..."
    # -i: in-place
    echo "$FILES" | xargs cmake-format -i
    echo "✅ CMake formatting complete."
else
    echo "🧪 Checking CMake formatting (Dry Run)..."
    
    # --check: return exit code 1 if files need formatting
    if echo "$FILES" | xargs cmake-format --check; then
        echo "✨ CMake files are perfectly formatted!"
        exit 0
    else
        echo "------------------------------------------------------"
        echo "❌ ERROR: CMake formatting violations found."
        echo "Run this script with --apply to fix them."
        echo "------------------------------------------------------"
        exit 1
    fi
fi