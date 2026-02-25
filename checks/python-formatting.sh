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

echo "🔍 Finding Python files (respecting .gitignore and submodules)..."

# Find tracked python files
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILES=$(git ls-files | grep -E '\.py$' | grep -vEf "$SCRIPT_DIR/format_ignore.txt")

if [[ -z "$FILES" ]]; then
    echo "No Python files found."
    exit 0
fi

if [ "$APPLY_FIXES" = true ]; then
    echo "🛠  Running black to fix formatting..."
    # Black formats in-place by default
    echo "$FILES" | xargs black
    echo "✅ Python formatting complete."
else
    echo "🧪 Checking Python formatting (Dry Run)..."
    
    # --check: return exit code 1 if files need formatting
    # --diff: show what would change
    # --quiet: reduce noise, we only want the diffs/errors
    if echo "$FILES" | xargs black --check --diff; then
        echo "✨ Python code is looking sharp!"
        exit 0
    else
        echo "------------------------------------------------------"
        echo "❌ ERROR: Python formatting violations found."
        echo "Run this script with --apply to fix them."
        echo "------------------------------------------------------"
        exit 1
    fi
fi