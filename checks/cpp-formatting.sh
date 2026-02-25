#!/bin/bash

APPLY_FIXES=false
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --apply) APPLY_FIXES=true ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Find files tracked by git, excluding submodules and ignored files
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILES=$(git ls-files | grep -E '\.(h(xx)?|c(xx?))$' | grep -vEf "$SCRIPT_DIR/format_ignore.txt")

if [[ -z "$FILES" ]]; then
    echo "No matching files found."
    exit 0
fi

if [ "$APPLY_FIXES" = true ]; then
    echo "🛠  Applying clang-format fixes in-place..."
    echo "$FILES" | xargs clang-format -i
    echo "✅ Formatting complete."
else
    echo "🔍 Checking formatting (Dry Run)..."
    
    # --dry-run: Don't change files
    # --Werror: Exit with non-zero if formatting is needed
    # We redirect stderr to a temp file to show you what's wrong
    if echo "$FILES" | xargs clang-format --dry-run --Werror 2>/tmp/fmt_err; then
        echo "✨ Everything is perfectly formatted!"
        exit 0
    else
        echo "❌ Formatting errors found!"
        cat /tmp/fmt_err
        echo ""
        echo "Run this script with --apply to fix these issues."
        exit 1
    fi
fi