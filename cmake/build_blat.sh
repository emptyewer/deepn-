#!/bin/bash
set -e
BLAT_DIR="$1"
NPROC="$2"

# Find libpng via brew or fall back to /opt/homebrew
PNG_PREFIX="${3:-$(brew --prefix libpng 2>/dev/null || echo /opt/homebrew)}"

make -C "$BLAT_DIR" -j"$NPROC" \
    "CFLAGS=-D_STATIC -I${PNG_PREFIX}/include" \
    "L=-L${PNG_PREFIX}/lib -lpng"
