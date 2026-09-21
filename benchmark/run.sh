#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$ROOT_DIR/benchmark/.build"
OUTPUT_FILE="${1:-$ROOT_DIR/benchmark/results.csv}"
JULIA_PROJECT_ARGS=(--project="$ROOT_DIR/julia")

DIM="${DIM:-64}"
SAMPLES="${SAMPLES:-20000}"
REPEATS="${REPEATS:-5}"

mkdir -p "$BUILD_DIR"
mkdir -p "$(dirname "$OUTPUT_FILE")"

printf 'language,dim,samples,repeats,setup_sec,avg_sample_sec,min_sample_sec,checksum\n' > "$OUTPUT_FILE"

run_julia() {
    if ! command -v julia >/dev/null 2>&1; then
        echo 'Skipping Julia: julia was not found.' >&2
        return 0
    fi

    julia "${JULIA_PROJECT_ARGS[@]}" "$ROOT_DIR/julia/benchmark.jl" \
        --dim "$DIM" --samples "$SAMPLES" --repeats "$REPEATS" >> "$OUTPUT_FILE"
}

run_julia_distributions() {
    if ! command -v julia >/dev/null 2>&1; then
        return 0
    fi
    if ! julia "${JULIA_PROJECT_ARGS[@]}" -e 'using Distributions' >/dev/null 2>&1; then
        echo 'Skipping Julia Distributions.jl: package was not found.' >&2
        return 0
    fi

    julia "${JULIA_PROJECT_ARGS[@]}" "$ROOT_DIR/julia/benchmark_distributions.jl" \
        --dim "$DIM" --samples "$SAMPLES" --repeats "$REPEATS" >> "$OUTPUT_FILE"
}

run_cxx() {
    local compiler="${CXX:-g++}"
    if ! command -v "$compiler" >/dev/null 2>&1; then
        echo "Skipping C++: $compiler was not found." >&2
        return 0
    fi

    "$compiler" -std=c++23 -O3 -DNDEBUG \
        -I "$ROOT_DIR/cxx" \
        "$ROOT_DIR/cxx/benchmark.cpp" \
        -o "$BUILD_DIR/mvnormal_cxx"

    "$BUILD_DIR/mvnormal_cxx" \
        --dim "$DIM" --samples "$SAMPLES" --repeats "$REPEATS" >> "$OUTPUT_FILE"
}

run_fortran() {
    local compiler="${FC:-gfortran}"
    if ! command -v "$compiler" >/dev/null 2>&1; then
        echo "Skipping Fortran: $compiler was not found." >&2
        return 0
    fi

    "$compiler" -std=f2023 -O3 \
        -J "$BUILD_DIR" \
        "$ROOT_DIR/fortran/mvnormal.f90" \
        "$ROOT_DIR/fortran/benchmark.f90" \
        -o "$BUILD_DIR/mvnormal_fortran"

    "$BUILD_DIR/mvnormal_fortran" \
        --dim "$DIM" --samples "$SAMPLES" --repeats "$REPEATS" >> "$OUTPUT_FILE"
}

run_rust() {
    if ! command -v cargo >/dev/null 2>&1; then
        echo 'Skipping Rust: cargo was not found.' >&2
        return 0
    fi

    CARGO_TARGET_DIR="$BUILD_DIR/rust-target" \
        cargo build --release --manifest-path "$ROOT_DIR/rust/Cargo.toml"

    "$BUILD_DIR/rust-target/release/mvnormal_benchmark" \
        --dim "$DIM" --samples "$SAMPLES" --repeats "$REPEATS" >> "$OUTPUT_FILE"
}

run_julia
run_julia_distributions
run_cxx
run_fortran
run_rust

echo "Benchmark results written to $OUTPUT_FILE"
