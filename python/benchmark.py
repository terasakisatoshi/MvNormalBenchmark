#!/usr/bin/env python3
"""Benchmark the Python MvNormal implementation under the common conditions.

Run with the project environment, for example:

    uv run --project python python/benchmark.py --dim 64 --samples 20000 --repeats 5
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from mvnormal import MvNormal, NormalAlgorithm, NormalRng

COMPARISON_NORMAL_SEED = 0x5EED2021
WARMUP_SEED = 0xABAD1DEA


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="MvNormal Python benchmark")
    parser.add_argument("--dim", type=int, default=10)
    parser.add_argument("--samples", type=int, default=10_000)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument(
        "--normal",
        choices=("polar", "ziggurat"),
        default="polar",
        help="standard normal algorithm (default: polar)",
    )
    arguments = parser.parse_args()
    if arguments.dim <= 0:
        parser.error("--dim must be positive")
    if arguments.samples <= 0:
        parser.error("--samples must be positive")
    if arguments.repeats <= 0:
        parser.error("--repeats must be positive")
    return arguments


def build_inputs(dimension: int) -> tuple[np.ndarray, np.ndarray]:
    indices = np.arange(dimension, dtype=np.float64)
    mean = 0.01 * indices
    covariance = 0.25 ** np.abs(indices[:, None] - indices[None, :])
    return mean, covariance


def main() -> None:
    arguments = parse_args()
    algorithm = NormalAlgorithm.MARSAGLIA_POLAR
    if arguments.normal == "ziggurat":
        algorithm = NormalAlgorithm.ZIGGURAT

    mean, covariance = build_inputs(arguments.dim)

    setup_start = time.perf_counter()
    distribution = MvNormal(mean, covariance)
    setup_sec = time.perf_counter() - setup_start

    # Compile the Numba kernels before timing and keep the comparison seed
    # untouched by using a separate generator for the warmup.
    warmup_rng = NormalRng(WARMUP_SEED, algorithm)
    distribution.sample_checksum(warmup_rng, min(arguments.samples, 64))

    rng = NormalRng(COMPARISON_NORMAL_SEED, algorithm)
    durations = []
    checksum = 0.0
    for _ in range(arguments.repeats):
        start = time.perf_counter()
        checksum += distribution.sample_checksum(rng, arguments.samples)
        durations.append(time.perf_counter() - start)

    avg_sample_sec = sum(durations) / arguments.repeats
    min_sample_sec = min(durations)
    print(
        f"python-{algorithm.label()},{arguments.dim},{arguments.samples},"
        f"{arguments.repeats},{setup_sec:.9g},{avg_sample_sec:.9g},"
        f"{min_sample_sec:.9g},{checksum!r}"
    )


if __name__ == "__main__":
    main()
