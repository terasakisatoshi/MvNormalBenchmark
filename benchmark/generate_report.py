#!/usr/bin/env python3
"""Generate a Markdown speed comparison from benchmark CSV output."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path


REQUIRED_COLUMNS = {
    "language",
    "dim",
    "samples",
    "repeats",
    "setup_sec",
    "avg_sample_sec",
    "min_sample_sec",
    "checksum",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a Markdown report from benchmark/results.csv"
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="benchmark/results.csv",
        help="benchmark CSV input (default: benchmark/results.csv)",
    )
    parser.add_argument(
        "output",
        nargs="?",
        default="benchmark/REPORT.md",
        help="Markdown output (default: benchmark/REPORT.md)",
    )
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        columns = set(reader.fieldnames or ())
        missing = REQUIRED_COLUMNS - columns
        if missing:
            missing_text = ", ".join(sorted(missing))
            raise ValueError(f"CSV is missing columns: {missing_text}")
        rows = list(reader)

    if not rows:
        raise ValueError("CSV contains no benchmark rows")
    return rows


def seconds_to_ms(value: str) -> str:
    return f"{float(value) * 1_000.0:.3f}"


def generate_report(rows: list[dict[str, str]], source: Path) -> str:
    ordered = sorted(rows, key=lambda row: float(row["avg_sample_sec"]))
    fastest = float(ordered[0]["avg_sample_sec"])
    first = rows[0]
    generated_at = datetime.now(timezone.utc).isoformat(timespec="seconds")

    lines = [
        "# MvNormal benchmark report",
        "",
        f"Generated: `{generated_at}`",
        f"Input: `{source}`",
        "",
        "The primary comparison is `avg_sample_sec`: the elapsed time for all "
        "samples in one repeat, averaged over repeats. Lower is faster.",
        "",
        "## Conditions",
        "",
        f"- Dimension: `{first['dim']}`",
        f"- Samples per repeat: `{first['samples']}`",
        f"- Repeats: `{first['repeats']}`",
        "- Setup time includes covariance construction and Cholesky setup where "
        "the implementation performs them inside the measured setup section.",
        "",
        "## Speed comparison",
        "",
        "| Rank | Language | Setup (ms) | Avg sampling (ms) | Min sampling (ms) | Relative to fastest | Checksum |",
        "|---:|---|---:|---:|---:|---:|---:|",
    ]

    for rank, row in enumerate(ordered, start=1):
        average = float(row["avg_sample_sec"])
        relative = average / fastest if fastest > 0.0 else float("nan")
        lines.append(
            f"| {rank} | `{row['language']}` | "
            f"{seconds_to_ms(row['setup_sec'])} | "
            f"{seconds_to_ms(row['avg_sample_sec'])} | "
            f"{seconds_to_ms(row['min_sample_sec'])} | "
            f"{relative:.2f}x | `{row['checksum']}` |"
        )

    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- The custom `polar` and `ziggurat` paths use the common xorshift64 "
            "seed `0x5EED2021`; checksums can still differ because floating-point "
            "operations, Cholesky factors, and accumulation order differ.",
            "- Julia's default and `Distributions.jl` rows use Julia/Random and "
            "are not part of the common custom-RNG stream.",
            "- The Rust `rand-distr` and `statrs` rows use `StdRng` (and "
            "`rand_distr`'s Zignor Ziggurat inside `statrs`) and are not part of "
            "the common custom-RNG stream.",
            "- The Fortran `stdlib` row uses `stdlib_stats_distribution_normal` "
            "and `stdlib_linalg`'s Cholesky and is not part of the common "
            "custom-RNG stream.",
            "- The Julia `*-batch` rows generate every sample of a repeat in one "
            "batched call (`lmul!`/BLAS with `BLAS.set_num_threads(1)`), so their "
            "timing is a batched rather than per-sample comparison.",
            "- Compilation time is not included in the sampling timing.",
            "- The Julia official row is included when `Distributions.jl` is "
            "available in `julia/Project.toml`.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    arguments = parse_args()
    source = Path(arguments.input)
    destination = Path(arguments.output)
    report = generate_report(read_rows(source), source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(report, encoding="utf-8")
    print(f"Markdown report written to {destination}")


if __name__ == "__main__":
    main()
