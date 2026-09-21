# MvNormal benchmark report

Generated: `2026-09-21T11:46:30+00:00`
Input: `benchmark/results.csv`

The primary comparison is `avg_sample_sec`: the elapsed time for all samples in one repeat, averaged over repeats. Lower is faster.

## Conditions

- Dimension: `64`
- Samples per repeat: `20000`
- Repeats: `5`
- Setup time includes covariance construction and Cholesky setup where the implementation performs them inside the measured setup section.

## Speed comparison

| Rank | Language | Setup (ms) | Avg sampling (ms) | Min sampling (ms) | Relative to fastest | Checksum |
|---:|---|---:|---:|---:|---:|---:|
| 1 | `julia-distributions-batch` | 0.049 | 6.354 | 6.278 | 1.00x | `2.0127889682539785e6` |
| 2 | `julia-batch` | 0.079 | 7.555 | 7.456 | 1.19x | `2.0127889682539785e6` |
| 3 | `julia-ziggurat-batch` | 0.067 | 8.395 | 8.320 | 1.32x | `2.0213094795035578e6` |
| 4 | `julia` | 0.075 | 9.327 | 9.234 | 1.47x | `2.0135557629042931e6` |
| 5 | `cxx-ziggurat` | 0.028 | 9.580 | 9.491 | 1.51x | `2021309.4795035999` |
| 6 | `rust-ziggurat` | 0.030 | 9.596 | 9.528 | 1.51x | `2021309.479503571056` |
| 7 | `julia-ziggurat` | 0.075 | 9.892 | 9.812 | 1.56x | `2.021309479503571e6` |
| 8 | `fortran-ziggurat` | 0.028 | 12.358 | 12.280 | 1.94x | `2.02130947950357106E+006` |
| 9 | `julia-polar-batch` | 0.069 | 12.501 | 12.419 | 1.97x | `2.015275870064887e6` |
| 10 | `fortran-stdlib` | 0.040 | 13.025 | 12.933 | 2.05x | `2.01738156312306272E+006` |
| 11 | `julia-polar` | 0.069 | 14.423 | 14.309 | 2.27x | `2.0152758700649042e6` |
| 12 | `cxx-polar` | 0.033 | 15.031 | 14.576 | 2.37x | `2015275.8700647992` |
| 13 | `python-polar` | 0.133 | 15.047 | 15.013 | 2.37x | `2015275.8700648781` |
| 14 | `python-ziggurat` | 0.133 | 15.629 | 15.416 | 2.46x | `2021309.4795035562` |
| 15 | `julia-distributions` | 0.041 | 16.027 | 15.931 | 2.52x | `2.0135557629042934e6` |
| 16 | `rust-polar` | 0.039 | 16.512 | 15.210 | 2.60x | `2015275.870064904215` |
| 17 | `rust-statrs` | 0.124 | 16.533 | 16.455 | 2.60x | `2022495.290060152067` |
| 18 | `fortran-polar` | 0.039 | 19.075 | 17.873 | 3.00x | `2.01527587006490421E+006` |
| 19 | `rust-rand-distr` | 0.027 | 25.526 | 25.445 | 4.02x | `2022495.290060150903` |

## Notes

- The custom `polar` and `ziggurat` paths use the common xorshift64 seed `0x5EED2021`; checksums can still differ because floating-point operations, Cholesky factors, and accumulation order differ.
- Julia's default and `Distributions.jl` rows use Julia/Random and are not part of the common custom-RNG stream.
- The Rust `rand-distr` and `statrs` rows use `StdRng` (and `rand_distr`'s Zignor Ziggurat inside `statrs`) and are not part of the common custom-RNG stream.
- The Fortran `stdlib` row uses `stdlib_stats_distribution_normal` and `stdlib_linalg`'s Cholesky and is not part of the common custom-RNG stream.
- The Julia `*-batch` rows generate every sample of a repeat in one batched call (`mul!`/BLAS with `BLAS.set_num_threads(1)`), so their timing is a batched rather than per-sample comparison.
- Compilation time is not included in the sampling timing.
- The Julia official row is included when `Distributions.jl` is available in `julia/Project.toml`.
