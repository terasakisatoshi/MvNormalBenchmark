# MvNormal benchmark report

Generated: `2026-09-21T11:51:11+00:00`
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
| 1 | `julia-distributions-batch` | 0.052 | 6.200 | 6.138 | 1.00x | `2.0127889682539785e6` |
| 2 | `julia-batch` | 0.135 | 6.343 | 6.250 | 1.02x | `2.0127889682539785e6` |
| 3 | `julia-ziggurat-batch` | 0.082 | 7.614 | 7.547 | 1.23x | `2.0213094795035578e6` |
| 4 | `julia` | 0.083 | 8.903 | 8.831 | 1.44x | `2.0135557629042931e6` |
| 5 | `rust-ziggurat` | 0.029 | 9.436 | 9.400 | 1.52x | `2021309.479503571056` |
| 6 | `cxx-ziggurat` | 0.026 | 9.624 | 9.495 | 1.55x | `2021309.4795035999` |
| 7 | `julia-ziggurat` | 0.072 | 9.633 | 9.583 | 1.55x | `2.021309479503571e6` |
| 8 | `julia-polar-batch` | 0.073 | 11.707 | 11.218 | 1.89x | `2.015275870064887e6` |
| 9 | `fortran-ziggurat` | 0.026 | 12.310 | 12.243 | 1.99x | `2.02130947950357106E+006` |
| 10 | `fortran-stdlib` | 0.053 | 12.903 | 12.870 | 2.08x | `2.01738156312306272E+006` |
| 11 | `julia-polar` | 0.070 | 14.344 | 14.289 | 2.31x | `2.0152758700649042e6` |
| 12 | `python-polar` | 0.139 | 14.901 | 14.775 | 2.40x | `2015275.8700648781` |
| 13 | `python-ziggurat` | 0.131 | 15.547 | 15.475 | 2.51x | `2021309.4795035562` |
| 14 | `rust-statrs` | 0.122 | 16.600 | 16.539 | 2.68x | `2022495.290060152067` |
| 15 | `rust-polar` | 0.063 | 17.803 | 15.050 | 2.87x | `2015275.870064904215` |
| 16 | `julia-distributions` | 0.040 | 18.206 | 17.842 | 2.94x | `2.0135557629042934e6` |
| 17 | `cxx-polar` | 0.075 | 18.219 | 14.673 | 2.94x | `2015275.8700647992` |
| 18 | `fortran-polar` | 0.053 | 20.566 | 17.773 | 3.32x | `2.01527587006490421E+006` |
| 19 | `rust-rand-distr` | 0.028 | 25.659 | 25.533 | 4.14x | `2022495.290060150903` |

## Notes

- The custom `polar` and `ziggurat` paths use the common xorshift64 seed `0x5EED2021`; checksums can still differ because floating-point operations, Cholesky factors, and accumulation order differ.
- Julia's default and `Distributions.jl` rows use Julia/Random and are not part of the common custom-RNG stream.
- The Rust `rand-distr` and `statrs` rows use `StdRng` (and `rand_distr`'s Zignor Ziggurat inside `statrs`) and are not part of the common custom-RNG stream.
- The Fortran `stdlib` row uses `stdlib_stats_distribution_normal` and `stdlib_linalg`'s Cholesky and is not part of the common custom-RNG stream.
- The Julia `*-batch` rows generate every sample of a repeat in one batched call (`mul!`/BLAS with `BLAS.set_num_threads(1)`), so their timing is a batched rather than per-sample comparison.
- Compilation time is not included in the sampling timing.
- The Julia official row is included when `Distributions.jl` is available in `julia/Project.toml`.
