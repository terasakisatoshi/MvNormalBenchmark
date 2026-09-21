# MvNormal benchmark report

Generated: `2026-09-21T11:36:52+00:00`
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
| 1 | `julia` | 0.074 | 9.211 | 9.027 | 1.00x | `2.0135557629042931e6` |
| 2 | `cxx-ziggurat` | 0.029 | 9.530 | 9.482 | 1.03x | `2021309.4795035999` |
| 3 | `rust-ziggurat` | 0.026 | 9.639 | 9.489 | 1.05x | `2021309.479503571056` |
| 4 | `julia-ziggurat` | 0.076 | 9.830 | 9.732 | 1.07x | `2.021309479503571e6` |
| 5 | `fortran-ziggurat` | 0.028 | 12.263 | 12.189 | 1.33x | `2.02130947950357106E+006` |
| 6 | `fortran-stdlib` | 0.039 | 12.878 | 12.811 | 1.40x | `2.01738156312306272E+006` |
| 7 | `julia-polar` | 0.083 | 14.437 | 14.412 | 1.57x | `2.0152758700649042e6` |
| 8 | `julia-distributions` | 0.060 | 15.001 | 14.629 | 1.63x | `2.0135557629042931e6` |
| 9 | `python-polar` | 0.136 | 15.005 | 14.934 | 1.63x | `2015275.8700648781` |
| 10 | `python-ziggurat` | 0.136 | 15.485 | 15.418 | 1.68x | `2021309.4795035562` |
| 11 | `rust-statrs` | 0.129 | 16.541 | 16.505 | 1.80x | `2022495.290060152067` |
| 12 | `cxx-polar` | 0.073 | 19.160 | 14.561 | 2.08x | `2015275.8700647992` |
| 13 | `rust-polar` | 0.086 | 20.714 | 15.094 | 2.25x | `2015275.870064904215` |
| 14 | `fortran-polar` | 0.063 | 21.922 | 17.795 | 2.38x | `2.01527587006490421E+006` |
| 15 | `rust-rand-distr` | 0.028 | 25.500 | 25.316 | 2.77x | `2022495.290060150903` |

## Notes

- The custom `polar` and `ziggurat` paths use the common xorshift64 seed `0x5EED2021`; checksums can still differ because floating-point operations, Cholesky factors, and accumulation order differ.
- Julia's default and `Distributions.jl` rows use Julia/Random and are not part of the common custom-RNG stream.
- The Rust `rand-distr` and `statrs` rows use `StdRng` (and `rand_distr`'s Zignor Ziggurat inside `statrs`) and are not part of the common custom-RNG stream.
- The Fortran `stdlib` row uses `stdlib_stats_distribution_normal` and `stdlib_linalg`'s Cholesky and is not part of the common custom-RNG stream.
- Compilation time is not included in the sampling timing.
- The Julia official row is included when `Distributions.jl` is available in `julia/Project.toml`.
