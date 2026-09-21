# MvNormal benchmark report

Generated: `2026-09-21T11:06:51+00:00`
Input: `/Users/terasaki/work/terasakisatoshi/MvNormalBenchmark/benchmark/results.csv`

The primary comparison is `avg_sample_sec`: the elapsed time for all samples in one repeat, averaged over repeats. Lower is faster.

## Conditions

- Dimension: `64`
- Samples per repeat: `20000`
- Repeats: `5`
- Setup time includes covariance construction and Cholesky setup where the implementation performs them inside the measured setup section.

## Speed comparison

| Rank | Language | Setup (ms) | Avg sampling (ms) | Min sampling (ms) | Relative to fastest | Checksum |
|---:|---|---:|---:|---:|---:|---:|
| 1 | `julia` | 0.070 | 9.376 | 9.342 | 1.00x | `2.0135557629042931e6` |
| 2 | `rust-ziggurat` | 0.029 | 9.619 | 9.563 | 1.03x | `2021309.479503571056` |
| 3 | `cxx-ziggurat` | 0.026 | 9.662 | 9.555 | 1.03x | `2021309.4795035999` |
| 4 | `julia-ziggurat` | 0.107 | 9.861 | 9.769 | 1.05x | `2.021309479503571e6` |
| 5 | `fortran-ziggurat` | 0.028 | 12.372 | 12.281 | 1.32x | `2.02130947950357106E+006` |
| 6 | `julia-polar` | 0.087 | 14.432 | 14.336 | 1.54x | `2.0152758700649042e6` |
| 7 | `julia-distributions` | 0.045 | 15.075 | 14.873 | 1.61x | `2.0135557629042931e6` |
| 8 | `python-polar` | 0.132 | 15.085 | 14.978 | 1.61x | `2015275.8700648781` |
| 9 | `python-ziggurat` | 0.133 | 15.657 | 15.602 | 1.67x | `2021309.4795035562` |
| 10 | `cxx-polar` | 0.043 | 16.213 | 14.635 | 1.73x | `2015275.8700647992` |
| 11 | `rust-polar` | 0.026 | 16.329 | 15.055 | 1.74x | `2015275.870064904215` |
| 12 | `rust-statrs` | 0.117 | 16.815 | 16.737 | 1.79x | `2022495.290060152067` |
| 13 | `fortran-polar` | 0.038 | 18.730 | 17.173 | 2.00x | `2.01527587006490421E+006` |
| 14 | `rust-rand-distr` | 0.025 | 25.633 | 25.514 | 2.73x | `2022495.290060150903` |

## Notes

- The custom `polar` and `ziggurat` paths use the common xorshift64 seed `0x5EED2021`; checksums can still differ because floating-point operations, Cholesky factors, and accumulation order differ.
- Julia's default and `Distributions.jl` rows use Julia/Random and are not part of the common custom-RNG stream.
- The Rust `rand-distr` and `statrs` rows use `StdRng` (and `rand_distr`'s Zignor Ziggurat inside `statrs`) and are not part of the common custom-RNG stream.
- Compilation time is not included in the sampling timing.
- The Julia official row is included when `Distributions.jl` is available in `julia/Project.toml`.
