# MvNormal benchmark report

Generated: `2026-09-21T08:41:49+00:00`
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
| 1 | `julia` | 0.068 | 9.285 | 9.244 | 1.00x | `2.0135557629042931e6` |
| 2 | `rust-ziggurat` | 0.030 | 9.381 | 8.998 | 1.01x | `2021309.479503571056` |
| 3 | `cxx-ziggurat` | 0.027 | 9.576 | 9.514 | 1.03x | `2021309.4795035999` |
| 4 | `julia-ziggurat` | 0.079 | 9.772 | 9.696 | 1.05x | `2.021309479503571e6` |
| 5 | `fortran-ziggurat` | 0.025 | 12.188 | 11.953 | 1.31x | `2.02130947950357106E+006` |
| 6 | `julia-polar` | 0.078 | 14.370 | 14.257 | 1.55x | `2.0152758700649042e6` |
| 7 | `julia-distributions` | 0.079 | 15.186 | 14.740 | 1.64x | `2.0135557629042931e6` |
| 8 | `cxx-polar` | 0.027 | 16.191 | 14.668 | 1.74x | `2015275.8700647992` |
| 9 | `rust-polar` | 0.074 | 18.619 | 14.960 | 2.01x | `2015275.870064904215` |
| 10 | `fortran-polar` | 0.052 | 19.973 | 17.780 | 2.15x | `2.01527587006490421E+006` |

## Notes

- The custom `polar` and `ziggurat` paths use the common xorshift64 seed `0x5EED2021`; checksums can still differ because floating-point operations, Cholesky factors, and accumulation order differ.
- Julia's default and `Distributions.jl` rows use Julia/Random and are not part of the common custom-RNG stream.
- Compilation time is not included in the sampling timing.
- The Julia official row is included when `Distributions.jl` is available in `julia/Project.toml`.
