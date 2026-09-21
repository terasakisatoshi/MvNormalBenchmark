# MvNormal benchmark report

Generated: `2026-09-21T08:03:04+00:00`
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
| 1 | `rust-ziggurat` | 0.025 | 9.188 | 8.824 | 1.00x | `2021309.479503571056` |
| 2 | `julia` | 0.076 | 9.254 | 9.237 | 1.01x | `2.0135557629042931e6` |
| 3 | `cxx-ziggurat` | 0.026 | 9.295 | 8.953 | 1.01x | `2021309.4795035999` |
| 4 | `julia-ziggurat` | 0.077 | 9.768 | 9.723 | 1.06x | `2.021309479503571e6` |
| 5 | `julia-polar` | 0.077 | 14.294 | 14.218 | 1.56x | `2.0152758700649042e6` |
| 6 | `fortran-ziggurat` | 0.025 | 14.324 | 14.169 | 1.56x | `2.02130947950357106E+006` |
| 7 | `julia-distributions` | 0.042 | 15.350 | 14.893 | 1.67x | `2.0135557629042931e6` |
| 8 | `cxx-polar` | 0.055 | 17.586 | 14.688 | 1.91x | `2015275.8700647992` |
| 9 | `rust-polar` | 0.090 | 20.884 | 14.912 | 2.27x | `2015275.870064904215` |
| 10 | `fortran-polar` | 0.068 | 21.923 | 17.007 | 2.39x | `2.01527587006490421E+006` |

## Notes

- The custom `polar` and `ziggurat` paths use the common xorshift64 seed `0x5EED2021`; checksums can still differ because floating-point operations, Cholesky factors, and accumulation order differ.
- Julia's default and `Distributions.jl` rows use Julia/Random and are not part of the common custom-RNG stream.
- Compilation time is not included in the sampling timing.
- The Julia official row is included when `Distributions.jl` is available in `julia/Project.toml`.
