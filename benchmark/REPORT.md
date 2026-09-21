# MvNormal benchmark report

Generated: `2026-09-21T08:17:37+00:00`
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
| 1 | `julia` | 0.834 | 9.215 | 9.175 | 1.00x | `2.0135557629042931e6` |
| 2 | `rust-ziggurat` | 0.026 | 9.262 | 8.889 | 1.01x | `2021309.479503571056` |
| 3 | `cxx-ziggurat` | 0.030 | 9.294 | 9.051 | 1.01x | `2021309.4795035999` |
| 4 | `julia-ziggurat` | 0.096 | 9.733 | 9.647 | 1.06x | `2.021309479503571e6` |
| 5 | `fortran-ziggurat` | 0.025 | 13.406 | 12.980 | 1.45x | `2.02130947950357106E+006` |
| 6 | `julia-polar` | 0.078 | 14.344 | 14.280 | 1.56x | `2.0152758700649042e6` |
| 7 | `julia-distributions` | 0.040 | 15.260 | 14.616 | 1.66x | `2.0135557629042931e6` |
| 8 | `rust-polar` | 0.070 | 19.487 | 14.761 | 2.11x | `2015275.870064904215` |
| 9 | `cxx-polar` | 0.082 | 19.823 | 14.366 | 2.15x | `2015275.8700647992` |
| 10 | `fortran-polar` | 0.062 | 21.061 | 17.368 | 2.29x | `2.01527587006490421E+006` |

## Notes

- The custom `polar` and `ziggurat` paths use the common xorshift64 seed `0x5EED2021`; checksums can still differ because floating-point operations, Cholesky factors, and accumulation order differ.
- Julia's default and `Distributions.jl` rows use Julia/Random and are not part of the common custom-RNG stream.
- Compilation time is not included in the sampling timing.
- The Julia official row is included when `Distributions.jl` is available in `julia/Project.toml`.
