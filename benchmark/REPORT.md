# MvNormal benchmark report

Generated: `2026-09-21T08:40:38+00:00`
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
| 1 | `julia` | 0.083 | 9.270 | 9.229 | 1.00x | `2.0135557629042931e6` |
| 2 | `cxx-ziggurat` | 0.026 | 9.569 | 9.534 | 1.03x | `2021309.4795035999` |
| 3 | `rust-ziggurat` | 0.042 | 10.127 | 9.851 | 1.09x | `2021309.479503571056` |
| 4 | `fortran-ziggurat` | 0.080 | 12.684 | 12.340 | 1.37x | `2.02130947950357106E+006` |
| 5 | `julia-ziggurat` | 0.090 | 13.916 | 9.808 | 1.50x | `2.021309479503571e6` |
| 6 | `julia-polar` | 0.080 | 15.293 | 14.987 | 1.65x | `2.0152758700649042e6` |
| 7 | `julia-distributions` | 0.043 | 16.323 | 15.445 | 1.76x | `2.0135557629042931e6` |
| 8 | `cxx-polar` | 0.053 | 17.680 | 14.753 | 1.91x | `2015275.8700647992` |
| 9 | `fortran-polar` | 0.043 | 19.849 | 18.477 | 2.14x | `2.01527587006490421E+006` |
| 10 | `rust-polar` | 0.038 | 22.164 | 17.025 | 2.39x | `2015275.870064904215` |

## Notes

- The custom `polar` and `ziggurat` paths use the common xorshift64 seed `0x5EED2021`; checksums can still differ because floating-point operations, Cholesky factors, and accumulation order differ.
- Julia's default and `Distributions.jl` rows use Julia/Random and are not part of the common custom-RNG stream.
- Compilation time is not included in the sampling timing.
- The Julia official row is included when `Distributions.jl` is available in `julia/Project.toml`.
