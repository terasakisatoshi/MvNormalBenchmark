# MvNormal benchmark report

Generated: `2026-09-21T06:32:47+00:00`
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
| 1 | `julia` | 1.093 | 9.061 | 9.003 | 1.00x | `-52654.41119174457` |
| 2 | `julia-distributions` | 0.026 | 15.112 | 14.663 | 1.67x | `-52654.411191744526` |
| 3 | `rust` | 0.048 | 17.434 | 15.284 | 1.92x | `2010703.010879529640` |
| 4 | `cxx` | 0.041 | 19.028 | 16.895 | 2.10x | `-2609.6076536198825` |
| 5 | `fortran` | 0.061 | 24.680 | 20.336 | 2.72x | `3.15234645558872772E+006` |

## Notes

- Random-number generators differ by language, so checksums are not expected to match across implementations.
- Compilation time is not included in the sampling timing.
- The Julia official row is included when `Distributions.jl` is available in `julia/Project.toml`.
