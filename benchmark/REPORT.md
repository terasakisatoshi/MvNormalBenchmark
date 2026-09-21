# MvNormal benchmark report

Generated: `2026-09-21T05:44:32+00:00`
Input: `/Users/terasaki/work/terasakisatoshi/MvNormal/benchmark/results.csv`

The primary comparison is `avg_sample_sec`: the elapsed time for all samples in one repeat, averaged over repeats. Lower is faster.

## Conditions

- Dimension: `64`
- Samples per repeat: `20000`
- Repeats: `5`
- Setup time includes covariance construction and Cholesky setup where the implementation performs them inside the measured setup section.

## Speed comparison

| Rank | Language | Setup (ms) | Avg sampling (ms) | Min sampling (ms) | Relative to fastest | Checksum |
|---:|---|---:|---:|---:|---:|---:|
| 1 | `julia` | 0.064 | 15.578 | 15.556 | 1.00x | `558932.2140488301` |
| 2 | `julia-distributions` | 0.024 | 15.706 | 15.526 | 1.01x | `558932.2140488295` |
| 3 | `cxx` | 0.073 | 21.320 | 16.038 | 1.37x | `-2609.6076536198825` |
| 4 | `fortran` | 0.087 | 31.560 | 26.024 | 2.03x | `3.14834553375710314E+006` |
| 5 | `rust` | 0.078 | 35.881 | 30.163 | 2.30x | `2008496.332629454555` |

## Notes

- Random-number generators differ by language, so checksums are not expected to match across implementations.
- Compilation time is not included in the sampling timing.
- The Julia official row is included when `Distributions.jl` is available in `julia/Project.toml`.
