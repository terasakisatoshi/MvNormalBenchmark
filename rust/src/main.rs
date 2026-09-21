use std::env;
use std::process;
use std::time::Instant;

use mvnormal_benchmark::{MvNormal, NormalAlgorithm, StandardRng};
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::{Distribution, StandardNormal};
use statrs::distribution::MultivariateNormal as StatrsMvNormal;

const COMPARISON_NORMAL_SEED: u64 = 0x5EED_2021;

#[derive(Clone, Copy, PartialEq, Eq)]
enum NormalChoice {
    Polar,
    Ziggurat,
    RandDistr,
    Statrs,
}

impl NormalChoice {
    fn label(self) -> &'static str {
        match self {
            Self::Polar => "polar",
            Self::Ziggurat => "ziggurat",
            Self::RandDistr => "rand-distr",
            Self::Statrs => "statrs",
        }
    }
}

struct Config {
    dimension: usize,
    samples: usize,
    repeats: usize,
    normal_algorithm: NormalChoice,
}

fn usage() -> &'static str {
    "usage: mvnormal_benchmark --dim <dimension> --samples <samples> --repeats <repeats> --normal <polar|ziggurat|rand-distr|statrs>"
}

fn parse_positive(value: &str, name: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("{name} must be a positive integer, got {value:?}"))?;
    if parsed == 0 {
        return Err(format!("{name} must be greater than zero"));
    }
    Ok(parsed)
}

fn parse_args() -> Result<Config, String> {
    let mut dimension = None;
    let mut samples = None;
    let mut repeats = None;
    let mut normal_algorithm = NormalChoice::Polar;
    let mut arguments = env::args().skip(1);

    while let Some(flag) = arguments.next() {
        let value = arguments
            .next()
            .ok_or_else(|| format!("missing value for {flag}\n{}", usage()))?;
        match flag.as_str() {
            "--dim" => dimension = Some(parse_positive(&value, "--dim")?),
            "--samples" => samples = Some(parse_positive(&value, "--samples")?),
            "--repeats" => repeats = Some(parse_positive(&value, "--repeats")?),
            "--normal" => {
                normal_algorithm = match value.as_str() {
                    "polar" => NormalChoice::Polar,
                    "ziggurat" => NormalChoice::Ziggurat,
                    "rand-distr" | "rand_distr" => NormalChoice::RandDistr,
                    "statrs" => NormalChoice::Statrs,
                    _ => {
                        return Err(format!(
                            "--normal must be polar, ziggurat, rand-distr, or statrs, got {value:?}"
                        ));
                    }
                }
            }
            _ => return Err(format!("unknown argument {flag}\n{}", usage())),
        }
    }

    Ok(Config {
        dimension: dimension.ok_or_else(|| format!("missing --dim\n{}", usage()))?,
        samples: samples.ok_or_else(|| format!("missing --samples\n{}", usage()))?,
        repeats: repeats.ok_or_else(|| format!("missing --repeats\n{}", usage()))?,
        normal_algorithm,
    })
}

fn distance(row: usize, column: usize) -> i32 {
    row.abs_diff(column) as i32
}

fn run(config: Config) -> Result<(), String> {
    let dimension = config.dimension;
    let mean = (0..dimension)
        .map(|index| (index as f64) * 0.01)
        .collect::<Vec<_>>();
    // Row-major flattened covariance. The matrix is symmetric, so the same
    // vector works for both the row-major and column-major consumers.
    let covariance_flat = (0..dimension)
        .flat_map(|row| {
            (0..dimension).map(move |column| 2.0_f64.powi(-2 * distance(row, column)))
        })
        .collect::<Vec<_>>();

    let mut output = vec![0.0; dimension];
    let mut durations = Vec::with_capacity(config.repeats);
    let mut checksum = 0.0;
    let setup_sec;

    if config.normal_algorithm == NormalChoice::Statrs {
        // `statrs::distribution::MultivariateNormal` caches the Cholesky
        // factor and samples `L * z + mu`. It is the closest Rust equivalent
        // of `Distributions.jl`'s `MvNormal`.
        let setup_start = Instant::now();
        let distribution =
            StatrsMvNormal::new(mean, covariance_flat).map_err(|error| error.to_string())?;
        setup_sec = setup_start.elapsed().as_secs_f64();

        let mut rng = StdRng::seed_from_u64(COMPARISON_NORMAL_SEED);
        for _ in 0..config.repeats {
            let start = Instant::now();
            for _ in 0..config.samples {
                let sample = distribution.sample(&mut rng);
                checksum += sample.iter().sum::<f64>();
            }
            durations.push(start.elapsed().as_secs_f64());
        }
    } else {
        let covariance = (0..dimension)
            .map(|row| {
                (0..dimension)
                    .map(|column| 2.0_f64.powi(-2 * distance(row, column)))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        let setup_start = Instant::now();
        let distribution = MvNormal::new(mean, covariance).map_err(|error| error.to_string())?;
        setup_sec = setup_start.elapsed().as_secs_f64();

        match config.normal_algorithm {
            NormalChoice::Polar | NormalChoice::Ziggurat => {
                let algorithm = if config.normal_algorithm == NormalChoice::Ziggurat {
                    NormalAlgorithm::Ziggurat
                } else {
                    NormalAlgorithm::MarsagliaPolar
                };
                let mut rng = StandardRng::with_algorithm(COMPARISON_NORMAL_SEED, algorithm);
                for _ in 0..config.repeats {
                    let start = Instant::now();
                    for _ in 0..config.samples {
                        distribution
                            .sample_inplace(&mut rng, &mut output)
                            .map_err(|error| error.to_string())?;
                        checksum += output.iter().sum::<f64>();
                    }
                    durations.push(start.elapsed().as_secs_f64());
                }
            }
            NormalChoice::RandDistr => {
                // `rand_distr::StandardNormal` uses the ZIGNOR Ziggurat method.
                // It owns its own random stream, so it is reported separately
                // from the shared xorshift64 comparison paths.
                let mut rng = StdRng::seed_from_u64(COMPARISON_NORMAL_SEED);
                let normal = StandardNormal;
                for _ in 0..config.repeats {
                    let start = Instant::now();
                    for _ in 0..config.samples {
                        distribution
                            .sample_inplace_with(|| normal.sample(&mut rng), &mut output)
                            .map_err(|error| error.to_string())?;
                        checksum += output.iter().sum::<f64>();
                    }
                    durations.push(start.elapsed().as_secs_f64());
                }
            }
            NormalChoice::Statrs => unreachable!("statrs is handled before constructing MvNormal"),
        }
    }

    let avg_sample_sec = durations.iter().sum::<f64>() / config.repeats as f64;
    let min_sample_sec = durations.iter().copied().fold(f64::INFINITY, f64::min);
    println!(
        "rust-{},{},{},{},{:.9},{:.9},{:.9},{:.12}",
        config.normal_algorithm.label(),
        config.dimension,
        config.samples,
        config.repeats,
        setup_sec,
        avg_sample_sec,
        min_sample_sec,
        checksum
    );
    Ok(())
}

fn main() {
    match parse_args().and_then(run) {
        Ok(()) => {}
        Err(error) => {
            eprintln!("error: {error}");
            process::exit(2);
        }
    }
}
