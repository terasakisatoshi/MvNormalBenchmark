use std::env;
use std::process;
use std::time::Instant;

use mvnormal_benchmark::{MvNormal, StandardRng};

struct Config {
    dimension: usize,
    samples: usize,
    repeats: usize,
}

fn usage() -> &'static str {
    "usage: mvnormal_benchmark --dim <dimension> --samples <samples> --repeats <repeats>"
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
    let mut arguments = env::args().skip(1);

    while let Some(flag) = arguments.next() {
        let value = arguments
            .next()
            .ok_or_else(|| format!("missing value for {flag}\n{}", usage()))?;
        match flag.as_str() {
            "--dim" => dimension = Some(parse_positive(&value, "--dim")?),
            "--samples" => samples = Some(parse_positive(&value, "--samples")?),
            "--repeats" => repeats = Some(parse_positive(&value, "--repeats")?),
            _ => return Err(format!("unknown argument {flag}\n{}", usage())),
        }
    }

    Ok(Config {
        dimension: dimension.ok_or_else(|| format!("missing --dim\n{}", usage()))?,
        samples: samples.ok_or_else(|| format!("missing --samples\n{}", usage()))?,
        repeats: repeats.ok_or_else(|| format!("missing --repeats\n{}", usage()))?,
    })
}

fn run(config: Config) -> Result<(), String> {
    let covariance = (0..config.dimension)
        .map(|row| {
            (0..config.dimension)
                .map(|column| if row == column { 2.0 } else { 0.1 })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mean = (0..config.dimension)
        .map(|index| (index as f64) * 0.01)
        .collect::<Vec<_>>();

    let setup_start = Instant::now();
    let distribution = MvNormal::new(mean, covariance).map_err(|error| error.to_string())?;
    let setup_sec = setup_start.elapsed().as_secs_f64();

    let mut rng = StandardRng::new(0x5EED_2021);
    let mut scratch = vec![0.0; config.dimension];
    let mut output = vec![0.0; config.dimension];
    let mut durations = Vec::with_capacity(config.repeats);
    let mut checksum = 0.0;

    for _ in 0..config.repeats {
        let start = Instant::now();
        for _ in 0..config.samples {
            distribution
                .sample_into(&mut rng, &mut scratch, &mut output)
                .map_err(|error| error.to_string())?;
            checksum += output.iter().sum::<f64>();
        }
        durations.push(start.elapsed().as_secs_f64());
    }

    let avg_sample_sec = durations.iter().sum::<f64>() / config.repeats as f64;
    let min_sample_sec = durations.iter().copied().fold(f64::INFINITY, f64::min);
    println!(
        "rust,{},{},{},{:.9},{:.9},{:.9},{:.12}",
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
