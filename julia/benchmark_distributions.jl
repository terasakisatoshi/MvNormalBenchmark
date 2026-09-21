using Distributions
using LinearAlgebra
using Random

function parse_args(args)
    values = Dict{String, Int}(
        "dim" => 10,
        "samples" => 10_000,
        "repeats" => 5,
    )
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || throw(ArgumentError("unexpected argument: $arg"))

        key, value = if occursin('=', arg)
            split(arg[3:end], "=", limit=2) |> Tuple
        else
            i < length(args) || throw(ArgumentError("missing value for $arg"))
            (arg[3:end], args[i + 1])
        end
        haskey(values, key) || throw(ArgumentError("unknown option: --$key"))
        parsed = try
            parse(Int, value)
        catch
            throw(ArgumentError("--$key must be an integer"))
        end
        parsed >= 0 || throw(ArgumentError("--$key must be non-negative"))
        values[key] = parsed
        i += occursin('=', arg) ? 1 : 2
    end

    values["dim"] > 0 || throw(ArgumentError("--dim must be positive"))
    values["samples"] > 0 || throw(ArgumentError("--samples must be positive"))
    values["repeats"] > 0 || throw(ArgumentError("--repeats must be positive"))
    return values
end

function main(args)
    options = parse_args(args)
    dim = options["dim"]
    nsamples = options["samples"]
    repeats = options["repeats"]

    setup_rng = MersenneTwister(0x5eed)
    μ = randn(setup_rng, dim)
    A = randn(setup_rng, dim, dim)
    Σ = A * A' + dim * I

    setup_start = time_ns()
    d = MvNormal(μ, Σ)
    setup_sec = (time_ns() - setup_start) / 1.0e9

    # Compile the official package's sampling path before timing it.
    warmup_rng = MersenneTwister(0xabad1dea)
    out = Vector{Float64}(undef, dim)
    warmup_checksum = 0.0
    for _ in 1:nsamples
        rand!(warmup_rng, d, out)
        warmup_checksum += sum(out)
    end

    sample_times = Vector{Float64}(undef, repeats)
    checksum = 0.0
    rng = MersenneTwister(0xc0ffee)
    for repeat in 1:repeats
        start = time_ns()
        for _ in 1:nsamples
            rand!(rng, d, out)
            checksum += sum(out)
        end
        sample_times[repeat] = (time_ns() - start) / 1.0e9
    end

    # Keep the warmup result observable without including it in the reported checksum.
    checksum += 0.0 * warmup_checksum
    avg_sample_sec = sum(sample_times) / repeats
    min_sample_sec = minimum(sample_times)
    println("julia-distributions,$dim,$nsamples,$repeats,$setup_sec,$avg_sample_sec,$min_sample_sec,$checksum")
end

main(ARGS)
