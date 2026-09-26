include(joinpath(@__DIR__, "mvnormal.jl"))

const COMPARISON_NORMAL_SEED = 0x5EED2021

function parse_args(args)
    values = Dict{String, Int}(
        "dim" => 10,
        "samples" => 10_000,
        "repeats" => 5,
    )
    normal = :julia
    batch = false
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || throw(ArgumentError("unexpected argument: $arg"))

        if arg == "--batch"
            batch = true
            i += 1
            continue
        end

        key, value = if occursin('=', arg)
            split(arg[3:end], "=", limit=2) |> Tuple
        else
            i < length(args) || throw(ArgumentError("missing value for $arg"))
            (arg[3:end], args[i + 1])
        end
        if key == "normal"
            value in ("julia", "polar", "ziggurat") ||
                throw(ArgumentError("--normal must be julia, polar, or ziggurat"))
            normal = Symbol(value)
            i += occursin('=', arg) ? 1 : 2
            continue
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
    return (values=values, normal=normal, batch=batch)
end

function main(args)
    # Keep the comparison single-threaded even though Julia's `mul!` would
    # otherwise use a multithreaded BLAS for the batched path.
    BLAS.set_num_threads(1)

    parsed = parse_args(args)
    options = parsed.values
    normal = parsed.normal
    batch = parsed.batch
    dim = options["dim"]
    nsamples = options["samples"]
    repeats = options["repeats"]

    normal === :julia && Random.seed!(0xc0ffee)
    μ = [0.01 * (index - 1) for index in 1:dim]
    Σ = [ldexp(1.0, -2 * abs(row - column)) for row in 1:dim, column in 1:dim]

    setup_start = time_ns()
    d = MvNormal(μ, Σ)
    setup_sec = (time_ns() - setup_start) / 1.0e9

    rng = if normal === :julia
        Random.default_rng()
    elseif normal === :polar
        MarsagliaPolarRNG(COMPARISON_NORMAL_SEED)
    else
        ZigguratRNG(COMPARISON_NORMAL_SEED)
    end
    warmup_rng = if normal === :julia
        Random.default_rng()
    elseif normal === :polar
        MarsagliaPolarRNG(0xabad1dea)
    else
        ZigguratRNG(0xabad1dea)
    end

    sample_times = Vector{Float64}(undef, repeats)
    checksum = 0.0

    if batch
        out = Matrix{Float64}(undef, dim, nsamples)

        # Compile the batched sampling path before timing it.
        sample!(warmup_rng, d, out)
        warmup_checksum = sum(out)

        for repeat in 1:repeats
            start = time_ns()
            sample!(rng, d, out)
            sample_times[repeat] = (time_ns() - start) / 1.0e9
            checksum += sum(out)
        end
    else
        # Compile the per-sample sampling path before timing it.
        out = Vector{Float64}(undef, dim)
        warmup_checksum = 0.0
        for _ in 1:nsamples
            sample!(warmup_rng, d, out)
            warmup_checksum += sum(out)
        end

        for repeat in 1:repeats
            start = time_ns()
            for _ in 1:nsamples
                sample!(rng, d, out)
                checksum += sum(out)
            end
            sample_times[repeat] = (time_ns() - start) / 1.0e9
        end
    end

    # Keep the warmup result observable without including it in the reported checksum.
    checksum += 0.0 * warmup_checksum
    avg_sample_sec = sum(sample_times) / repeats
    min_sample_sec = minimum(sample_times)
    label = (normal === :julia ? "julia" : "julia-$normal") * (batch ? "-batch" : "")
    println("$label,$dim,$nsamples,$repeats,$setup_sec,$avg_sample_sec,$min_sample_sec,$checksum")
end

main(ARGS)
