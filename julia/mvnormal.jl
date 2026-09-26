using LinearAlgebra
using Random
import Statistics: mean

const ZIGGURAT_LAYERS = 256
const ZIGGURAT_R = 3.6541528853610088
const ZIGGURAT_V = 0.004928673233974658
const TWO_POW_63 = 9.223372036854775808e18
const SIGNED_MAGNITUDE_MASK = UInt64(0x7fffffffffffffff)

abstract type AbstractNormalRNG end

mutable struct NormalRNG <: AbstractNormalRNG
    state::UInt64
    algorithm::Symbol
    spare_normal::Float64
    has_spare::Bool
end

function NormalRNG(seed::Integer, algorithm::Symbol)
    algorithm in (:polar, :ziggurat) ||
        throw(ArgumentError("normal algorithm must be :polar or :ziggurat"))
    state = UInt64(seed)
    state == 0 && (state = UInt64(0x9e3779b97f4a7c15))
    return NormalRNG(state, algorithm, 0.0, false)
end

MarsagliaPolarRNG(seed::Integer) = NormalRNG(seed, :polar)
ZigguratRNG(seed::Integer) = NormalRNG(seed, :ziggurat)

@inline function next_u64!(rng::NormalRNG)
    value = rng.state
    value = xor(value, value << 13)
    value = xor(value, value >> 7)
    value = xor(value, value << 17)
    rng.state = value
    return value
end

@inline function uniform_open01(rng::NormalRNG)
    return (Float64(next_u64!(rng) >> 11) + 0.5) * (1.0 / 9007199254740992.0)
end

const ZIGGURAT_TABLES = let
    x = zeros(Float64, ZIGGURAT_LAYERS)
    k = zeros(UInt64, ZIGGURAT_LAYERS)
    w = zeros(Float64, ZIGGURAT_LAYERS)
    f = zeros(Float64, ZIGGURAT_LAYERS)
    tail_density = exp(-0.5 * ZIGGURAT_R * ZIGGURAT_R)
    q = ZIGGURAT_V / tail_density

    x[ZIGGURAT_LAYERS] = ZIGGURAT_R
    f[1] = 1.0
    f[ZIGGURAT_LAYERS] = tail_density
    w[1] = q / TWO_POW_63
    w[ZIGGURAT_LAYERS] = ZIGGURAT_R / TWO_POW_63
    k[1] = UInt64(floor(ZIGGURAT_R / q * TWO_POW_63))
    k[2] = UInt64(0)

    for index in (ZIGGURAT_LAYERS - 1):-1:2
        x[index] = sqrt(-2.0 * log(ZIGGURAT_V / x[index + 1] + f[index + 1]))
        f[index] = exp(-0.5 * x[index] * x[index])
        w[index] = x[index] / TWO_POW_63
        k[index + 1] = UInt64(floor(x[index] / x[index + 1] * TWO_POW_63))
    end
    (k, w, f)
end

const ZIGGURAT_K = ZIGGURAT_TABLES[1]
const ZIGGURAT_W = ZIGGURAT_TABLES[2]
const ZIGGURAT_F = ZIGGURAT_TABLES[3]

@inline function standard_normal_pair!(rng::NormalRNG)
    while true
        u = 2.0 * uniform_open01(rng) - 1.0
        v = 2.0 * uniform_open01(rng) - 1.0
        radius_squared = u * u + v * v
        if 0.0 < radius_squared < 1.0
            scale = sqrt(-2.0 * log(radius_squared) / radius_squared)
            return u * scale, v * scale
        end
    end
end

@inline function standard_normal_ziggurat!(rng::NormalRNG)
    while true
        bits = next_u64!(rng)
        index = Int(bits & UInt64(0xff)) + 1
        sign = (bits & (UInt64(1) << 63)) == 0 ? 1.0 : -1.0
        magnitude = bits & SIGNED_MAGNITUDE_MASK
        x = Float64(magnitude) * ZIGGURAT_W[index]
        if magnitude < ZIGGURAT_K[index]
            return sign * x
        end

        if index == 1
            while true
                tail_x = -log(uniform_open01(rng)) / ZIGGURAT_R
                tail_y = -log(uniform_open01(rng))
                if 2.0 * tail_y >= tail_x * tail_x
                    return sign * (ZIGGURAT_R + tail_x)
                end
            end
        end

        y = ZIGGURAT_F[index] +
            uniform_open01(rng) * (ZIGGURAT_F[index - 1] - ZIGGURAT_F[index])
        if y < exp(-0.5 * x * x)
            return sign * x
        end
    end
end

@inline function standard_normal!(rng::NormalRNG)
    if rng.algorithm === :ziggurat
        return standard_normal_ziggurat!(rng)
    end
    if rng.has_spare
        rng.has_spare = false
        return rng.spare_normal
    end
    first, second = standard_normal_pair!(rng)
    rng.spare_normal = second
    rng.has_spare = true
    return first
end

function fill_standard_normals!(rng::NormalRNG, output::AbstractVector{Float64})
    if rng.algorithm === :ziggurat
        @inbounds for index in eachindex(output)
            output[index] = standard_normal_ziggurat!(rng)
        end
        return output
    end

    index = firstindex(output)
    if rng.has_spare && !isempty(output)
        output[index] = standard_normal!(rng)
        index += 1
    end
    @inbounds while index + 1 <= lastindex(output)
        first, second = standard_normal_pair!(rng)
        output[index] = first
        output[index + 1] = second
        index += 2
    end
    if index <= lastindex(output)
        output[index] = standard_normal!(rng)
    end
    return output
end

const MvNormalRNG = Union{AbstractRNG, NormalRNG}

"""
    MvNormal(μ, Σ)

A multivariate normal distribution with mean vector `μ` and covariance
matrix `Σ`.  The covariance matrix is factorized once at construction time
as `Σ = L * L'`, where `L` is lower triangular.
"""
struct MvNormal
    μ::Vector{Float64}
    L::Matrix{Float64}

    # Keep the two-argument constructor reserved for validation below.
    MvNormal(μ::Vector{Float64}, L::Matrix{Float64}, ::Val{:validated}) = new(μ, L)
end

function MvNormal(μ::AbstractVector{<:Real}, Σ::AbstractMatrix{<:Real})
    n = length(μ)
    size(Σ, 1) == size(Σ, 2) ||
        throw(DimensionMismatch("covariance matrix must be square"))
    size(Σ, 1) == n ||
        throw(DimensionMismatch("mean and covariance dimensions must agree"))
    n > 0 || throw(ArgumentError("dimension must be positive"))

    μf = Float64.(μ)
    Σf = Matrix{Float64}(Σ)
    all(isfinite, μf) || throw(ArgumentError("mean must contain only finite values"))
    all(isfinite, Σf) ||
        throw(ArgumentError("covariance matrix must contain only finite values"))
    issymmetric(Σf) || throw(ArgumentError("covariance matrix must be symmetric"))

    factor = try
        cholesky(Symmetric(Σf); check=true)
    catch err
        throw(ArgumentError(
            "covariance matrix must be positive definite: $(sprint(showerror, err))"))
    end

    # Cholesky stores the factor in a dense buffer whose unused upper
    # triangle may still contain entries from the input matrix.  Keep only
    # the mathematical lower-triangular factor required by x = μ + L*z.
    L = tril(Matrix(factor.L))
    return MvNormal(μf, L, Val(:validated))
end

"""Construct a zero-mean multivariate normal distribution."""
MvNormal(Σ::AbstractMatrix{<:Real}) =
    MvNormal(zeros(Float64, size(Σ, 1)), Σ)

"""Return the dimensionality of the distribution."""
dimension(d::MvNormal) = length(d.μ)
Base.length(d::MvNormal) = dimension(d)

"""Return a copy of the mean vector."""
mean(d::MvNormal) = copy(d.μ)

"""Return the covariance matrix reconstructed from the Cholesky factor."""
covariance(d::MvNormal) = d.L * d.L'

"""Return a copy of the lower-triangular Cholesky factor."""
cholesky_factor(d::MvNormal) = copy(d.L)

"""
    sample!(rng, d, out)

Fill `out` with one sample from `d` and return `out`.  The calculation is
`out = μ + L * z`, with `z` a vector of independent standard normal values.
"""
function fill_standard_normals!(rng::AbstractRNG, output::AbstractVector{Float64})
    return randn!(rng, output)
end

function sample!(rng::MvNormalRNG, d::MvNormal, out::AbstractVector{Float64})
    length(out) == dimension(d) ||
        throw(DimensionMismatch("output vector dimension must agree with distribution"))

    # Generate z directly in the output buffer.  Process columns from right
    # to left so that out[j] still contains z[j] when column j is used.  The
    # entries at indices i >= j have already become partial outputs, so this
    # order permits contiguous column-major access without another buffer.
    fill_standard_normals!(rng, out)
    n = dimension(d)
    @inbounds for j in n:-1:1
        z_j = out[j]
        out[j] = d.μ[j]
        i = j
        while i + 1 <= n
            out[i] = muladd(d.L[i, j], z_j, out[i])
            out[i + 1] = muladd(d.L[i + 1, j], z_j, out[i + 1])
            i += 2
        end
        while i <= n
            out[i] = muladd(d.L[i, j], z_j, out[i])
            i += 1
        end
    end
    return out
end

"""Draw one sample using the supplied random-number generator."""
function sample(rng::MvNormalRNG, d::MvNormal)
    out = Vector{Float64}(undef, dimension(d))
    return sample!(rng, d, out)
end

"""Draw one sample using Julia's default random-number generator."""
sample(d::MvNormal) = sample(Random.default_rng(), d)

"""
    sample(rng, d, nsamples)

Draw `nsamples` samples.  Each column of the returned matrix is one sample.
"""
function sample(rng::MvNormalRNG, d::MvNormal, nsamples::Integer)
    nsamples >= 0 || throw(ArgumentError("number of samples must be non-negative"))
    out = Matrix{Float64}(undef, dimension(d), nsamples)
    return sample!(rng, d, out)
end

"""
    sample!(rng, d, out)

Fill the preallocated matrix `out` with one sample per column from `d`.
Standard normals are generated directly into `out`, the lower-triangular
Cholesky factor is applied in place with `lmul!`, and the mean is added.
Using `out` as the standard-normal buffer as well avoids a separate scratch
matrix and lets the in-place triangular multiply run slightly faster than an
out-of-place `mul!` into a second buffer.
"""
function sample!(rng::MvNormalRNG, d::MvNormal, out::StridedMatrix{Float64})
    n = dimension(d)
    size(out, 1) == n ||
        throw(DimensionMismatch("output row dimension must agree with distribution"))

    fill_standard_normals!(rng, vec(out))
    lmul!(LowerTriangular(d.L), out)
    @inbounds for column in axes(out, 2)
        for row in axes(out, 1)
            out[row, column] += d.μ[row]
        end
    end
    return out
end

sample(d::MvNormal, nsamples::Integer) =
    sample(Random.default_rng(), d, nsamples)
