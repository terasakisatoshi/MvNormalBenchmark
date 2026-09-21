using LinearAlgebra
using Random
import Statistics: mean

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
function sample!(rng::AbstractRNG, d::MvNormal, out::AbstractVector{Float64})
    length(out) == dimension(d) ||
        throw(DimensionMismatch("output vector dimension must agree with distribution"))

    # Generate z directly in the output buffer.  Since L is lower
    # triangular, processing rows from bottom to top lets us overwrite
    # out[i] after all z[j] needed by that row (j < i) have been read.
    randn!(rng, out)
    @inbounds for i in dimension(d):-1:1
        value = d.μ[i]
        for j in 1:i
            value += d.L[i, j] * out[j]
        end
        out[i] = value
    end
    return out
end

"""Draw one sample using the supplied random-number generator."""
function sample(rng::AbstractRNG, d::MvNormal)
    out = Vector{Float64}(undef, dimension(d))
    return sample!(rng, d, out)
end

"""Draw one sample using Julia's default random-number generator."""
sample(d::MvNormal) = sample(Random.default_rng(), d)

"""
    sample(rng, d, nsamples)

Draw `nsamples` samples.  Each column of the returned matrix is one sample.
"""
function sample(rng::AbstractRNG, d::MvNormal, nsamples::Integer)
    nsamples >= 0 || throw(ArgumentError("number of samples must be non-negative"))
    z = randn(rng, dimension(d), nsamples)
    out = d.L * z
    out .+= d.μ
    return out
end

sample(d::MvNormal, nsamples::Integer) =
    sample(Random.default_rng(), d, nsamples)
