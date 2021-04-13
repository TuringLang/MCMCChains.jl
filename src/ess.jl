# methods
abstract type AbstractESSMethod end

"""
    ESSMethod <: AbstractESSMethod

The `ESSMethod` uses a standard algorithm for estimating the
effective sample size of MCMC chains.

It is is based on the discussion by [^Vehtari2019] and uses the
biased estimator of the autocovariance, as discussed by [^Geyer1992].
In contrast to Geyer, the divisor `n - 1` is used in the estimation of
the autocovariance to obtain the unbiased estimator of the variance for lag 0.

[^Geyer1992]: Geyer, C. J. (1992). Practical Markov Chain Monte Carlo. Statistical Science, 473-483. <https://projecteuclid.org/euclid.ss/1177011137>.

[^Vehtari2019]: Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P. C. (2021). Rank-normalization, folding, and localization: An improved ``\\widehat {R}`` for assessing convergence of MCMC. Bayesian Analysis. <https://arxiv.org/pdf/1903.08008.pdf>.
"""
struct ESSMethod <: AbstractESSMethod end

"""
    FFTESSMethod <: AbstractESSMethod

The `FFTESSMethod` uses a standard algorithm for estimating
the effective sample size of MCMC chains.

It is is based on the discussion by [^Vehtari2019] and uses the
biased estimator of the autocovariance, as discussed by [^Geyer1992].
In contrast to Geyer, the divisor `n - 1` is used in the estimation of
the autocovariance to obtain the unbiased estimator of the variance for lag 0.

In contrast to [`ESSMethod`](@ref), this method uses fast Fourier transforms
(FFTs) for estimating the autocorrelation.
"""
struct FFTESSMethod <: AbstractESSMethod end

"""
    BDAESSMethod <: AbstractESSMethod

The `BDAESSMethod` uses a standard algorithm for estimating the effective sample size of
MCMC chains.

It is is based on the discussion by [^Vehtari2019] and uses the
variogram estimator of the autocorrelation function discussed in [^Gelman2013].

[^Gelman2013]: Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). Bayesian data analysis. CRC press.
"""
struct BDAESSMethod <: AbstractESSMethod end

# caches
struct ESSCache{T,S}
    samples::Matrix{T}
    chain_var::Vector{S}
end

struct FFTESSCache{T,S,C,P,I}
    samples::Matrix{T}
    chain_var::Vector{S}
    samples_cache::C
    plan::P
    invplan::I
end

mutable struct BDAESSCache{T,S,M}
    samples::Matrix{T}
    chain_var::Vector{S}
    mean_chain_var::M
end

function build_cache(::ESSMethod, samples::Matrix, var::Vector)
    # check arguments
    niter, nchains = size(samples)
    length(var) == nchains || throw(DimensionMismatch())

    return ESSCache(samples, var)
end

function build_cache(::FFTESSMethod, samples::Matrix, var::Vector)
    # check arguments
    niter, nchains = size(samples)
    length(var) == nchains || throw(DimensionMismatch())

    # create cache for FFT
    T = complex(eltype(samples))
    n = nextprod([2, 3], 2 * niter - 1)
    samples_cache = Matrix{T}(undef, n, nchains)

    # create plans of FFTs
    fft_plan = plan_fft!(samples_cache, 1)
    ifft_plan = plan_ifft!(samples_cache, 1)

    return FFTESSCache(samples, var, samples_cache, fft_plan, ifft_plan)
end

function build_cache(::BDAESSMethod, samples::Matrix, var::Vector)
    # check arguments
    nchains = size(samples, 2)
    length(var) == nchains || throw(DimensionMismatch())

    return BDAESSCache(samples, var, mean(var))
end

update!(cache::ESSCache) = nothing

function update!(cache::FFTESSCache)
    # copy samples and add zero padding
    samples = cache.samples
    samples_cache = cache.samples_cache
    niter, nchains = size(samples)
    n = size(samples_cache, 1)
    T = eltype(samples_cache)
    @inbounds for j in 1:nchains
        for i in 1:niter
            samples_cache[i, j] = samples[i, j]
        end
        for i in (niter + 1):n
            samples_cache[i, j] = zero(T)
        end
    end

    # compute unnormalized autocovariance
    cache.plan * samples_cache
    @. samples_cache = abs2(samples_cache)
    cache.invplan * samples_cache

    nothing
end

function update!(cache::BDAESSCache)
    # recompute mean of within-chain variances
    cache.mean_chain_var = mean(cache.chain_var)

    return
end

function mean_autocov(k::Int, cache::ESSCache)
    # check arguments
    samples = cache.samples
    niter, nchains = size(samples)
    0 ≤ k < niter || throw(ArgumentError("only lags ≥ 0 and < $niter are supported"))

    # compute mean of unnormalized autocovariance estimates
    firstrange = 1:(niter - k)
    lastrange = (k + 1):niter
    s = mean(1:nchains) do i
        if eltype(samples) isa BlasReal
            # call into BLAS if possible
            x = dot(samples, firstrange, samples, lastrange)
            firstrange = firstrange .+ niter
            lastrange = lastrange .+ niter
            return x
        else
            # otherwise use views
            return dot(view(samples, firstrange, i), view(samples, lastrange, i))
        end
    end

    # normalize autocovariance estimators by `niter - 1` instead
    # of `niter - k` to obtain
    # - unbiased estimators of the variance for lag 0
    # - biased but more stable estimators for all other lags as discussed by
    #   Geyer (1992)
    return s / (niter - 1)
end

function mean_autocov(k::Int, cache::FFTESSCache)
    # check arguments
    niter, nchains = size(cache.samples)
    0 ≤ k < niter || throw(ArgumentError("only lags ≥ 0 and < $niter are supported"))

    # compute mean autocovariance
    # we use biased but more stable estimators as discussed by Geyer (1992)
    samples_cache = cache.samples_cache
    chain_var = cache.chain_var
    return mean(1:nchains) do i
        real(samples_cache[k + 1, i]) / real(samples_cache[1, i]) * chain_var[i]
    end
end

function mean_autocov(k::Int, cache::BDAESSCache)
    # check arguments
    samples = cache.samples
    niter, nchains = size(samples)
    0 ≤ k < niter || throw(ArgumentError("only lags ≥ 0 and < $niter are supported"))

    # compute mean autocovariance
    n = niter - k
    idxs = 1:n
    s = mean(1:nchains) do j
        return sum(idxs) do i
            abs2(samples[i, j] - samples[k + i, j])
        end
    end

    return cache.mean_chain_var - s / (2 * n)
end

"""
    ess(chains::Chains; kwargs...)

Estimate the effective sample size and the potential scale reduction.
"""
function ess(
    chains::Chains;
    sections = _default_sections(chains),
    kwargs...
)
    # subset the chain
    _chains = Chains(chains, _clean_sections(chains, sections))

    # estimate the effective sample size and rhat
    ess, rhat = ess_rhat(_chains.value.data; kwargs...)

    # Calculate ESS/minute if available
    dur = wall_duration(chains)

    # convert to namedtuple
    nt = if dur === missing
        merge((parameters = names(_chains),), (ess = ess, rhat = rhat))
    else
        merge((parameters = names(_chains),), (ess = ess, rhat = rhat, ess_per_sec=ess/dur))
    end

    return ChainDataFrame("ESS", nt)
end

function ess_rhat(
    chains::AbstractArray{<:Union{Missing,Real},3};
    method::AbstractESSMethod = ESSMethod(),
    maxlag::Int = 250
)
    # compute size of matrices (each chain is split!)
    niter = size(chains, 1) ÷ 2
    nparams = size(chains, 2)
    nchains = 2 * size(chains, 3)
    ntotal = niter * nchains

    # do not compute estimates if there is only one sample or lag
    maxlag = min(maxlag, niter - 1)
    maxlag > 0 || return fill(missing, nparams), fill(missing, nparams)

    # define caches for mean and variance
    U = typeof(zero(eltype(chains)) / 1)
    T = promote_type(eltype(chains), typeof(zero(eltype(chains)) / 1))
    chain_mean = Array{T}(undef, 1, nchains)
    chain_var = Array{T}(undef, nchains)
    samples = Array{T}(undef, niter, nchains)

    # compute correction factor
    correctionfactor = (niter - 1) / niter

    # define cache for the computation of the autocorrelation
    esscache = build_cache(method, samples, chain_var)
 
    # define output arrays
    ess = Vector{T}(undef, nparams)
    rhat = Vector{T}(undef, nparams)

    # for each parameter
    for (i, chains_slice) in enumerate(eachslice(chains; dims = 2))
        # check that no values are missing
        if any(x -> x === missing, chains_slice)
            rhat[i] = missing
            ess[i] = missing
            continue
        end

        # split chains
        copyto_split!(samples, chains_slice)

        # calculate mean of chains
        mean!(chain_mean, samples)

        # calculate within-chain variance
        @inbounds for j in 1:nchains
            chain_var[j] = var(view(samples, :, j); mean = chain_mean[j], corrected = true)
        end
        W = mean(chain_var)

        # compute variance estimator var₊, which accounts for between-chain variance as well
        var₊ = correctionfactor * W + var(chain_mean; corrected = true)
        inv_var₊ = inv(var₊)

        # estimate the potential scale reduction
        rhat[i] = sqrt(var₊ / W)

        # center the data around 0
        samples .-= chain_mean

        # update cache
        update!(esscache)

        # compute the first two autocorrelation estimates
        # by combining autocorrelation (or rather autocovariance) estimates of each chain
        ρ_odd = 1 - inv_var₊ * (W - mean_autocov(1, esscache))
        ρ_even = one(ρ_odd) # estimate at lag 0 is known

        # sum correlation estimates
        pₜ = ρ_even + ρ_odd
        sum_pₜ = pₜ

        k = 2
        while k < maxlag
            # compute subsequent autocorrelation of all chains
            # by combining estimates of each chain
            ρ_even = 1 - inv_var₊ * (W - mean_autocov(k, esscache))
            ρ_odd = 1 - inv_var₊ * (W - mean_autocov(k + 1, esscache))

            # stop summation if p becomes non-positive
            Δ = ρ_even + ρ_odd
            Δ > zero(Δ) || break

            # generate a monotone sequence
            pₜ = min(Δ, pₜ)

            # update sum
            sum_pₜ += pₜ

            # update indices
            k += 2
        end

        # estimate the effective sample size
        τ = 2 * sum_pₜ - 1
        ess[i] = ntotal / τ
    end

    return ess, rhat
end

"""
	copyto_split!(out::AbstractMatrix, x::AbstractMatrix)

Copy the elements of matrix `x` to matrix `out`, in which each column of `x` is split.

If the number of rows in `x` is odd, the sample at index `(size(x, 1) + 1) / 2` is dropped.
"""
function copyto_split!(out::AbstractMatrix, x::AbstractMatrix)
    # check dimensions
    nrows_out, ncols_out = size(out)
    nrows_x, ncols_x = size(x)
    ncols_out == 2 * ncols_x ||
        throw(DimensionMismatch("the output matrix must have twice as many columns as the input matrix"))
    nrows_out == nrows_x ÷ 2 ||
        throw(DimensionMismatch("the output matrix must have half as many rows as as the input matrix"))

    jout = 0
    offset = iseven(nrows_x) ? nrows_out : nrows_out + 1
    @inbounds for j in 1:ncols_x
        jout += 1
        for i in 1:nrows_out
            out[i, jout] = x[i, j]
        end

        jout += 1
        ix = offset
        for i in 1:nrows_out
            ix += 1
            out[i, jout] = x[ix, j]
        end
    end
    
    out
end
