# methods
abstract type AbstractESSMethod end

"""
    ESSMethod <: AbstractESSMethod

The `ESSMethod` uses a standard algorithm for estimating the
effective sample size of MCMC chains.

It is is based on the discussion by
[Vehtari et al. (2019)](https://arxiv.org/pdf/1903.08008.pdf) and uses the
plug-in estimator of the autocorrelation function discussed by
[Geyer (1992)](https://projecteuclid.org/euclid.ss/1177011137).
"""
struct ESSMethod <: AbstractESSMethod end

"""
    FFTESSMethod <: AbstractESSMethod

The `FFTESSMethod` uses a standard algorithm for estimating
the effective sample size of MCMC chains.

It is is based on the discussion by
[Vehtari et al. (2019)](https://arxiv.org/pdf/1903.08008.pdf) and uses the
plug-in estimator of the autocorrelation function discussed by
[Geyer (1992)](https://projecteuclid.org/euclid.ss/1177011137). In contrast to
`ESSMethod`, it uses fast Fourier transforms (FFTs) for
estimating the autocorrelation function.
"""
struct FFTESSMethod <: AbstractESSMethod end

"""
    BDAESSMethod <: AbstractESSMethod

The `BDAESSMethod` uses a standard algorithm for estimating the effective sample size of
MCMC chains.

It is is based on the discussion by
[Vehtari et al. (2019)](https://arxiv.org/pdf/1903.08008.pdf) and uses the
estimator of the autocorrelation function discussed in 
[Bayesian Data Analysis (2013)](https://www.taylorfrancis.com/books/9780429113079).
"""
struct BDAESSMethod <: AbstractESSMethod end

# caches
mutable struct ESSCache{T}
    samples::Matrix{T}
    var::Vector{T}
end

mutable struct FFTESSCache{T,P,I}
    A::T
    plan::P
    invplan::I
    niter::Int
end

mutable struct BDAESSCache{T}
    samples::Matrix{T}
    var::Vector{T}
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
    A = Matrix{T}(undef, n, nchains)

    # create plans of FFTs
    fft_plan = plan_fft!(A, 1)
    ifft_plan = plan_ifft!(A, 1)

    return FFTESSCache(A, fft_plan, ifft_plan, niter)
end

function build_cache(::BDAESSMethod, samples::Matrix, var::Vector)
    # check arguments
    nchains = size(samples, 2)
    length(var) == nchains || throw(DimensionMismatch())

    return BDAESSCache(samples, var)
end

update_cache!(cache, samples::Matrix, var::Vector) = nothing

function update_cache!(cache::FFTESSCache, samples::Matrix, var::Vector)
    # check arguments
    niter, nchains = size(samples)
    niter == cache.niter || throw(DimensionMismatch())
    A = cache.A
    nchains == size(A, 2) || throw(DimensionMismatch())
	
    # copy samples and add zero padding
    n = size(A, 1)
    T = eltype(A)
    @inbounds for j in 1:nchains
        for i in 1:niter
            A[i, j] = samples[i, j]
        end
        for i in (niter + 1):n
            A[i, j] = zero(T)
        end
    end

    # compute unnormalized autocorrelation
    cache.plan * A
    @. A = abs2(A)
    cache.invplan * A

    nothing
end

# estimation of the autocorrelation function
function mean_autocorr(k::Int, cache::ESSCache)
    # check arguments
    samples = cache.samples
    niter, nchains = size(samples)
    0 ≤ k < niter || throw(ArgumentError("only lags ≥ 0 and < $niter are supported"))

    # compute mean autocorrelation
    var = cache.var
    s = zero(eltype(var))
    firstrange = 1:(niter - k)
    lastrange = (k + 1):niter
    @inbounds for i in 1:nchains
        # increment unnormalized correlation estimates
        if eltype(samples) isa LinearAlgebra.BlasReal
            # call into BLAS if possible
            s += dot(samples, firstrange, samples, lastrange) / var[i]
            firstrange = firstrange .+ niter
            lastrange = lastrange .+ niter
        else
            # otherwise use views
            s += dot(view(samples, firstrange, i), view(samples, lastrange, i)) / var[i]
        end
    end

    return s / length(samples)
end

function mean_autocorr(k::Int, cache::FFTESSCache)
    # check arguments
    niter = cache.niter
    0 ≤ k < niter || throw(ArgumentError("only lags ≥ 0 and < $niter are supported"))

    # compute mean autocorrelation
    A = cache.A
    nchains = size(A, 2)
    s = zero(real(eltype(A)))
    @inbounds for i in 1:nchains
        s += real(A[k + 1, i]) / real(A[1, i])
    end

    return s / nchains
end

function mean_autocorr(k::Int, cache::BDAESSCache)
    # check arguments
    samples = cache.samples
    niter, nchains = size(samples)
    0 ≤ k < niter || throw(ArgumentError("only lags ≥ 0 and < $niter are supported"))

    # compute mean autocorrelation
    var = cache.var
    s = zero(eltype(var))
    @inbounds for j in 1:nchains
        sj = zero(s)
        for i in 1:(niter - k)
            sj += abs2(samples[i, j] - samples[k + i, j])
        end
        s += sj / var[j]
    end

    return 1 - s / (2 * length(samples))
end

"""
    ess(chains::Chains; sections = _default_sections(chains), kwargs...)

Estimate the effective sample size and the potential scale reduction.
"""
function ess(
    chains::Chains;
    sections = _default_sections(chains), # perhaps directly nothing by default - or no sections argument?
    kwargs...
)
    # subset the chain
    _chains = Chains(chains, _clean_sections(chains, sections)) # may be unnecessary, already applied upstream

    # estimate the effective sample size and rhat
    ess, rhat = ess_rhat(_chains.value.data; kwargs...)

    # convert to namedtuple
    nt = merge((parameters = names(_chains),), (ess = ess, rhat = rhat))

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
    correctionfactor = niter / (niter - 1)

    # define cache for the computation of the autocorrelation
    esscache = build_cache(method, samples, chain_var)
 
    # define output arrays
    ess = Vector{T}(undef, nparams)
    rhat = Vector{T}(undef, nparams)

    # for each parameter
    for (i, chains_slice) in enumerate(eachcol(chains))
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
            chain_var[j] = var(view(samples, :, j); mean = chain_mean[j], corrected = false)
        end
        mean_chain_var = mean(chain_var)
        W = correctionfactor * mean_chain_var

        # compute variance estimator var₊, which accounts for between-chain variance as well
        var₊ = mean_chain_var + var(chain_mean; corrected = true)
        inv_var₊ = inv(var₊)

        # estimate the potential scale reduction
        rhat[i] = sqrt(var₊ / W)

        # center the data around 0
        samples .-= chain_mean

        # update cache
        update_cache!(esscache, samples, chain_var)

        # compute the first two autocorrelation terms
        mean_ρ = mean_autocorr(1, esscache)
        ρ_odd = 1 - inv_var₊ * (W - mean_ρ)
        ρ_even = one(ρ_odd)

        # sum correlation estimates
        pₜ = ρ_even + ρ_odd
        sum_pₜ = pₜ

        k = 2
        while k < maxlag
            # compute and combine autocorrelation of all chains
            mean_ρ = mean_autocorr(k, esscache)
            ρ_even = 1 - inv_var₊ * (W - mean_ρ)

            mean_ρ = mean_autocorr(k + 1, esscache)
            ρ_odd = 1 - inv_var₊ * (W - mean_ρ)

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
