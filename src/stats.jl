#################### Posterior Statistics ####################

const DEFAULT_CI_PROB = 0.94

"""
    autocor(
        chains;
        append_chains = true,
        demean = true,
        [lags,]
        kwargs...,
    )

Compute the autocorrelation of each parameter for the chain.

The default `lags` are `[1, 5, 10, 50]`, upper-bounded by `n - 1` where `n` is the number of samples used in the estimation.

Setting `append_chains=false` will return a vector of dataframes containing the autocorrelations for each chain.
"""
function autocor(
    chains::Chains;
    sections = _default_sections(chains),
    append_chains::Bool = true,
    demean::Bool = true,
    lags::AbstractVector{<:Integer} = _default_lags(chains, append_chains),
    var_names = nothing,
    kwargs...,
)
    chn = Chains(chains, _clean_sections(chains, sections))

    # Obtain names of parameters.
    names_of_params = var_names === nothing ? names(chn) : var_names

    # Construct column names for lags.
    col_names = Symbol.("lag", lags)

    # avoids using summarize directly to support simultaneously computing a large number of
    # lags without constructing a huge NamedTuple
    if append_chains
        data = _permutedims_diagnostics(chn.value.data)
        vals = stack(map(eachslice(data; dims = 3)) do x
            return autocor(vec(x), lags; demean = demean)
        end)
        table = Tables.table(vals'; header = col_names)
        return SummaryStats(table; name = "Autocorrelation", labels = names_of_params)
    else
        data = to_vector_of_matrices(chn)
        return map(enumerate(data)) do (i, x)
            name_chain = "Autocorrelation (Chain $i)"
            vals = stack(map(eachslice(x; dims = 2)) do xi
                return autocor(xi, lags; demean = demean)
            end)
            table = Tables.table(vals'; header = col_names)
            return SummaryStats(table; name = name_chain, labels = names_of_params)
        end
    end
end

"""
    _default_lags(chains::Chains, append_chains::Bool)

Compute the vector of default lags for estimating the autocorrelation of the samples in `chains`.

The default lags are `[1, 5, 10, 50]`, upper-bounded by `n - 1` where `n` is the number of samples used in the estimation.
I.e., `n = size(chains, 1)` if `append_chains = false`, and `n = size(chains, 1) * size(chains, 3)` otherwise.
"""
function _default_lags(chains::Chains, append_chains::Bool)
    # Number of samples used for estimating the autocorrelation
    n = append_chains ? size(chains, 1) * size(chains, 3) : size(chains, 1)

    return [lag for lag in (1, 5, 10, 50) if lag < n]
end

"""
    cor(chains[; sections, append_chains = true])

Compute the Pearson correlation matrix for the chain.

Setting `append_chains=false` will return a vector of dataframes containing a correlation
matrix for each chain.
"""
function cor(
    chains::Chains;
    sections = _default_sections(chains),
    append_chains = true,
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Obstain names of parameters.
    names_of_params = names(_chains)

    if append_chains
        df = summarystats_cor("Correlation", names_of_params, to_matrix(_chains))
        return df
    else
        vector_of_df = [
            summarystats_cor("Correlation - Chain $i", names_of_params, data) for
            (i, data) in enumerate(to_vector_of_matrices(_chains))
        ]
        return vector_of_df
    end
end

function summarystats_cor(name, names_of_params, chains::AbstractMatrix)
    # Compute the correlation matrix.
    cormat = cor(chains)

    # Summarize the results in a dict
    dict = OrderedCollections.OrderedDict(
        k => v for (k, v) in zip(names_of_params, eachcol(cormat))
    )

    # Create a SummaryStats.
    return SummaryStats(dict; name = name, labels = names_of_params)
end

"""
    changerate(chains[; sections, append_chains = true])

Compute the change rate for the chain.

Setting `append_chains=false` will return a vector of dataframes containing the change
rates for each chain.
"""
function changerate(
    chains::Chains{<:Real};
    sections = _default_sections(chains),
    append_chains = true,
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Obstain names of parameters.
    names_of_params = names(_chains)

    if append_chains
        stats = summarystats_changerate("Change Rate", names_of_params, _chains.value.data)
        return stats
    else
        vector_of_stats = [
            summarystats_changerate("Change Rate - Chain $i", names_of_params, data) for
            (i, data) in enumerate(to_vector_of_matrices(_chains))
        ]
        return vector_of_stats
    end
end

function summarystats_changerate(name, names_of_params, chains)
    # Compute the change rates.
    changerates, mvchangerate = changerate(chains)

    # Summarize the results in a named tuple.
    nt = (; label = names_of_params, changerate = changerates)

    # Create a SummaryStats.
    return SummaryStats(nt; name = name), mvchangerate
end

changerate(chains::AbstractMatrix{<:Real}) = changerate(reshape(chains, Val(3)))
function changerate(chains::AbstractArray{<:Real,3})
    niters, nparams, nchains = size(chains)

    changerates = zeros(nparams)
    mvchangerate = 0.0

    for chain = 1:nchains, iter = 2:niters
        isanychanged = false

        for param = 1:nparams
            # update if the sample is different from the one in the previous iteration
            if chains[iter-1, param, chain] != chains[iter, param, chain]
                changerates[param] += 1
                isanychanged = true
            end
        end

        mvchangerate += isanychanged
    end

    factor = nchains * (niters - 1)
    changerates ./= factor
    mvchangerate /= factor

    changerates, mvchangerate
end

"""
    describe(io, chains[;
             q = [0.025, 0.25, 0.5, 0.75, 0.975],
             kwargs...])
Print chain metadata, summary statistics, and quantiles. Use `describe(chains)` for REPL output to `stdout`, or specify `io` for other streams (e.g., file output).
"""
function DataAPI.describe(
    io::IO,
    chains::Chains;
    q = [0.025, 0.25, 0.5, 0.75, 0.975],
    kwargs...,
)
    print(io, "Chains ", chains, ":\n\n", header(chains))

    summstats = summarystats(chains; kwargs...)
    println(io)
    show(io, MIME("text/plain"), summstats)

    qs = quantile(chains; q = q, kwargs...)
    println(io)
    println(io)
    show(io, MIME("text/plain"), qs)
end

# Convenience method for default IO
DataAPI.describe(chains::Chains; kwargs...) = DataAPI.describe(stdout, chains; kwargs...)

"""
    eti(chn::Chains; prob::Real=$DEFAULT_CI_PROB, kwargs...)

Return the equal-tailed interval (ETI) representing `prob` probability mass.

The bounds of the ETI are the symmetric quantiles so that the interval contains `prob`
probability mass.

Remaining keyword arguments are forwarded to [`summarize`](@ref).

See also [`quantile`](@ref), [`hdi`](@ref)

# Examples

```jldoctest
julia> using StableRNGs; rng = StableRNG(42);

julia> val = rand(rng, 500, 2, 3);

julia> chn = Chains(val, [:a, :b]);

julia> eti(chn)
ETI
    eti94           
 a  0.0337 .. 0.971
 b  0.0300 .. 0.962
```
"""
function PosteriorStats.eti(chn::Chains; prob::Real = DEFAULT_CI_PROB, kwargs...)
    eti_name = Symbol("eti$(_prob_to_string(prob))")
    return summarize(chn, eti_name => (x -> eti(x; prob)); name = "ETI", kwargs...)
end

"""
    hdi(chn::Chains; prob::Real=$DEFAULT_CI_PROB, method=:unimodal, kwargs...)

Return the highest density interval (HDI) representing `prob` probability mass.

Note that for the default (`method=:unimodal`), this will return a single interval.
For multiple intervals for discontinuous regions, use `method=:multimodal`.
See [`PosteriorStats.hdi`](@extref) for more details.

Remaining keyword arguments are forwarded to [`summarize`](@ref).

See also [`eti`](@ref)

# Examples

```jldoctest
julia> using StableRNGs; rng = StableRNG(42);

julia> val = rand(rng, 500, 2, 3);

julia> chn = Chains(val, [:a, :b]);

julia> hdi(chn)
HDI
    hdi94           
 a  0.0630 .. 0.994
 b  0.0404 .. 0.968
```
"""
function PosteriorStats.hdi(chn::Chains; prob::Real = DEFAULT_CI_PROB, kwargs...)
    hdi_name = Symbol("hdi$(_prob_to_string(prob))")
    return summarize(chn, hdi_name => (x -> hdi(x; prob)); name = "HDI", kwargs...)
end

_prob_to_string(prob; digits = 2) =
    replace(string(round(100 * prob; digits)), r"\.0+$" => "")

@deprecate hpd(chn::Chains; alpha::Real = 0.05, kwargs...) hdi(
    chn;
    prob = 1 - alpha,
    kwargs...,
)

"""
    quantile(chains[; q = (0.025, 0.25, 0.5, 0.75, 0.975), append_chains = true, kwargs...])

Compute the quantiles for each parameter in the chain.

Setting `append_chains=false` will return a vector of dataframes containing the quantiles
for each chain.

For intervals defined by symmetric quantiles, see [`eti`](@ref).
"""
function quantile(
    chains::Chains;
    q::Union{Tuple,AbstractVector} = (0.025, 0.25, 0.5, 0.75, 0.975),
    kwargs...,
)
    # compute quantiles
    func_names = Tuple(Symbol.(100 .* q, :%))
    return summarize(
        chains,
        func_names => (Base.Fix2(quantile, q) ∘ cskip);
        name = "Quantiles",
        kwargs...,
    )
end


"""
    summarystats(chains; kwargs...)

Compute default summary statistics from the `chains`.

`kwargs` are forwarded to [`summarize`](@ref). To customize the summary statistics, see
`summarize`.
"""
function summarystats(chains::Chains; name = "Summary Statistics", kwargs...)
    return summarize(chains; name, kwargs...)
end

"""
    mean(chains[, params; kwargs...])

Calculate the mean of a chain.
"""
function mean(chains::Chains; kwargs...)
    return summarize(chains, :mean => mean ∘ cskip; name = "Mean", kwargs...)
end

mean(chn::Chains, syms) = mean(chn[:, syms, :])
# resolve method ambiguity with `mean(f, ::AbstractArray)`
mean(chn::Chains, syms::AbstractVector) = mean(chn[:, syms, :])
