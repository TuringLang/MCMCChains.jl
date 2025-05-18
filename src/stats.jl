#################### Posterior Statistics ####################

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
    append_chains = true,
    demean::Bool = true,
    lags::AbstractVector{<:Integer} = _default_lags(chains, append_chains),
    kwargs...
)
    funs = Function[]
    func_names = @. Symbol("lag ", lags)
    for i in lags
        push!(funs, x -> autocor(x, [i], demean=demean)[1])
    end

    return summarize(
        chains, funs...;
        func_names = func_names,
        append_chains = append_chains,
        name = "Autocorrelation",
        kwargs...
    )
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
    cor(chains[; sections, append_chains = true, kwargs...])

Compute the Pearson correlation matrix for the chain.

Setting `append_chains=false` will return a vector of dataframes containing a correlation
matrix for each chain.
"""
function cor(
    chains::Chains;
    sections = _default_sections(chains),
    append_chains = true,
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Obstain names of parameters.
    names_of_params = names(_chains)

    if append_chains
        df = chaindataframe_cor("Correlation", names_of_params, to_matrix(_chains))
        return df
    else
        vector_of_df = [
            chaindataframe_cor(
                "Correlation - Chain $i", names_of_params, data
            )
            for (i, data) in enumerate(to_vector_of_matrices(_chains))
        ]
        return vector_of_df
    end
end

function chaindataframe_cor(name, names_of_params, chains::AbstractMatrix; kwargs...)
    # Compute the correlation matrix.
    cormat = cor(chains)

    # Summarize the results in a named tuple.
    nt = (; parameters = names_of_params,
          zip(names_of_params, (cormat[:, i] for i in axes(cormat, 2)))...)

    # Create a ChainDataFrame.
    return ChainDataFrame(name, nt; kwargs...)
end

"""
    changerate(chains[; sections, append_chains = true, kwargs...])

Compute the change rate for the chain.

Setting `append_chains=false` will return a vector of dataframes containing the change
rates for each chain.
"""
function changerate(
    chains::Chains{<:Real};
    sections = _default_sections(chains),
    append_chains = true,
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Obstain names of parameters.
    names_of_params = names(_chains)

    if append_chains
        df = chaindataframe_changerate("Change Rate", names_of_params, _chains.value.data)
        return df
    else
        vector_of_df = [
            chaindataframe_changerate(
                "Change Rate - Chain $i", names_of_params, data
            )
            for (i, data) in enumerate(to_vector_of_matrices(_chains))
        ]
        return vector_of_df
    end
end

function chaindataframe_changerate(name, names_of_params, chains; kwargs...)
    # Compute the change rates.
    changerates, mvchangerate = changerate(chains)

    # Summarize the results in a named tuple.
    nt = (; zip(names_of_params, changerates)..., multivariate = mvchangerate)

    # Create a ChainDataFrame.
    return ChainDataFrame(name, nt; kwargs...)
end

changerate(chains::AbstractMatrix{<:Real}) = changerate(reshape(chains, Val(3)))
function changerate(chains::AbstractArray{<:Real,3})
    niters, nparams, nchains = size(chains)

    changerates = zeros(nparams)
    mvchangerate = 0.0

    for chain in 1:nchains, iter in 2:niters
        isanychanged = false

        for param in 1:nparams
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
             etype = :bm,
             kwargs...])
Print chain metadata, summary statistics, and quantiles. Use `describe(chains)` for REPL output to `stdout`, or specify `io` for other streams (e.g., file output).
"""
function DataAPI.describe(
    io::IO,
    chains::Chains;
    q = [0.025, 0.25, 0.5, 0.75, 0.975],
    etype = :bm,
    kwargs...
)
    print(io, "Chains ", chains, ":\n\n", header(chains))
    
    summstats = summarystats(chains; etype = etype, kwargs...)
    println(io)
    show(io, MIME("text/plain"), summstats)
    
    qs = quantile(chains; q = q, kwargs...)
    println(io)
    show(io, MIME("text/plain"), qs)
end

# Convenience method for default IO
DataAPI.describe(chains::Chains; kwargs...) = DataAPI.describe(stdout, chains; kwargs...)

function _hpd(x::AbstractVector{<:Real}; alpha::Real=0.05)
    n = length(x)
    m = max(1, ceil(Int, alpha * n))

    y = sort(x)
    a = y[1:m]
    b = y[(n - m + 1):n]
    _, i = findmin(b - a)

    return [a[i], b[i]]
end

"""
    hpd(chn::Chains; alpha::Real=0.05, kwargs...)

Return the highest posterior density interval representing `1-alpha` probability mass.

Note that this will return a single interval and will not return multiple intervals for discontinuous regions.

# Examples

```julia-repl
julia> val = rand(500, 2, 3);
julia> chn = Chains(val, [:a, :b]);

julia> hpd(chn)
HPD
  parameters     lower     upper 
      Symbol   Float64   Float64 

           a    0.0554    0.9944
           b    0.0114    0.9460
```
"""
function hpd(chn::Chains; alpha::Real=0.05, kwargs...)
    labels = [:lower, :upper]
    l(x) = _hpd(x, alpha=alpha)[1]
    u(x) = _hpd(x, alpha=alpha)[2]
    return summarize(chn, l, u; name = "HPD", func_names = labels, kwargs...)
end

"""
    quantile(chains[; q = [0.025, 0.25, 0.5, 0.75, 0.975], append_chains = true, kwargs...])

Compute the quantiles for each parameter in the chain.

Setting `append_chains=false` will return a vector of dataframes containing the quantiles
for each chain.
"""
function quantile(
    chains::Chains;
    q::AbstractVector = [0.025, 0.25, 0.5, 0.75, 0.975],
    append_chains = true,
    kwargs...
)
    # compute quantiles
    funs = Function[]
    func_names = @. Symbol(100 * q, :%)
    for i in q
        push!(funs, x -> quantile(cskip(x), i))
    end

    return summarize(
        chains, funs...;
        func_names = func_names,
        append_chains = append_chains,
        name = "Quantiles",
        kwargs...
    )
end


"""
    function summarystats(
        chains;
        sections = _default_sections(chains),
        append_chains= true,
        autocov_method::AbstractAutocovMethod = AutocovMethod(),
        maxlag = 250,
        kwargs...
    )

Compute the mean, standard deviation, Monte Carlo standard error, bulk- and tail- effective
sample size, and ``\\widehat{R}`` diagnostic for each parameter in the chain.

Setting `append_chains=false` will return a vector of dataframes containing the summary
statistics for each chain.

When estimating the effective sample size, autocorrelations are computed for at most `maxlag` lags.
"""
function summarystats(
    chains::Chains;
    sections = _default_sections(chains),
    append_chains::Bool = true,
    autocov_method::MCMCDiagnosticTools.AbstractAutocovMethod = AutocovMethod(),
    maxlag = 250,
    name = "Summary Statistics",
    kwargs...
)
    # Store everything.
    funs = [mean∘cskip, std∘cskip]
    func_names = [:mean, :std]

    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Calculate MCSE and ESS/R-hat separately.
    nt_additional = NamedTuple()
    try
        mcse_df = MCMCDiagnosticTools.mcse(
            _chains; sections = nothing, autocov_method = autocov_method, maxlag = maxlag,
        )
        nt_additional = merge(nt_additional, (; mcse=mcse_df.nt.mcse))
    catch e
        @warn "MCSE calculation failed: $e"
    end

    try
        ess_tail_df = MCMCDiagnosticTools.ess(
            _chains; sections = nothing, autocov_method = autocov_method, maxlag = maxlag, kind=:tail
        )
        nt_additional = merge(nt_additional, (ess_tail=ess_tail_df.nt.ess,))
    catch e
        @warn "Tail ESS calculation failed: $e"
    end

    try
        ess_rhat_rank_df = MCMCDiagnosticTools.ess_rhat(
            _chains; sections = nothing, autocov_method = autocov_method, maxlag = maxlag, kind=:rank
        )
        nt_ess_rhat_rank = (
            ess_bulk=ess_rhat_rank_df.nt.ess,
            rhat=ess_rhat_rank_df.nt.rhat,
            ess_per_sec=ess_rhat_rank_df.nt.ess_per_sec
        )
        nt_additional = merge(nt_additional, nt_ess_rhat_rank)
    catch e
        @warn "Bulk ESS/R-hat calculation failed: $e"
    end

    # Possibly re-order the columns to stay backwards-compatible.
    additional_keys = (:mcse, :ess_bulk, :ess_tail, :rhat, :ess_per_sec)
    additional_df = ChainDataFrame("Additional", (; ((k, nt_additional[k]) for k in additional_keys if k ∈ keys(nt_additional))...))

    # Summarize.
    summary_df = summarize(
        _chains, funs...;
        func_names,
        append_chains,
        additional_df,
        name,
        sections = nothing
    )

    return summary_df
end

"""
    mean(chains[, params; kwargs...])

Calculate the mean of a chain.
"""
function mean(chains::Chains; kwargs...)
    # Store everything.
    funs = [mean∘cskip]
    func_names = [:mean]

    # Summarize.
    summary_df = summarize(
        chains, funs...;
        func_names = func_names,
        name = "Mean",
        kwargs...
    )

    return summary_df
end

mean(chn::Chains, syms) = mean(chn[:, syms, :])
# resolve method ambiguity with `mean(f, ::AbstractArray)`
mean(chn::Chains, syms::AbstractVector) = mean(chn[:, syms, :])