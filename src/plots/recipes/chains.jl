@recipe f(c::Chains, s::Symbol) = c, [s]

@recipe function f(
    chains::Chains,
    i::Int;
    colordim = :chain,
    barbounds = (-Inf, Inf),
    maxlag = nothing,
    append_chains = false,
)
    st = get(plotattributes, :seriestype, :traceplot)
    c = append_chains || st == :pooleddensity ? pool_chain(chains) : chains

    if colordim == :parameter
        title --> "Chain $(MCMCChains.chains(c)[i])"
        label --> permutedims(map(string, names(c)))
        val = c.value[:, :, i]
        _actual_indices_for_violin = 1:size(c, 2)
    elseif colordim == :chain
        title --> string(names(c)[i])
        label --> permutedims(map(x -> "Chain $x", MCMCChains.chains(c)))
        val = c.value[:, i, :]
        _actual_indices_for_violin = MCMCChains.chains(c)
    else
        throw(ArgumentError("`colordim` must be one of `:chain` or `:parameter`"))
    end

    if st == :mixeddensity || st == :pooleddensity
        discrete = indiscretesupport(c, barbounds)
        st = if colordim == :chain
            discrete[i] ? :histogram : :density
        else
            # NOTE: It might make sense to overlay histograms and density plots here.
            :density
        end
        seriestype := st
    end

    if st == :autocorplot
        lags = 0:(maxlag === nothing ? round(Int, 10 * log10(length(range(c)))) : maxlag)
        # Chains are already appended in `c` if desired, hence we use `append_chains=false`
        ac = autocor(c; sections = nothing, lags = lags, append_chains = false)
        ac_mat = convert(Array{Float64}, ac)
        val = colordim == :parameter ? ac_mat[:, :, i]' : ac_mat[i, :, :]
        _AutocorPlot(lags, val)
    elseif st ∈ (:violinplot, :violin)
        show_boxplot_kw = get(plotattributes, :show_boxplot, true)
        # Passes `_actual_indices_for_violin`` to `_ViolinPlot`` to ensure correct x-axis labels (parameter names or chain numbers).
        return _ViolinPlot(c, val, _actual_indices_for_violin, show_boxplot_kw, colordim)
    elseif st ∈ keys(translationdict)
        translationdict[st](c, val)
    else
        range(c), val
    end
end
