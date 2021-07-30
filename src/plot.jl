@shorthands meanplot
@shorthands autocorplot
@shorthands mixeddensity
@shorthands pooleddensity
@shorthands traceplot
@shorthands corner
@shorthands ridgeline

struct _TracePlot; c; val; end
struct _MeanPlot; c; val;  end
struct _DensityPlot; c; val;  end
struct _HistogramPlot; c; val;  end
struct _AutocorPlot; lags; val;  end
struct _RidgelinePlot; par; val; end

# define alias functions for old syntax
const translationdict = Dict(
                        :traceplot => _TracePlot,
                        :meanplot => _MeanPlot,
                        :density => _DensityPlot,
                        :histogram => _HistogramPlot,
                        :autocorplot => _AutocorPlot,
                        :pooleddensity => _DensityPlot,
                        :ridgeline => _RidgelinePlot
                      )

const supportedplots = push!(collect(keys(translationdict)), :mixeddensity, :corner)

@recipe f(c::Chains, s::Symbol) = c, [s]

@recipe function f(
    chains::Chains, i::Int;
    colordim = :chain,
    barbounds = (-Inf, Inf),
    maxlag = nothing,
    append_chains = false
)
    st = get(plotattributes, :seriestype, :traceplot)
    c = append_chains || st == :pooleddensity ? pool_chain(chains) : chains

    if colordim == :parameter
        title --> "Chain $(MCMCChains.chains(c)[i])"
        label --> string.(names(c))
        val = c.value[:, :, i]
    elseif colordim == :chain
        title --> string(names(c)[i])
        label --> map(x -> "Chain $x", MCMCChains.chains(c))
        val = c.value[:, i, :]
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
        ac = autocor(c; sections = nothing, lags = lags)
        ac_mat = convert(Array, ac)
        val = colordim == :parameter ? ac_mat[:, :, i]' : ac_mat[i, :, :]
        _AutocorPlot(lags, val)
    elseif st ∈ supportedplots
        translationdict[st](c, val)
    else
        range(c), val
    end
end

@recipe function f(p::_DensityPlot)
    xaxis --> "Sample value"
    yaxis --> "Density"
    trim --> true
    [collect(skipmissing(p.val[:,k])) for k in 1:size(p.val, 2)]
end

@recipe function f(p::_HistogramPlot)
    xaxis --> "Sample value"
    yaxis --> "Frequency"
    fillalpha --> 0.7
    bins --> 25
    trim --> true
    [collect(skipmissing(p.val[:,k])) for k in 1:size(p.val, 2)]
end

@recipe function f(p::_MeanPlot)
    seriestype := :path
    xaxis --> "Iteration"
    yaxis --> "Mean"
    range(p.c), cummean(p.val)
end

@recipe function f(p::_AutocorPlot)
    seriestype := :path
    xaxis --> "Lag"
    yaxis --> "Autocorrelation"
    p.lags, p.val
end

@recipe function f(p::_TracePlot)
    seriestype := :path
    xaxis --> "Iteration"
    yaxis --> "Sample value"
    range(p.c), p.val
end

@recipe function f(
    chains::Chains,
    parameters::AbstractVector{Symbol};
    colordim = :chain
)
    colordim != :chain &&
        error("Symbol names are interpreted as parameter names, only compatible with ",
              "`colordim = :chain`")

    ret = indexin(parameters, names(chains))
    any(y === nothing for y in ret) && error("Parameter not found")

    return chains, Int.(ret)
end

@recipe function f(
    chains::Chains,
    parameters::AbstractVector{<:Integer} = Int[];
    sections = _default_sections(chains),
    width = 500,
    height = 250,
    colordim = :chain,
    append_chains = false
)
    _chains = isempty(parameters) ? Chains(chains, _clean_sections(chains, sections)) : chains
    c = append_chains ? pool_chain(_chains) : _chains
    ptypes = get(plotattributes, :seriestype, (:traceplot, :mixeddensity))
    ptypes = ptypes isa Symbol ? (ptypes,) : ptypes
    @assert all(ptype -> ptype ∈ supportedplots, ptypes)
    ntypes = length(ptypes)
    nrows, nvars, nchains = size(c)
    isempty(parameters) && (parameters = colordim == :chain ? (1:nvars) : (1:nchains))
    N = length(parameters)

    if :corner ∉ ptypes
        size --> (ntypes*width, N*height)
        legend --> false

        multiple_plots = N * ntypes > 1
        if multiple_plots
            layout := (N, ntypes)
        end

        i = 0
        for par in parameters
            for ptype in ptypes
                i += 1

                @series begin
                    if multiple_plots
                        subplot := i
                    end
                    colordim := colordim
                    seriestype := ptype
                    c, par
                end
            end
        end
    else
        ntypes > 1 && error(":corner is not compatible with multiple seriestypes")
        Corner(c, names(c)[parameters])
    end
end

struct Corner
    c
    parameters
end

@recipe function f(corner::Corner)
    label --> permutedims(corner.parameters)
    compact --> true
    size --> (600, 600)
    ar = collect(Array(corner.c.value[:, corner.parameters,i]) for i in chains(corner.c))
    RecipesBase.recipetype(:cornerplot, vcat(ar...))
end

@recipe function f(
    chains::Chains,
    par_names = chains.name_map[:parameters],
    q = [0.1, 0.9];
    riser = 0.1,
    spacer = 0.4,
    barbounds = (-Inf, Inf),
    #,
    #show_mean = true,
    #show_median = true,
    #show_interval = true,
    fill_q = true
    )

    st = get(plotattributes, :seriestype, :traceplot)
    n_iter, n_par, n_chains = size(chains)
    chain_arr = Array(chains, append_chains = true)
    discrete = indiscretesupport(chains, barbounds)
    s_type = [discrete[i] ? :histogram : :path for i in 1:length(par_names)]
    if st == :ridgeline
        k = 0
        for i in 1:length(par_names)
            h = riser + spacer*(i-1)
            qs = quantile(chain_arr[:,i], q)
            #chain_qs = filter(x -> q[1] <= x <= q[2], chain_arr[:,i])
            k_density = KernelDensity.kde(chain_arr[:,i])
            par = (fill_q ? filter(x -> qs[1] <= x <= qs[2], k_density.x) : k_density.x)
            val = (fill_q ? (pdf(k_density, par) .+ h) : k_density.density .+ h)
            chain_med = median(chain_arr[:,i])
            chain_mean = mean(chain_arr[:,i])
            min = minimum(k_density.density .+ h)
            k = 1
            @series begin
                subplot := k
                seriestype := :path
                #yaxis --> "$(par_names[i])"
                fillrange := min
                #fill --> (0, 0.5)
                label := nothing
                yguide --> "$(par_names[i])"
                par, val
            end
            @series begin
                subplot := k
                seriestype := :path
                label := nothing
                k_density.x, k_density.density .+ h
            end
            @series begin
                subplot := k
                seriestype := :hline
                label := nothing
                linecolor --> "#BBBBBB"
                linewidth --> 1.1
                [h]
            end

            @series begin
                subplot := k
                seriestype := :path
                linecolor --> "dark red"
                linewidth --> 1
                label --> (i == 1 ? "Mean" : nothing)
                #ylims := (minimum(val), maximum(val))
                [chain_mean, chain_mean], [min, min + pdf(k_density, chain_mean)]
            end

            @series begin
                subplot := k
                seriestype := :path
                linecolor --> "black"
                linewidth --> 1
                label --> (i == 1 ? "Median" : nothing)
                #ylims := (minimum(val), maximum(val))
                [chain_med, chain_med], [min, min + pdf(k_density, chain_med)]
            end

        end

    else

    end
end

