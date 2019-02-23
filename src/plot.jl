using RecipesBase
import RecipesBase: plot

export plot,
    traceplot,
    meanplot,
    densityplot,
    histogramplot,
    mixeddensityplot,
    autcorplot

@shorthands meanplot
@shorthands autocorplot
@shorthands mixeddensity
@shorthands traceplot
@shorthands corner

struct _TracePlot; c; val; end
struct _MeanPlot; c; val;  end
struct _DensityPlot; c; val;  end
struct _HistogramPlot; c; val;  end
struct _AutocorPlot; lags; val;  end

# define alias functions for old syntax
const translationdict = Dict(
                        :traceplot => _TracePlot,
                        :meanplot => _MeanPlot,
                        :density => _DensityPlot,
                        :histogram => _HistogramPlot,
                        :autocorplot => _AutocorPlot
                      )

const supportedplots = push!(collect(keys(translationdict)), :mixeddensity, :corner)

@recipe f(c::AbstractChains, s::Symbol) = c, [s]

@recipe function f(c::AbstractChains, i::Int; colordim = :chain, barbounds = (0, Inf), maxlag = nothing, section = :parameters)
    st = get(plotattributes, :seriestype, :traceplot)

    if colordim == :parameter
        title := "Chain $(chains(c)[i])"
        labels := c
        val = c.value[:, :, i]
    elseif colordim == :chain
        title := names(c)[i]
        labels := map(k -> "Chain $(chains(c)[k])", 1:size(c)[3])
        val = c.value[:, i, :]
    else
        throw(ArgumentError("`colordim` must be one of `:chain` or `:parameter`"))
    end

    if st == :mixeddensity
        discrete = MCMCChains.indiscretesupport(c, barbounds)
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
        ac = MCMCChains.autocor(c, lags=collect(lags); showall=true).summaries[1]
        val = colordim == :parameter ? ac.value[:, :, i]' : ac.value[i, :, :]
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
    [collect(skipmissing(p.val[:,k])) for k in 1:size(p.val, 2)]
end

@recipe function f(p::_HistogramPlot)
    xaxis --> "Sample value"
    yaxis --> "Frequency"
    fillalpha --> 0.7
    [collect(skipmissing(p.val[:,k])) for k in 1:size(p.val, 2)]
end

@recipe function f(p::_MeanPlot)
    seriestype := :path
    xaxis --> "Iteration"
    yaxis --> "Mean"
    range(p.c), MCMCChains.cummean(p.val)
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

@recipe function f(chn::MCMCChains.AbstractChains, parameters::AbstractVector{Symbol};
        colordim = :chain, section = :parameters)
    c = Chains(chn, section)
    colordim != :chain && error("Symbol names are interpreted as parameter names, only compatible with `colordim = :chain`")
    ret = indexin(parameters, Symbol.(keys(c)))
    any(y -> y == nothing, ret) && error("Parameter not found")
    c, Int.(ret)
end

@recipe function f(chn::MCMCChains.AbstractChains,
                   parameters::AbstractVector{<:Integer} = Int[];
                   width = 500,
                   height = 250,
                   colordim = :chain,
                   section = :parameters
                  )
    c = isempty(parameters) ? Chains(chn, section; sorted=true) : sort(chn)
    ptypes = get(plotattributes, :seriestype, (:traceplot, :mixeddensity))
    ptypes = ptypes isa AbstractVector || ptypes isa Tuple ? ptypes : (ptypes,)
    @assert all(map(ptype -> ptype ∈ supportedplots, ptypes))
    ntypes = length(ptypes)
    nrows, nvars, nchains = size(c)
    isempty(parameters) && (parameters = colordim == :chain ? (1:nvars) : (1:nchains))
    N = length(parameters)

    if !(:corner ∈ ptypes)
        layout := (N, ntypes)
        size --> (ntypes*width, N*height)
        indices = reshape(1:N*ntypes, ntypes, N)'
        legend --> false
        for (j, ptype) in enumerate(ptypes)
            for (i, par) in enumerate(parameters)
                @series begin
                    subplot := indices[i, j]
                    colordim := colordim
                    seriestype := ptype
                    c, par
                end
            end
        end
    else
        length(ptypes) > 1 && error(":corner is not compatible with multiple seriestypes")
        Corner(c, Symbol.(keys(c)[parameters]))
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
    RecipesBase.recipetype(:cornerplot, reduce(hcat, corner.c[:,s,:] for s in corner.parameters))
end
