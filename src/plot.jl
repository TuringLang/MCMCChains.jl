export MeanPlot, AutocorPlot, HistogramPlot, DensityPlot, MixedDensityPlot, TracePlot
export plot

@userplot MeanPlot
@userplot AutocorPlot
@userplot HistogramPlot
@userplot DensityPlot
@userplot MixedDensityPlot
@userplot TracePlot

supportedplots = [
                  MeanPlot, 
                  AutocorPlot, 
                  HistogramPlot,
                  DensityPlot,
                  MixedDensityPlot,
                  TracePlot
]
colorparamters = [:chain, :parameter]

@recipe function f(p::MeanPlot; colordim = :chain)
    c, i = p.args
    seriestype := :line
    xaxis := "Iteration"
    yaxis := "Mean"
    
    @assert colordim ∈ colorparamters

    if colordim == :parameter
        title := "Chain $(c.chains[i])"
        labels := c.names
        c.range, MCMCChain.cummean(c.value[:, :, i])
    else
        title := c.names[i]
        labels := map(k -> "Chain $(c.chains[k])", 1:size(c)[3])
        c.range, MCMCChain.cummean(c.value[:, i, :])
    end
end

@recipe function f(p::AutocorPlot; maxlag = nothing, colordim = :chain)
    c, i = p.args
    lags = 0:(maxlag === nothing ? round(Int, 10 * log10(length(c.range))) : maxlag)
    ac = MCMCChain.autocor(c, lags=collect(lags))
    xaxis := "Lag"
    yaxis := "Autocorrelation"

    @assert colordim ∈ colorparamters
    if colordim == :parameter
        title := "Chain $(c.chains[i])"
        labels := c.names
        lags, ac.value[:, :, i]'
    else
        title := c.names[i]
        labels := map(k -> "Chain $(c.chains[k])", 1:size(c)[3])
        lags, ac.value[i, :, :]
    end
end

@recipe function f(p::HistogramPlot; colordim = :chain)
    c, i = p.args
    seriestype := :histogram
    xaxis := "Sample value"
    yaxis := "Frequency"
    fillalpha := 0.7

    @assert colordim ∈ colorparamters
    
    if colordim == :parameter
        title := "Chain $(c.chains[i])"
        labels := c.names
        c.value[:, :, i]
    else
        title := c.names[i]
        labels := map(k -> "Chain $(c.chains[k])", 1:size(c)[3])
        c.value[:, i, :]
    end
end

@recipe function f(p::DensityPlot; colordim = :chain)
    c, i = p.args
    seriestype := :density
    xaxis := "Sample value"
    yaxis := "Density"
    
    @assert colordim ∈ colorparamters
    
    if colordim == :parameter
        title := "Chain $(c.chains[i])"
        labels := c.names
        c.value[:, :, i]
    else
        title := c.names[i]
        labels := map(k -> "Chain $(c.chains[k])", 1:size(c)[3])
        c.value[:, i, :]
    end
end

@recipe function f(p::TracePlot; colordim = :chain)
    c, i = p.args
    seriestype := :line
    xaxis := "Iteration"
    yaxis := "Sample value"
    @assert colordim ∈ colorparamters
    
    if colordim == :parameter
        title := "Chain $(c.chains[i])"
        labels := c.names
        c.range, c.value[:, :, i]
    else
        title := c.names[i]
        labels := map(k -> "Chain $(c.chains[k])", 1:size(c)[3])
        c.range, c.value[:, i, :]
    end
end

@recipe function f(p::MixedDensityPlot; barbounds = (0, Inf), colordim = :chain)
    c, i = p.args
    discrete = MCMCChain.indiscretesupport(c, barbounds)
    
    seriestype := :density
    xaxis := "Sample value"
    yaxis := "Density"
    
    @assert colordim ∈ colorparamters
    
    if colordim == :parameter
        title := "Chain $(c.chains[i])"
        labels := c.names
        c.value[:, :, i]
    else
        seriestype := discrete[i] ? :histogram : :density
        yaxis := discrete[i] ? "Frequency" : "Density"
        fillalpha := discrete[i] ? 0.7 : 1.0
        title := c.names[i]
        labels := map(k -> "Chain $(c.chains[k])", 1:size(c)[3])
        c.value[:, i, :]
    end
end

@recipe function f(c::MCMCChain.AbstractChains;
                   ptypes = [TracePlot, MixedDensityPlot], 
                   width = 500, 
                   height = 250,
                   colordim = :chain
                  )

    # sanity checks
    @assert all(map(ptype -> ptype ∈ supportedplots, ptypes))
    @assert colordim ∈ colorparamters

    nrows, nvars, nchains = size(c.value)
    ntypes = length(ptypes)
    N = colordim == :chain ? nvars : nchains
    layout := (N, ntypes)
    size := (ntypes*width, N*height)
    indices = reshape(1:N*ntypes, ntypes, N)'

    legend --> false

    for (j, ptype) in enumerate(ptypes)
        for i in 1:N
            @series begin
                subplot := indices[i, j]
                colordim := colordim
                ptype([c, i])
            end
        end
    end
end

"""
    plot(c::AbstractChains, ptype::DataType)

Plot the summary of a MCMC simulation stored in `c` by using the plotting type specified
in ptype.

Optional Arguments:
- `width = 500`
- `height = 250`
- `colordim = :chain` ... either :chain or :parameter

"""
@recipe function f(c::AbstractChains,
              ptype::DataType;
              width = 500,
              height = 250,
              colordim = :chain
             )
    
    @assert ptype ∈ supportedplots
    @assert colordim ∈ colorparamters
    
    nrows, nvars, nchains = size(c.value)
    N = colordim == :chain ? nvars : nchains
    layout := (N, 1)
    size := (width, N*height)
    indices = reshape(1:N, 1, N)'

    legend --> false

    for i in 1:N
        @series begin
            subplot := indices[i, 1]
            colordim := colordim
            ptype([c, i])
        end
    end
end

# define alias functions for old syntax
translationdict = Dict(
                        :trace => TracePlot,
                        :mean => MeanPlot,
                        :density => DensityPlot,
                        :histogram => HistogramPlot,
                        :mixeddensity => MixedDensityPlot,
                        :autocor => AutocorPlot
                      )

function plot(c::AbstractChains, psyms::Vector{Symbol}; args...)
    @assert all(map(psym -> haskey(translationdict, psym), psyms))
    @warn "This syntax is deprecated, please use plot(c; ptypes = [TracePlot]) instead."
    return plot(c; ptypes = map(psym -> translationdict[psym], psyms))
end

function plot(c::AbstractChains, psym::Symbol; args...)
    @assert haskey(translationdict, psym)
    @warn "This syntax is deprecated, please use plot(c, $(translationdict[psym])) instead"
    return plot(c, translationdict[psym]; args...)
end
