export MeanPlot, AutocorPlot, HistogramPlot, DensityPlot, MixedDensityPlot, TracePlot

@userplot MeanPlot
@userplot AutocorPlot
@userplot HistogramPlot
@userplot DensityPlot
@userplot MixedDensityPlot
@userplot TracePlot

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
    
    # value field chains x paramters x iterations
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
    
    # value field chains x paramters x iterations
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
    
    # value field chains x paramters x iterations
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
    
    # value field chains x paramters x iterations
    if colordim == :parameter
        title := "Chain $(c.chains[i])"
        labels := c.names
        c.value[:, :, i]
    else
        seriestype := discrete[i] ? :histogram : :density
        yaxis := discrete[i] ? "Frequency" : "Density"
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
    nrows, nvars, nchains = size(c.value)
    ntypes = length(ptypes)
    N = colordim == :chain ? nvars : nchains
    layout := (N, ntypes)
    size := (ntypes*width, N*height)
    indices = reshape(1:N*ntypes, ntypes, N)'

    @assert colordim ∈ colorparamters

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
