using RecipesBase
using StatPlots
import Plots.plot
#################### Posterior Plot Recipies ####################

function plot(c::Chain.AbstractChains, ptype::Symbol; args...)
  #ptype == :autocor      ? autocorplot(c; legend=legend, args...) :
  #ptype == :bar          ? barplot(c; legend=legend, args...) :
  #ptype == :contour      ? contourplot(c; args...) :
  #ptype == :density      ? densityplot(c; legend=legend, args...) :
  ptype == :mean         ? meanplot(c; args...) :
  ptype == :mixeddensity ? mixeddensityplot(c; args...) :
  ptype == :trace        ? traceplot(c; args...) :
    throw(ArgumentError("unsupported plot type $ptype"))
end

struct MeanPlot end

@recipe function plot(ct::Tuple{Chain.AbstractChains, MeanPlot}; barbounds=(0, Inf))
    
    c, t = ct
    
    nrows, nvars, nchains = size(c.value)
    
    xaxis --> "Iteration"
    yaxis --> "Running Mean"
    layout := (nvars,1)

    x = repeat(collect(c.range), outer=[nchains])
    
    for i in 1:nvars
        @series begin
            seriestype := :line
            subplot := i
            title := c.names[i]
            labels := map(k -> "Chain $k", 1:nchains)
            Chain.cummean(c.value[:, i, :])
        end
    end
end

"""
    meanplot(c::Chain.Chains; args...)

Plot the running mean for each parameter.
"""
meanplot(c::Chain.AbstractChains; args...) = plot( (c, MeanPlot()); args...)

struct MixedDensityPlot end

@recipe function plot(ct::Tuple{Chain.AbstractChains, MixedDensityPlot}; barbounds=(0, Inf))
    
    c, t = ct
    
    nrows, nvars, nchains = size(c.value)
    
    xaxis --> "Value"
    layout := (nvars,1)

    x = repeat(collect(c.range), outer=[nchains])
    
    discrete = Chain.indiscretesupport(c, barbounds)
    for i in 1:nvars
        @series begin
            seriestype := discrete[i] ? :histogram : :density
            subplot := i
            title := c.names[i]
            labels := map(k -> "Chain $k", 1:nchains)
            c.value[:, i, :]
        end
    end
end

"""
    mixeddensityplot(c::Chain.Chains; args...)

Plot a density plot / histogram of each parameter.
"""
mixeddensityplot(c::Chain.AbstractChains; args...) = plot( (c, MixedDensityPlot()); args...)

struct TracePlot end

@recipe function plot(ct::Tuple{Chain.AbstractChains, TracePlot})
    
    c, t = ct
    
    nrows, nvars, nchains = size(c.value)
    
    xaxis --> "Iteration"
    yaxis --> "Value"
    layout := (nvars,1)

    x = repeat(collect(c.range), outer=[nchains])

    for i in 1:nvars
        @series begin
            seriestype := :line
            subplot := i
            title := c.names[i]
            labels := map(k -> "Chain $k", 1:nchains)
            c.value[:, i, :]
        end
    end
end

"""
    traceplot(c::Chain.Chains; args...)

Plot the trace of each parameter.
"""
traceplot(c::Chain.AbstractChains; args...) = plot( (c, TracePlot()); args...)