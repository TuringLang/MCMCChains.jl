using RecipesBase
using StatPlots
import Plots.plot
#################### Posterior Plot Recipies ####################

@recipe function plot(c::MCMCChain.AbstractChains, ptypes::Vector{Symbol}; 
                      maxlag = round(Int, 10 * log10(length(c.range))),
  barbounds = (0, Inf))

  nrows, nvars, nchains = size(c.value)
  ntypes = length(ptypes)
  layout := (nvars, ntypes)
  indices = reshape(1:nvars*ntypes, ntypes, nvars)'
  lags = 0:maxlag

  legend --> false

  for (j, ptype) in enumerate(ptypes)

    if ptype == :mean
      for i in 1:nvars
        @series begin
          seriestype := :line
          xaxis := "Iteration"
          yaxis := "Mean"
          subplot := indices[i, j]
          title := c.names[i]
          labels := map(k -> "Chain $k", 1:nchains)
          MCMCChain.cummean(c.value[:, i, :])
        end
      end
    elseif ptype == :autocor
      ac = MCMCChain.autocor(c, lags=collect(lags))
      for i in 1:nvars
        @series begin
          seriestype := :line
          xaxis := "Lag"
          yaxis := "Autocorrelation"
          subplot := indices[i, j]
          title := c.names[i]
          labels := map(k -> "Chain $k", 1:nchains)
          ac.value[i, :, :]
        end
      end
    elseif ptype == :mixeddensity
      discrete = MCMCChain.indiscretesupport(c, barbounds)
      for i in 1:nvars
        @series begin
          seriestype := discrete[i] ? :histogram : :density
          xaxis := "Value"
          yaxis := "Density"
          subplot := indices[i, j]
          title := c.names[i]
          labels := map(k -> "Chain $k", 1:nchains)
          c.value[:, i, :]
        end
      end
    elseif ptype in [:density, :histogram, :trace]
      for i in 1:nvars
        @series begin
          seriestype := ptype == :trace ? :line : ptype
          xaxis := ptype == :trace ? "Iteration" : "Value"
          yaxis := ptype == :trace ? "Value" : "Density"
          subplot := indices[i, j]
          title := c.names[i]
          labels := map(k -> "Chain $k", 1:nchains)
          c.value[:, i, :]
        end
      end
    end
  end
end

plot(c::MCMCChain.AbstractChains, ptype::Symbol; args...) = plot(c, [ptype]; args...)

meanplot(c::MCMCChain.AbstractChains; args...) = plot(c, [:mean]; args...)
autocorplot(c::MCMCChain.AbstractChains; args...) = plot(c, [:autocor]; args...)
histogramplot(c::MCMCChain.AbstractChains; args...) = plot(c, [:histogram]; args...)
densityplot(c::MCMCChain.AbstractChains; args...) = plot(c, [:density]; args...)
mixeddensityplot(c::MCMCChain.AbstractChains; args...) = plot(c, [:mixeddensity]; args...)
traceplot(c::MCMCChain.AbstractChains; args...) = plot(c, [:trace]; args...)
