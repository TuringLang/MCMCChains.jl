export MeanPlot, AutocorPlot, HistogramPlot, DensityPlot, MixedDensityPlot, TracePlot

@userplot MeanPlot
@userplot AutocorPlot
@userplot HistogramPlot
@userplot DensityPlot
@userplot MixedDensityPlot
@userplot TracePlot

@recipe function f(p::MeanPlot)
  c, i = p.args
  seriestype := :line
  xaxis := "Iteration"
  yaxis := "Mean"
  MCMCChain.cummean(c.value[:, i, :])
end

@recipe function f(p::AutocorPlot; maxlag = nothing)
  c, i = p.args
  lags = 0:(maxlag === nothing ? round(Int, 10 * log10(length(c.range))) : maxlag)
  ac = MCMCChain.autocor(c, lags=collect(lags))
  xaxis := "Lag"
  yaxis := "Autocorrelation"
  ac.value[i, :, :]
end

@recipe function f(p::HistogramPlot)
  c, i = p.args
  seriestype := :histogram
  xaxis := "Value"
  yaxis := "Value"
  c.value[:, i, :]
end

@recipe function f(p::DensityPlot)
  c, i = p.args
  seriestype := :density
  xaxis := "Value"
  yaxis := "Density"
  c.value[:, i, :]
end

@recipe function f(p::TracePlot)
  c, i = p.args
  seriestype := :line
  xaxis := "Iteration"
  yaxis := "Density"
  c.value[:, i, :]
end

@recipe function f(p::MixedDensityPlot, barbounds = (0, Inf))
  c, i = p.args
  discrete = MCMCChain.indiscretesupport(c, barbounds)
  seriestype := discrete[i] ? :histogram : :density
  xaxis := "Value"
  yaxis := "Density"
  c.value[:, i, :]
end

@recipe function f(c::MCMCChain.AbstractChains;
                   ptypes = [TracePlot, MeanPlot], 
                   width = 500, 
                   height = 250
                  )
  nrows, nvars, nchains = size(c.value)
  ntypes = length(ptypes)
  layout := (nvars, ntypes)
  size := (ntypes*width, nvars*height)
  indices = reshape(1:nvars*ntypes, ntypes, nvars)'

  legend --> false

  for (j, ptype) in enumerate(ptypes)
    for i in 1:nvars
      @series begin
          subplot := indices[i, j]
          title := c.names[i]
          labels := map(k -> "Chain $k", 1:nchains)
          ptype([c, i])
      end
    end
  end
end
