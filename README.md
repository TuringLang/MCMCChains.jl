# MCMCChain.jl
[![Build Status](https://travis-ci.org/TuringLang/MCMCChain.jl.svg?branch=master)](https://travis-ci.org/TuringLang/MCMCChain.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/1av8osv0099nqw8m/branch/master?svg=true)](https://ci.appveyor.com/project/trappmartin/mcmcchain-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/TuringLang/MCMCChain.jl/badge.svg?branch=master)](https://coveralls.io/github/TuringLang/MCMCChain.jl?branch=master)

Implementation of Julia types for summarizing MCMC simulations and utility functions for diagnostics and visualizations.

## Example
The following simple example illustrates how to use Chain to visually summarize a MCMC simulation:
```julia
using MCMCChain
using StatPlots

theme(:ggplot2);

# Define the experiment
n_iter = 500;
n_name = 3;
n_chain = 2;

# experiment results
val = randn(n_iter, n_name, n_chain) .+ [1, 2, 3]';
val = hcat(val, rand(1:2, n_iter, 1, n_chain));

# construct a Chains object
chn = Chains(val);

# visualize the MCMC simulation results
p1 = plot(chn)
p2 = plot(chn, colordim = :parameter)

# save to a png file
savefig(p1, "demo-plot-parameters.png")
savefig(p2, "demo-plot-chains.png")

```
This code results in the visualizations shown below. Note that the plot function takes the additional arguments described in the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package.

Summarize parameters | Summarize chains
:-------------------------:|:-------------------------:
`plot(chn; colordim = :chain)` | `plot(chn; colordim = :parameter)`
![p1](https://user-images.githubusercontent.com/7974003/45822242-f0009180-bce2-11e8-8fa0-a97c8732400f.png)  |  ![p2](https://user-images.githubusercontent.com/7974003/45822249-f131be80-bce2-11e8-8dd3-42db7d58abd9.png)



## Manual
### Chains type
```julia
# construction of a Chains object
Chains(iterations::Int, params::Int;
    start = 1, thin = 1, chains = 1,
    names = String[])

# construction of a Chains object using an
# iteration * params * chains
# array (values).
Chains(values::Array{T, 3};
    start = 1, thin = 1, chains = 1,
    names = String[])

# Indexing a Chains object
chn = Chains(...)
chn_param1 = chn[:,2,:] # returns a new Chains object for parameter 2
chn[:,2,:] = ... # set values for parameter 2
```

### Convergence Diagnostics functions
#### Discrete Diagnostic
Options for method are  `[:weiss, :hangartner, :DARBOOT, MCBOOT, :billinsgley, :billingsleyBOOT]`

```julia
discretediag(c::AbstractChains; frac=0.3, method=:weiss, nsim=1000)
```

#### Gelman, Rubin, and Brooks Diagnostics
```julia
gelmandiag(c::AbstractChains; alpha=0.05, mpsrf=false, transform=false)
```

#### Geweke Diagnostic
```julia
gewekediag(c::AbstractChains; first=0.1, last=0.5, etype=:imse)
```

#### Heidelberger and Welch Diagnostics
```julia
heideldiag(c::AbstractChains; alpha=0.05, eps=0.1, etype=:imse)
```

#### Raftery and Lewis Diagnostic
```julia
rafterydiag(c::AbstractChains; q=0.025, r=0.005, s=0.95, eps=0.001)
```

### Plotting
```julia
# construct a plot
plot(c::AbstractChains, seriestype = (:traceplot, :mixeddensity))
plot(c::AbstractChains; ptypes = [TracePlot, MixedDensityPlot]) # deprecated
plot(c::AbstractChains; [:trace, :mixeddensity]) # deprecated

# construct trace plots
traceplot(c::AbstractChains)
plot(c::AbstractChains, seriestype = :traceplot)
plot(c::AbstractChains, TracePlot) # deprecated
plot(c::AbstractChains, :trace) # deprecated

# construct running average plots
meanplot(c::AbstractChains)
plot(c::AbstractChains, seriestype = :meanplot)
plot(c::AbstractChains, MeanPlot) # deprecated
plot(c::AbstractChains, :mean) # deprecated

# construct density plots
density(c::AbstractChains)
plot(c::AbstractChains, seriestype = :density)
plot(c::AbstractChains, DensityPlot) # deprecated
plot(c::AbstractChains, :density) # deprecated

# construct histogram plots
histogram(c::AbstractChains)
plot(c::AbstractChains, seriestype = :histogram)
plot(c::AbstractChains, HistogramPlot) # deprecated
plot(c::AbstractChains, :histogram) # deprecated

# construct mixed density plots
mixeddensity(c::AbstractChains)
plot(c::AbstractChains, seriestype = :mixeddensity)
plot(c::AbstractChains, MixedDensityPlot) # deprecated
plot(c::AbstractChains, :mixeddensity) # deprecated

# construct autocorrelation plots
autocorplot(c::AbstractChains)
plot(c::AbstractChains, seriestype = :autocorplot)
plot(c::AbstractChains, AutocorPlot) # deprecated
plot(c::AbstractChains, :autocor) # deprecated
```

## License Notice
Note that this package heavily uses and adapts code from the Mamba.jl package licensed under MIT License, see License.md.
