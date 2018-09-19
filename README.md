# MCMCChain.jl
[![Build Status](https://travis-ci.org/TuringLang/MCMCChain.jl.svg?branch=master)](https://travis-ci.org/TuringLang/MCMCChain.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/1av8osv0099nqw8m/branch/master?svg=true)](https://ci.appveyor.com/project/trappmartin/mcmcchain-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/TuringLang/MCMCChain.jl/badge.svg?branch=master)](https://coveralls.io/github/TuringLang/MCMCChain.jl?branch=master)

Implementation of Julia types for summarizing MCMC simulations and utility functions for diagnostics and visualizations. 

## Example
The following simple example illustrates how to use Chain to visually summarize a MCMC simulation:
```julia
using MCMCChain
using Plots, StatPlots

# Define the experiment
n_iter = 500
n_name = 3
n_chain = 2

# experiment results
val = randn(n_iter, n_name, n_chain) .+ [1, 2, 3]'
val = hcat(val, rand(1:2, n_iter, 1, n_chain))

# construct a Chains object
chn = Chains(val)

# visualize the MCMC simulation results 
plot(chn)

# save to a png file
savefig("demo-plot.png")
```
This code results in the visualization shown below. Note that the plot function takes the additional arguments described in the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package.

![demo_plot](https://user-images.githubusercontent.com/7974003/45752721-5798dd00-bc0e-11e8-817f-0634f8243c87.png)

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
plot(c::AbstractChains; ptypes = [TracePlot, MeanPlot])

# construct trace plots
plot(c::AbstractChains, ptypes = [TracePlot])
traceplot(c::AbstractChains, variable::Int)

# construct running average plots
plot(c::AbstractChains, ptypes = [MeanPlot])
meanplot(c::AbstractChains, variable::Int)

# construct density plots
plot(c::AbstractChains, ptypes = [DensityPlot])
densityplot(c::AbstractChains, variable::Int)

# construct histogram plots
plot(c::AbstractChains, ptypes = [HistogramPlot])
histogramplot(c::AbstractChains, variable::Int)

# construct mixed density plots
plot(c::AbstractChains, ptypes = [MixedDensityPlot])
mixeddensityplot(c::AbstractChains, variable::Int)

# construct autocorrelation plots
plot(c::AbstractChains, ptypes = [AutocorPlot])
autocorplot(c::AbstractChains, variable::Int)
```

## License Notice
Note that this package heavily uses and adapts code from the Mamba.jl package licensed under MIT License, see License.md.
