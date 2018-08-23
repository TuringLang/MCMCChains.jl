# MCMCChain.jl
[![Build Status](https://travis-ci.org/TuringLang/Chain.jl.svg?branch=master)](https://travis-ci.org/TuringLang/Chain.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/f9bs9jjpakyp1t59/branch/master?svg=true)](https://ci.appveyor.com/project/trappmartin/chain-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/TuringLang/Chain.jl/badge.svg?branch=master)](https://coveralls.io/github/TuringLang/Chain.jl?branch=master)


Implementation of Julia types for summarizing MCMC simulations and utility functions for diagnostics and visualizations. 

## Example
The following simple example illustrates how to use Chain to visually summarize a MCMC simulation:
```julia
using MCMCChain
using Plots

# Define the experiment
n_iter = 500
n_name = 3
n_chain = 2

# experiment results
val = randn(n_iter, n_name, n_chain) .+ [1, 2, 3]'
val = hcat(val, rand(1:2, n_iter, 1, n_chain))

# construct a Chains object
chn = Chains(val)

# visualize a density plot / histogram plot, autocorrelation plot and a running average plot
plot(chn, [:mixeddensity, :autocor, :mean])

# save to a png file
savefig("demo-plot.png")
```
This code results in the visualization shown below. Note that the plot function takes the additional arguments described in the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package.

![demo_plot](https://user-images.githubusercontent.com/7974003/44415798-325e7380-a569-11e8-82e7-74acf7b1f359.png)

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
plot(c::AbstractChains, type::Symbol)

# construct trace plots
plot(c::AbstractChains, :trace)
traceplot(c::AbstractChains)

# construct running average plots
plot(c::AbstractChains, :mean)
meanplot(c::AbstractChains)

# construct density plots
plot(c::AbstractChains, :density)
densityplot(c::AbstractChains)

# construct histogram plots
plot(c::AbstractChains, :histogram)
histogramplot(c::AbstractChains)

# construct mixed density plots
plot(c::AbstractChains, :mixeddensity)
mixeddensityplot(c::AbstractChains)

# construct autocorrelation plots
plot(c::AbstractChains, :autocor)
autocorplot(c::AbstractChains)

# combine different kinds of plots, e.g.
plot(c::AbstractChain, [:trace, :density])
```

## License Notice
Note that this package heavily uses and adapts code from the Mamba.jl package licensed under MIT License, see License.md.
