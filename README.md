# Chain.jl
[![Build Status](https://travis-ci.org/TuringLang/Chain.jl.svg?branch=master)](https://travis-ci.org/TuringLang/Chain.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/f9bs9jjpakyp1t59/branch/master?svg=true)](https://ci.appveyor.com/project/trappmartin/chain-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/TuringLang/Chain.jl/badge.svg?branch=master)](https://coveralls.io/github/TuringLang/Chain.jl?branch=master)


Implementation of Julia types for summarizing MCMC simulations and utility functions for diagnostics and visualizations. 

## Manual
####Plotting
```julia
# construct trace plots
plot(c::AbstractChain, :trace)
traceplot(c::AbstractChain)

# construct running average plots
plot(c::AbstractChain, :mean)
meanplot(c::AbstractChain)

# construct density plots
plot(c::AbstractChain, :density)
densityplot(c::AbstractChain)

# construct histogram plots
plot(c::AbstractChain, :histogram)
histogramplot(c::AbstractChain)

# construct mixed density plots
plot(c::AbstractChain, :mixeddensity)
mixeddensity plot(c::AbstractChain)

# construct autocorrelation plots
plot(c::AbstractChain, :autocor)
autocorplot(c::AbstractChain)

# combine different kinds of plots, e.g.
plot(c::AbstractChain, [:trace, :density])
```


## Example
The following simple example illustrates how to use Chain to visually summarize a MCMC simulation:
```julia
using Chain
using Plots

# Define the experiment
n_iter = 500
n_name = 3
n_chain = 2

# some sample experiment results
val = randn(n_iter, n_name, n_chain) .+ [1, 2, 3]'
val = hcat(val, rand(1:2, n_iter, 1, n_chain))

# construct a Chains object
chn = Chains(val)

# visualize a density plot / histogram plot, autocorrelation plot and a running average plot
plot(chn, [:mixeddensity, :autocor, :mean])

# save to a png file
savefig("demo-plot.png")
```
The above code results in the following visualization which can be adjusted using the additional arguments described in the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package.



## License Notice
Note that this package heavily uses and adapts code from the Mamba.jl package licensed under MIT License, see License.md.
