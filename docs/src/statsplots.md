# StatsPlots.jl

MCMCChains implements many functions for plotting via [StatsPlots.jl](https://github.com/JuliaPlots/StatsPlots.jl).

## Simple example 

The following simple example illustrates how to use Chain to visually summarize a MCMC simulation:

```@example statsplots
using MCMCChains
using StatsPlots; gr()
Plots.reset_defaults() # hide
# Thanks to https://github.com/JuliaPlots/Plots.jl/issues/897. # hide
upscale = 8 # hide
default(size=(840,600)) # hide
# To avoid vectorized graphics for these crowded images. # hide
default(fmt = :png) # hide

theme(:ggplot2)

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
plot!(size=(840, 600), fmt = :png) # hide
```

\

```@example statsplots
p2 = plot(chn, colordim = :parameter)
plot!(size=(840, 400), fmt = :png) # hide
```

\
Note that the plot function takes the additional arguments described in the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package.

## Mixed density

```@example statsplots
plot(chn, seriestype = :mixeddensity)
```

Or, for all seriestypes, use the alternative shorthand syntax:

```@example statsplots
mixeddensity(chn)
```

## Trace

```@example statsplots
plot(chn, seriestype = :traceplot)
```

```@example statsplots
traceplot(chn)
```

## Running average

```@example statsplots
meanplot(chn)
```

## Density

```@example statsplots
density(chn)
```

## Histogram

```@example statsplots
histogram(chn)
```

## Autocorrelation

```@example statsplots
autocorplot(chn)
```

## Corner

**TODO:** Need to set the right params for this one.

```@example statsplots
# corner(chn, ["param_1", "param_2"])
```
