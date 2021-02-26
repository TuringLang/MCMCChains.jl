# StatsPlots.jl

MCMCChains implements many functions for plotting via [StatsPlots.jl](https://github.com/JuliaPlots/StatsPlots.jl).

## Simple example 

The following simple example illustrates how to use Chain to visually summarize a MCMC simulation:

```@example statsplots
using MCMCChains
using StatsPlots

theme(:ggplot2)

# Define the experiment
n_iter = 500
n_name = 3
n_chain = 2

# experiment results
val = randn(n_iter, n_name, n_chain) .+ [1, 2, 3]'
val = hcat(val, rand(1:2, n_iter, 1, n_chain))

# construct a Chains object
chn = Chains(val, [:A, :B, :C, :D])

# visualize the MCMC simulation results
plot(chn; size=(840, 600))
# This output is used in README.md too. # hide
filename = "default_plot.svg" # hide
savefig(filename); nothing # hide
```

![Default plot for Chains](default_plot.svg)
\

```@example statsplots
p2 = plot(chn, colordim = :parameter)
# Using png to avoid vectorized graphics for images with many vectors. # hide
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
plot!(fmt = :png) # hide
```

```@example statsplots
traceplot(chn)
plot!(fmt = :png) # hide
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
