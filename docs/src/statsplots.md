# StatsPlots.jl

MCMCChains implements many functions for plotting via [StatsPlots.jl](https://github.com/JuliaPlots/StatsPlots.jl).

## Simple example

The following simple example illustrates how to use Chain to visually summarize a MCMC simulation:

```@example statsplots
using MCMCChains
using StatsPlots

# Define the experiment
n_iter = 100
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

```@example statsplots
plot(chn, colordim = :parameter; size=(840, 400))
```

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

```@example statsplots
corner(chn)
```

For plotting multiple parameters, ridgeline, forest and caterpillar plots can be useful.
Please see the docstrings for a detailed description of these functions.

## Ridgeline

```@example statsplots
ridgelineplot(chn, [:C, :B, :A])
```

## Forest

```@example statsplots
forestplot(chn, [:C, :B, :A], hpd_val = [0.05, 0.15, 0.25])
```

## Caterpillar

```@example statsplots
forestplot(chn, chn.name_map[:parameters], hpd_val = [0.05, 0.15, 0.25], ordered = true)
```
