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
\

```@example statsplots
plot(chn, colordim = :parameter; size=(840, 400))
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

```@example statsplots
corner(chn)
```

For plotting multiple parameters, Ridgeline, Forest and Caterpillar plots can be useful.
Please see the docstrings for a detailed description of these functions.

## Ridgeline

```@example statsplots
ridgelineplot(chn, [:C, :B, :A])
```
""""
    ridgelineplot(chains::Chains, par_names::Vector{Symbol}; hpd_val = [0.05, 0.2],
    q = [0.1, 0.9], spacer = 0.5, _riser = 0.2, show_mean = true, show_median = true,
    show_qi = false, show_hpdi = true, fill_q = true, fill_hpd = false, ordered = false)

Given a `chains` object, returns a Ridgeline plot for the sampled parameters specified on
`par_names`.

For ridgeline plots, the following parameters are defined:

** (a) Fill **
Fill area below the curve can be determined by quantiles interval (`fill_q = true`) or
hdpi interval (`fill_hpd = true`). Default options are `fill_hpd = true` and `fill_q = false`.
If both `fill_q = false` and `fill_hpd = false`, then the whole area below the curve will be
filled. If no fill color is desired, it should be specified with series attributes. These
fill options are mutually exclusive.

** (b) Mean and median **
A vertical line can be plotted repesenting the mean (`show_mean = true`) or median
(`show_median = true`) of the density (kde) distribution. Both options can be plotted at the
 same time.

** (c) Intervals **
At the bottom of each density plot, a quantile interval (`show_qi = true`) or HPD interval
(`show_hdpi = true`) can be plotted. These options are mutually exclusive. Default options
are `show_qi = false` and `show_hpdi = true`.To plot quantile intervals, the values specified
as `q` will be taken, and for HPD intervals, only the smaller value specified in `hpd_val`
will be used.

Note: When one parameter is given, it will be plotted as a density plot with all the elements
described above.
"""

## Forest

```@example statsplots
forestplot(chn, [:C, :B, :A], hpd_val = [0.05, 0.15, 0.25])
```

"""
    forestplot(chains::Chains, par_names::Vector{Symbol}; hpd_val = [0.05, 0.2],
    q = [0.1, 0.9], spacer = 0.5, _riser = 0.2, show_mean = true, show_median = true,
    show_qi = false, show_hpdi = true, fill_q = true, fill_hpd = false, ordered = false)

Given a `chains` object, returns a Forest plot for the sampled parameters specified on
`par_names`.

Both `forest` and `caterpillar` plots are called using `forestplot` shorthands.
If `ordered = false`, then a `forest` plot will be returned, and if `ordered = true`,
a `caterpillar` plot will be returned.

For both plot types the following elemets can be plotted:

**High posterior density intervals (HPDI)** determined by the number of elements in `hpd_val`.
All the values in `hpd_val` will be used to construct the intervals with `MCMCChains.hpd.

**Quantile intervals** determined by the 2-element vector `q`.

**Mean and/or median.** Plotted as points with different `markershape..

"""

## Caterpillar

```@example statsplots
forestplot(chn, chn.name_map[:parameters], hpd_val = [0.05, 0.15, 0.25], ordered = true)
```
