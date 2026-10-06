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

Every chain of a parameter is binned over the same edges, so the bars line up and the counts can be compared.

The bin count comes from the draws, taking whichever of [Freedman and Diaconis (1981)](https://doi.org/10.1007/BF01025868) and Sturges asks for more.
Pass `bins` to override it, as a rule, a count, or the edges themselves.

```@example statsplots
histogram(chn, bins = :sturges)
```

```@example statsplots
histogram(chn, bins = 60, fillalpha = 0.7)
```

## Autocorrelation

```@example statsplots
autocorplot(chn)
```

## Rank

A rank plot histograms the ranks of the draws pooled over all chains, one line per chain.
If every chain is sampling the same posterior then the ranks are uniform, so each chain should stay near the dashed line at the expected count.
Rank plots keep working on long chains, where a trace plot turns into a band of ink ([Vehtari et al. 2021](https://doi.org/10.1214/20-BA1221)).

```@example statsplots
rankplot(chn)
```

Use `nbins` to set how many bins the ranks are collected into.

```@example statsplots
rankplot(chn, nbins = 10)
```

## Empirical CDF

```@example statsplots
plot(chn, seriestype = :ecdfplot)
```

Two chains that agree lie on top of each other here, which is easier to judge than two density curves, because the eye compares positions rather than areas.
The name `ecdfplot` belongs to StatsPlots and takes a vector, so a chain goes through `plot` with the series type.

## Convergence diagnostics

Most of the time the question is only whether anything went wrong, and `diagnosticsplot` answers that in one grid.
Each cell is coloured by whether the value is fine, worth a look, or bad, and grey where the diagnostic does not apply.

```@example statsplots
diagnosticsplot(chn)
```

`essplot`, `rhatplot` and `mcseplot` draw one dot per parameter for the diagnostics that [`ess`](@ref), [`rhat`](@ref) and [`mcse`](@ref) already report.
The dashed line is the threshold recommended by [Vehtari et al. (2021)](https://doi.org/10.1214/20-BA1221), an effective sample size of 400 and an R-hat of 1.01.

```@example statsplots
essplot(chn)
```

`kind` picks which effective sample size to estimate.

```@example statsplots
essplot(chn, kind = :tail)
```

`relative = true` divides by the number of draws, which keeps the axis readable when chains are long.

```@example statsplots
essplot(chn, relative = true)
```

```@example statsplots
rhatplot(chn)
```

```@example statsplots
mcseplot(chn, relative = true)
```

A single number says how the run ended.
`evolutionplot` says whether it was still improving, which is what tells you that running for longer would help.

```@example statsplots
evolutionplot(chn)
```

```@example statsplots
evolutionplot(chn, diagnostic = :rhat)
```

## Parallel coordinates

`parallelplot` draws one line per draw across all parameters, standardised so that parameters on different scales share one axis.

```@example statsplots
parallelplot(chn)
```

When the chain records divergent transitions they are drawn on top in a second colour, which is what makes this plot worth reading: it shows where in the parameter space the sampler is failing.

```@example statsplots
parallelplot(chn, num_draws = 100, standardise = false)
```

## Violin

Violin plots are similar to box plots but also show the probability density of the data at different values, smoothed by a kernel density estimator.

```@example statsplots
violinplot(chn) # Plotting parameter 1 across all chains
```

```@example statsplots
violinplot(chn, 1) # Plotting parameter 1 across all chains
```

```@example statsplots
violinplot(chn, :A) # Plotting a specific parameter across all chains
```

```@example statsplots
violinplot(chn, [:C, :B, :A]) # Plotting multiple specific parameters across all chains
```

```@example statsplots
violinplot(chn, 1, colordim = :parameter) # Plotting chain 1 across all parameters
```

```@example statsplots
violinplot(chn, show_boxplot = false) # Plotting all parameters without the inner boxplot
```

You can also aggregate (pool) samples from all chains for a given parameter by using `append_chains = true`. This is useful when you want to visualize the overall posterior distribution without distinguishing between individual chains.

```@example statsplots
violinplot(chn, :A, append_chains = true) # Single parameter, all chains appended
```

```@example statsplots
violinplot(chn, append_chains = true) # All parameters, all chains appended
```

You can also use the `plot` function with `seriestype = :violinplot` or `seriestype = :violin`

```@example statsplots
plot(chn, seriestype = :violin)
```

## Corner

```@example statsplots
corner(chn)
```

## Energy Plot

The energy plot is a diagnostic tool for HMC-based samplers (like NUTS) that helps diagnose sampling efficiency by visualizing the energy and energy transition distributions. This plot requires that the chain contains the internal sampler statistics `:hamiltonian_energy` and `:hamiltonian_energy_error`.

```@example statsplots
# First, we generate a chain that includes the required sampler parameters.
n_iter = 1000
n_chain = 4
val_params = randn(n_iter, 2, n_chain)
val_energy = randn(n_iter, 1, n_chain) .+ 20
val_energy_error = randn(n_iter, 1, n_chain) .* 0.5
full_val = hcat(val_params, val_energy, val_energy_error)

parameter_names = [:a, :b, :hamiltonian_energy, :hamiltonian_energy_error]
section_map = (
    parameters=[:a, :b],
    internals=[:hamiltonian_energy, :hamiltonian_energy_error],
)

chn_energy = Chains(full_val, parameter_names, section_map)

# Generate the energy plot (default is a density plot).
energyplot(chn_energy)
```

```@example statsplots
# The plot can also be generated as a histogram.
energyplot(chn_energy, kind=:histogram)
```

## Sampler diagnostics

`nutsplot` draws the statistics an HMC sampler records about its own behaviour, one series per chain.

```@example statsplots
# A chain carrying the sampler statistics NUTS records.
n_iter = 500
n_chain = 4
nuts_val = hcat(
    randn(n_iter, 1, n_chain),
    clamp.(0.9 .+ 0.08 .* randn(n_iter, 1, n_chain), 0, 1),
    0.35 .+ 0.02 .* randn(n_iter, 1, n_chain),
    float.(rand(2:5, n_iter, 1, n_chain)),
    float.(rand(n_iter, 1, n_chain) .< 0.03),
)
chn_nuts = Chains(
    nuts_val,
    [:a, :acceptance_rate, :step_size, :tree_depth, :numerical_error],
    (
        parameters = [:a],
        internals = [:acceptance_rate, :step_size, :tree_depth, :numerical_error],
    ),
)

nutsplot(chn_nuts)
```

```@example statsplots
nutsplot(chn_nuts, kind = :treedepth)
```

```@example statsplots
nutsplot(chn_nuts, kind = :divergence)
```

For plotting multiple parameters, ridgeline, forest and caterpillar plots can be useful.

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

## Posterior Predictive Checks (PPC)

Posterior Predictive Checks (PPC) are essential tools for Bayesian model validation. They compare observed data with samples from the posterior predictive distribution to assess whether the model can reproduce key features of the data. Prior Predictive Checks can also be performed to evaluate prior appropriateness before seeing the data.

```@example statsplots
using Random
Random.seed!(123)

# Generate posterior samples (parameters)
n_iter = 500
posterior_data = randn(n_iter, 2, 2)  # μ, σ parameters
posterior_chains = Chains(posterior_data, [:μ, :σ])

# Generate posterior predictive samples
n_obs = 20
pp_data = zeros(n_iter, n_obs, 2)
for i in 1:n_iter, j in 1:2
    μ = posterior_data[i, 1, j]
    σ = abs(posterior_data[i, 2, j]) + 0.5  # Ensure positive σ
    pp_data[i, :, j] = randn(n_obs) * σ .+ μ
end
pp_chains = Chains(pp_data)

# Generate observed data
Random.seed!(456)
observed = randn(n_obs) * 1.2 .+ 0.3

# Basic posterior predictive check (density overlay)
# Note: observed data is shown by default for posterior checks
ppcplot(posterior_chains, pp_chains, observed)
```

### Plot Types

Our PPC implementation supports four main plot types:

#### Density Plots (Default)
```@example statsplots
# Density overlay with customized transparency
ppcplot(posterior_chains, pp_chains, observed; 
        kind=:density, alpha=0.3, num_pp_samples=50)
```

#### Histogram Comparison
```@example statsplots
# Normalized histogram comparison
ppcplot(posterior_chains, pp_chains, observed; kind=:histogram)
```

#### Cumulative Distribution Functions
```@example statsplots
# Empirical CDFs comparison
ppcplot(posterior_chains, pp_chains, observed; kind=:cumulative)
```

#### Scatter Plots with Jitter
```@example statsplots
# Index-based scatter plot with automatic jitter for small samples
ppcplot(posterior_chains, pp_chains, observed; 
        kind=:scatter, num_pp_samples=8, jitter=0.3)
```

### Advanced Styling and Options

```@example statsplots
# Comprehensive customization example
ppcplot(posterior_chains, pp_chains, observed; 
        kind=:density,
        colors=[:steelblue, :darkred, :orange],  # [predictive, observed, mean]
        alpha=0.25,
        observed_rug=true,      # Add rug plot for observed data
        num_pp_samples=75,      # Limit predictive samples shown
        mean_pp=true,           # Show predictive mean
        legend=true,
        random_seed=42)         # Reproducible subsampling
```

### Prior Predictive Checks

Prior predictive checks assess whether priors generate reasonable data before observing actual data. The `ppc_group` parameter controls default behavior:

```@example statsplots
# Prior predictive check - observed data hidden by default
ppcplot(posterior_chains, pp_chains, observed; ppc_group=:prior)
```

```@example statsplots
# Prior check with observed data explicitly shown for comparison
ppcplot(posterior_chains, pp_chains, observed; 
        ppc_group=:prior, observed=true, alpha=0.4)
```

### Controlling Observed Data Display

You can explicitly control whether observed data is shown regardless of the check type:

```@example statsplots
# Posterior check without observed data
ppcplot(posterior_chains, pp_chains, observed; 
        ppc_group=:posterior, observed=false)
```

### Performance and Sampling Control

For large datasets or when you want to reduce visual clutter:

```@example statsplots
# Limit the number of predictive samples displayed
ppcplot(posterior_chains, pp_chains, observed; 
        num_pp_samples=25, 
        random_seed=123)  # Reproducible results
```

```julia
ppcplot(posterior_chains::Chains, posterior_predictive_chains::Chains, observed_data::Vector;
        kind=:density, alpha=nothing, num_pp_samples=nothing, mean_pp=true, observed=nothing,
        observed_rug=false, colors=[:steelblue, :black, :orange], jitter=nothing, 
        legend=true, random_seed=nothing, ppc_group=:posterior)
```

## Appearance

Every plot on this page is drawn with the same look: no grid, no box around the plot or the legend, grey axes and text, and chains coloured from the Okabe-Ito palette, which stays readable for the common forms of colour blindness.

`plot_style!` changes that for every plot drawn after it.

```@example statsplots
plot_style!(grid = :dots, framestyle = :box)
meanplot(chn)
```

```@example statsplots
plot_style!(grid = :none, framestyle = :axes, palette = ["#332288", "#117733", "#DDCC77", "#CC6677"])
meanplot(chn)
```

`background = :transparent` draws no background at all, which is what you want when the page behind the plot is not white.

Everything `plot_style!` sets is a default, so anything passed to a single plot still wins.

```@example statsplots
reset_plot_style!()
meanplot(chn, grid = true, legend = :outertopright)
```

By default only the first panel of a multi-panel plot names the chains, since the same legend on every panel is ink that says nothing new and lands on top of the data.
Set `legend_panel` to `:all` or `:none` to change that.

```@example statsplots
reset_plot_style!()
nothing # hide
```

## API

```@docs
energyplot
energyplot!
ppcplot
ppcplot!
ridgelineplot
ridgelineplot!
forestplot
forestplot!
diagnosticsplot
diagnosticsplot!
essplot
essplot!
rhatplot
rhatplot!
mcseplot
mcseplot!
evolutionplot
evolutionplot!
parallelplot
parallelplot!
nutsplot
nutsplot!
MCMCChains.CHAIN_PALETTE
plot_style
plot_style!
reset_plot_style!
```
