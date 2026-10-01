@shorthands meanplot
@shorthands autocorplot
@shorthands rankplot
@shorthands mixeddensity
@shorthands pooleddensity
@shorthands traceplot
@shorthands corner
@shorthands violinplot

"""
    energyplot(chains::Chains; kind=:density, kwargs...)

Generate an energy plot for the samples in `chains`.

The energy plot is a diagnostic tool for HMC-based samplers like NUTS. It displays the distributions of the Hamiltonian energy and the energy transition (error) to diagnose sampler efficiency and identify divergences.

This plot is only available for chains that contain the `:hamiltonian_energy` and `:hamiltonian_energy_error` statistics in their `:internals` section.

# Keywords
- `kind::Symbol` (default: `:density`): The type of plot to generate. Can be `:density` or `:histogram`.
"""
@userplot EnergyPlot

"""
    ppcplot(posterior_chains::Chains, posterior_predictive_chains::Chains, observed_data::Vector; kwargs...)

Generate a posterior/prior predictive check (PPC) plot comparing observed data with predictive samples.

PPC plots are a key tool for model validation in Bayesian analysis. They help assess whether the model 
can reproduce the key features of the observed data by comparing the observed data against samples from 
the posterior (or prior) predictive distribution.

# Arguments
- `posterior_chains::Chains`: MCMC samples from the posterior (or prior) distribution
- `posterior_predictive_chains::Chains`: Samples from the posterior (or prior) predictive distribution
- `observed_data::Vector`: The observed data values

# Keywords
- `kind::Symbol` (default: `:density`): Type of plot - `:density`, `:histogram`, `:scatter`, or `:cumulative`
- `alpha::Real` (default: `0.2` for density/cumulative, `0.7` for scatter): Transparency of predictive curves
- `num_pp_samples::Integer` (default: all samples): Number of predictive samples to plot
- `mean_pp::Bool` (default: `true`): Whether to plot the mean of predictive distribution
- `observed::Bool` (default: `true` for posterior, `false` for prior): Whether to plot observed data
- `observed_rug::Bool` (default: `false`): Whether to add a rug plot for observed data (kde/cumulative only)
- `colors::Vector` (default: `[:steelblue, :black, :orange]`): Colors for [predictive, observed, mean_pp]
- `jitter::Real` (default: `0.0`, `0.7` for scatter with ≤5 samples): Jitter amount for scatter plots
- `legend::Bool` (default: `true`): Whether to show legend
- `random_seed::Integer` (default: `nothing`): Random seed for reproducible subsampling
- `ppc_group::Symbol` (default: `:posterior`): Specify `:posterior` or `:prior` for appropriate defaults and labeling

# Examples
```julia
# Posterior Predictive Check
ppcplot(posterior_chains, posterior_predictive_chains, observed_data)

# Prior Predictive Check (observed data not shown by default)
ppcplot(prior_chains, prior_predictive_chains, observed_data; ppc_group=:prior)

# Histogram
ppcplot(chains, pp_chains, observed_data; kind=:histogram)

# Cumulative distribution
ppcplot(chains, pp_chains, observed_data; kind=:cumulative)

# Scatter plot with jitter
ppcplot(chains, pp_chains, observed_data; kind=:scatter, jitter=0.5)

# Prior check with observed data shown
ppcplot(prior_chains, pp_chains, observed_data; ppc_group=:prior, observed=true)

# Subset of predictive samples with custom colors
ppcplot(chains, pp_chains, observed_data; 
        num_pp_samples=20, 
        colors=[:blue, :red, :green], 
        random_seed=42)
```

# Notes
The `ppc_group` parameter controls default behavior:
- `:posterior`: Shows observed data by default, uses "Posterior Predictive Check" title
- `:prior`: Hides observed data by default, uses "Prior Predictive Check" title
"""
@userplot PPCPlot

"""
    ridgelineplot(chains::Chains[, params::Vector{Symbol}]; kwargs...)

Generate a ridgeline plot for the samples of the parameters `params` in `chains`.

By default, all parameters are plotted.

## Keyword arguments

The following options are available:

- `fill_q` (default: `false`) and `fill_hpd` (default: `true`):
  Fill the area below the curve in the quantiles interval (`fill_q = true`) or the highest posterior density (HPD) interval (`fill_hpd = true`).
  If both `fill_q = false` and `fill_hpd = false`, then the whole area below the curve is filled.
  If no fill color is desired, it should be specified with series attributes.
  These options are mutually exclusive.

- `show_mean` (default: `true`) and `show_median` (default: `true`):
  Plot a vertical line of the mean (`show_mean = true`) or median (`show_median = true`) of the posterior density estimate.
  If both options are set to `true`, both lines are plotted.

- `show_qi` (default: `false`) and `show_hpdi` (default: `true`):
  Plot a quantile interval (`show_qi = true`) or the largest HPD interval (`show_hpdi = true`) at the bottom of each density plot.
  These options are mutually exclusive.

- `q` (default: `[0.1, 0.9]`): The two quantiles used for plotting if `fill_q = true` or `show_qi = true`.

- `hpd_val` (default: `[0.05, 0.2]`): The complementary probability mass(es) of the highest posterior density intervals that are plotted if `fill_hpd = true` or `show_hpdi = true`.

!!! note
    If a single parameter is provided, the generated plot is a density plot with all the elements described above.
"""
@userplot RidgelinePlot

"""
    forestplot(chains::Chains[, params::Vector{Symbol}]; kwargs...)

Generate a forest or caterpillar plot for the samples of the parameters `params` in `chains`.

By default, all parameters are plotted.

## Keyword arguments

- `ordered` (default: `false`):
  If `ordered = false`, a forest plot is generated.
  If `ordered = true`, a caterpillar plot is generated.

- `fill_q` (default: `false`) and `fill_hpd` (default: `true`):
  Fill the area below the curve in the quantiles interval (`fill_q = true`) or the highest posterior density (HPD) interval (`fill_hpd = true`).
  If both `fill_q = false` and `fill_hpd = false`, then the whole area below the curve is filled.
  If no fill color is desired, it should be specified with series attributes.
  These options are mutually exclusive.

- `show_mean` (default: `true`) and `show_median` (default: `true`):
  Plot a vertical line of the mean (`show_mean = true`) or median (`show_median = true`) of the posterior density estimate.
  If both options are set to `true`, both lines are plotted.

- `show_qi` (default: `false`) and `show_hpdi` (default: `true`):
  Plot a quantile interval (`show_qi = true`) or the largest HPD interval (`show_hpdi = true`) at the bottom of each density plot.
  These options are mutually exclusive.

- `q` (default: `[0.1, 0.9]`): The two quantiles used for plotting if `fill_q = true` or `show_qi = true`.

- `hpd_val` (default: `[0.05, 0.2]`): The complementary probability mass(es) of the highest posterior density intervals that are plotted if `fill_hpd = true` or `show_hpdi = true`.
"""
@userplot ForestPlot

"""
    _interval_args(p, name)

Chain and parameter names for a ridgeline or forest plot.

`@userplot` leaves `p.args` untyped, so without this a missing parameter list surfaces as a
`BoundsError` from inside the recipe. Omitting the list plots every parameter.
"""
function _interval_args(p, name::AbstractString)
    if length(p.args) == 1
        chn = only(p.args)
        chn isa Chains ||
            throw(ArgumentError("$name expects a Chains as its first argument"))
        return chn, names(chn, :parameters)
    elseif length(p.args) == 2
        chn, par_names = p.args
        chn isa Chains ||
            throw(ArgumentError("$name expects a Chains as its first argument"))
        return chn, par_names
    end
    throw(
        ArgumentError(
            "$name expects a Chains and optionally a vector of parameter names, got $(length(p.args)) arguments",
        ),
    )
end
