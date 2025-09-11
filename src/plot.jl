@shorthands meanplot
@shorthands autocorplot
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

struct _TracePlot
    c::Any
    val::Any
end
struct _MeanPlot
    c::Any
    val::Any
end
struct _DensityPlot
    c::Any
    val::Any
end
struct _HistogramPlot
    c::Any
    val::Any
end
struct _AutocorPlot
    lags::Any
    val::Any
end
struct _ViolinPlot
    c::Any
    val::Any
    # param_indices: For accurate x-axis labeling (parameter names vs. chain indices).
    param_indices::Any
    # show_boxplot: To allow toggling the inner boxplot visibility.
    show_boxplot::Any
    # colordim: To guide x-axis label generation based on grouping (by chain or parameter).
    colordim::Any
end

# define alias functions for old syntax
const translationdict = Dict(
    :traceplot => _TracePlot,
    :meanplot => _MeanPlot,
    :density => _DensityPlot,
    :histogram => _HistogramPlot,
    :autocorplot => _AutocorPlot,
    :pooleddensity => _DensityPlot,
    :violinplot => _ViolinPlot,
    :violin => _ViolinPlot,
)

const supportedplots = push!(collect(keys(translationdict)), :mixeddensity, :corner)

@recipe f(c::Chains, s::Symbol) = c, [s]

@recipe function f(
    chains::Chains,
    i::Int;
    colordim = :chain,
    barbounds = (-Inf, Inf),
    maxlag = nothing,
    append_chains = false,
)
    st = get(plotattributes, :seriestype, :traceplot)
    c = append_chains || st == :pooleddensity ? pool_chain(chains) : chains

    if colordim == :parameter
        title --> "Chain $(MCMCChains.chains(c)[i])"
        label --> permutedims(map(string, names(c)))
        val = c.value[:, :, i]
        _actual_indices_for_violin = 1:size(c, 2)
    elseif colordim == :chain
        title --> string(names(c)[i])
        label --> permutedims(map(x -> "Chain $x", MCMCChains.chains(c)))
        val = c.value[:, i, :]
        _actual_indices_for_violin = MCMCChains.chains(c)
    else
        throw(ArgumentError("`colordim` must be one of `:chain` or `:parameter`"))
    end

    if st == :mixeddensity || st == :pooleddensity
        discrete = indiscretesupport(c, barbounds)
        st = if colordim == :chain
            discrete[i] ? :histogram : :density
        else
            # NOTE: It might make sense to overlay histograms and density plots here.
            :density
        end
        seriestype := st
    end

    if st == :autocorplot
        lags = 0:(maxlag === nothing ? round(Int, 10 * log10(length(range(c)))) : maxlag)
        # Chains are already appended in `c` if desired, hence we use `append_chains=false`
        ac = autocor(c; sections = nothing, lags = lags, append_chains = false)
        ac_mat = convert(Array{Float64}, ac)
        val = colordim == :parameter ? ac_mat[:, :, i]' : ac_mat[i, :, :]
        _AutocorPlot(lags, val)
    elseif st ∈ (:violinplot, :violin)
        show_boxplot_kw = get(plotattributes, :show_boxplot, true)
        # Passes `_actual_indices_for_violin`` to `_ViolinPlot`` to ensure correct x-axis labels (parameter names or chain numbers).
        return _ViolinPlot(c, val, _actual_indices_for_violin, show_boxplot_kw, colordim)
    elseif st ∈ keys(translationdict)
        translationdict[st](c, val)
    else
        range(c), val
    end
end

@recipe function f(p::_DensityPlot)
    xaxis --> "Sample value"
    yaxis --> "Density"
    trim --> true
    [collect(skipmissing(p.val[:, k])) for k = 1:size(p.val, 2)]
end

@recipe function f(p::_HistogramPlot)
    xaxis --> "Sample value"
    yaxis --> "Frequency"
    fillalpha --> 0.7
    bins --> 25
    trim --> true
    [collect(skipmissing(p.val[:, k])) for k = 1:size(p.val, 2)]
end

@recipe function f(p::_MeanPlot)
    seriestype := :path
    xaxis --> "Iteration"
    yaxis --> "Mean"
    range(p.c), cummean(p.val)
end

@recipe function f(p::_AutocorPlot)
    seriestype := :path
    xaxis --> "Lag"
    yaxis --> "Autocorrelation"
    p.lags, p.val
end

@recipe function f(p::_TracePlot)
    seriestype := :path
    xaxis --> "Iteration"
    yaxis --> "Sample value"
    range(p.c), p.val
end

@recipe function f(p::_ViolinPlot)
    num_series = size(p.val, 2)
    flat_data = vcat([collect(skipmissing(p.val[:, k])) for k = 1:num_series]...)

    plot_labels = String[]

    if p.colordim == :parameter
        plot_labels = string.(MCMCChains.names(p.c)[p.param_indices])
    elseif p.colordim == :chain
        plot_labels = ["Chain $(c_idx)" for c_idx in p.param_indices]
    else
        plot_labels = string.(1:num_series)
    end

    group_labels = repeat(1:num_series, inner = size(p.val, 1))

    xticks := (1:num_series, plot_labels)
    legend --> false

    @series begin
        seriestype := :violin
        x := group_labels
        y := flat_data
        group := group_labels
        ()
    end

    if p.show_boxplot
        @series begin
            seriestype := :boxplot
            bar_width := 0.1
            linewidth := 2
            fillalpha := 0.8
            x := group_labels
            y := flat_data
            group := group_labels
            ()
        end
    end
end

@recipe function f(p::EnergyPlot; kind = :density)
    chains = p.args[1]

    if kind ∉ (:density, :histogram)
        error("`kind` must be one of `:density` or `:histogram`")
    end

    internal_names = names(chains, :internals)
    required_params = [:hamiltonian_energy, :hamiltonian_energy_error]
    for param in required_params
        if param ∉ internal_names
            error(
                "`$param` not found in chain's internal parameters. Energy plots are only available for HMC/NUTS samplers.",
            )
        end
    end

    pooled = pool_chain(chains)
    energy = vec(pooled[:, :hamiltonian_energy, :])
    energy_error = vec(pooled[:, :hamiltonian_energy_error, :])

    mean_energy = mean(energy)
    std_energy = std(energy)
    centered_energy = (energy .- mean_energy) ./ std_energy
    scaled_energy_error = energy_error ./ std_energy

    title := "Energy Plot"
    xaxis := "Standardized Energy"
    yaxis := "Density"
    legend := :topright

    @series begin
        seriestype := kind
        label := "Marginal Energy"
        fillrange --> 0
        fillalpha --> 0.5
        normalize --> true
        bins --> 50
        centered_energy
    end

    @series begin
        seriestype := kind
        label := "Energy Transition"
        fillrange --> 0
        fillalpha --> 0.5
        normalize --> true
        bins --> 50
        scaled_energy_error
    end
end

@recipe function f(
    p::PPCPlot;
    kind = :density,
    alpha = nothing,
    num_pp_samples = nothing,
    mean_pp = true,
    observed = nothing,
    observed_rug = false,
    colors = [:steelblue, :black, :orange],
    jitter = nothing,
    legend = true,
    random_seed = nothing,
    ppc_group = :posterior,
)
    if length(p.args) < 3
        error("ppcplot requires at least 3 arguments: (posterior_chains, posterior_predictive_chains, observed_data)")
    end
    
    posterior_chains = p.args[1]
    pp_chains = p.args[2]
    observed_data = p.args[3]
    
    if !(posterior_chains isa Chains)
        error("First argument must be a Chains object (posterior chains)")
    end
    if !(pp_chains isa Chains) 
        error("Second argument must be a Chains object (posterior predictive chains)")
    end
    if !(observed_data isa AbstractVector)
        error("Third argument must be a vector (observed data)")
    end
    
    if kind ∉ (:density, :histogram, :scatter, :cumulative)
        error("`kind` must be one of `:density`, `:histogram`, `:scatter`, or `:cumulative`")
    end
    
    if ppc_group ∉ (:posterior, :prior)
        error("`ppc_group` must be one of `:posterior` or `:prior`")
    end
    
    if observed === nothing
        observed = (ppc_group == :posterior)
    end
    
    if length(colors) != 3
        error("`colors` must be a vector of length 3: [predictive_color, observed_color, mean_color]")
    end
    
    if alpha === nothing
        alpha = (kind == :scatter) ? 0.7 : 0.2
    end
    
    if jitter === nothing
        jitter = 0.0
    end
    
    pp_pooled = pool_chain(pp_chains)
    pp_data = Array(pp_pooled.value.data)
    pp_data = pp_data[:, :, 1]
    
    if random_seed !== nothing
        Random.seed!(random_seed)
    end
    
    total_pp_samples = size(pp_data, 1)
    if num_pp_samples === nothing
        if kind == :scatter
            num_pp_samples = min(5, total_pp_samples)
            if jitter == 0.0 && num_pp_samples <= 5
                jitter = 0.7
            end
        else
            num_pp_samples = total_pp_samples
        end
    end
    
    if num_pp_samples > total_pp_samples
        @warn "Requested $num_pp_samples samples but only $total_pp_samples available. Using all samples."
        num_pp_samples = total_pp_samples
    end
    
    if num_pp_samples < total_pp_samples
        sample_indices = Random.randperm(total_pp_samples)[1:num_pp_samples]
        pp_data = pp_data[sample_indices, :]
    end
    
    if ppc_group == :prior
        title := "Prior Predictive Check"
        predictive_label = "Prior Predictive"
        mean_label = "Prior Predictive Mean"
    else
        title := "Posterior Predictive Check" 
        predictive_label = "Posterior Predictive"
        mean_label = "Posterior Predictive Mean"
    end
    legend := legend
    
    if kind == :density
        xaxis := "Value"
        yaxis := "Density"
        
        for i in 1:size(pp_data, 1)
            @series begin
                seriestype := :density
                label := i == 1 ? predictive_label : ""
                color := colors[1]
                alpha := alpha
                linewidth := 1
                pp_data[i, :]
            end
        end
        
        if observed
            @series begin
                seriestype := :density
                label := "Observed"
                color := colors[2]
                linewidth := 2
                observed_data
            end
        end
        
        if mean_pp
            pp_mean = vec(mean(pp_data, dims=1))
            @series begin
                seriestype := :density
                label := mean_label
                color := colors[3]
                linewidth := 2
                linestyle := :dash
                pp_mean
            end
        end
        
        if observed_rug && observed
            y_min = 0
            @series begin
                seriestype := :scatter
                label := ""
                color := colors[2]
                markershape := :vline
                markersize := 2
                y := fill(y_min, length(observed_data))
                observed_data
            end
        end
        
    elseif kind == :histogram
        xaxis := "Value"
        yaxis := "Frequency"
        
        for i in 1:size(pp_data, 1)
            @series begin
                seriestype := :histogram
                label := i == 1 ? predictive_label : ""
                color := colors[1]
                alpha := alpha
                normalize := true
                bins := 30
                pp_data[i, :]
            end
        end
        
        if observed
            @series begin
                seriestype := :histogram
                label := "Observed"
                color := colors[2]
                alpha := 0.7
                normalize := true
                bins := 30
                observed_data
            end
        end
        
        if mean_pp
            pp_mean = vec(mean(pp_data, dims=1))
            @series begin
                seriestype := :histogram
                label := mean_label
                color := colors[3]
                alpha := 0.7
                normalize := true
                bins := 30
                pp_mean
            end
        end
        
    elseif kind == :cumulative
        xaxis := "Value"
        yaxis := "Cumulative Probability"
        
        all_data = vcat(vec(pp_data), observed_data)
        x_range = range(minimum(all_data), maximum(all_data), length=200)
        
        for i in 1:size(pp_data, 1)
            pp_ecdf = ecdf(pp_data[i, :])
            y_vals = pp_ecdf.(x_range)
            @series begin
                seriestype := :path
                label := i == 1 ? predictive_label : ""
                color := colors[1]
                alpha := alpha
                linewidth := 1
                x := x_range
                y := y_vals
                ()
            end
        end
        
        if observed
            obs_ecdf = ecdf(observed_data)
            obs_y_vals = obs_ecdf.(x_range)
            @series begin
                seriestype := :path
                label := "Observed"
                color := colors[2]
                linewidth := 2
                x := x_range
                y := obs_y_vals
                ()
            end
        end
        
        if mean_pp
            pp_mean = vec(mean(pp_data, dims=1))
            pp_mean_ecdf = ecdf(pp_mean)
            mean_y_vals = pp_mean_ecdf.(x_range)
            @series begin
                seriestype := :path
                label := mean_label
                color := colors[3]
                linewidth := 2
                linestyle := :dash
                x := x_range
                y := mean_y_vals
                ()
            end
        end
        
        if observed_rug && observed
            @series begin
                seriestype := :scatter
                label := ""
                color := colors[2]
                markershape := :vline
                markersize := 2
                x := observed_data
                y := fill(0, length(observed_data))
                ()
            end
        end
        
    elseif kind == :scatter
        xaxis := "Index"
        yaxis := "Value"
        
        if jitter > 0
            jitter_vals = jitter * (rand(size(pp_data, 2)) .- 0.5)
        else
            jitter_vals = zeros(size(pp_data, 2))
        end
        
        for i in 1:size(pp_data, 1)
            y_vals = pp_data[i, :] .+ (jitter > 0 ? jitter * (rand(length(pp_data[i, :])) .- 0.5) : 0)
            @series begin
                seriestype := :scatter
                label := i == 1 ? predictive_label : ""
                color := colors[1]
                alpha := alpha
                markersize := 3
                x := 1:length(y_vals)
                y := y_vals
                ()
            end
        end
        
        if observed
            obs_y = observed_data .+ jitter_vals[1:length(observed_data)]
            @series begin
                seriestype := :scatter
                label := "Observed"
                color := colors[2]
                markersize := 4
                markershape := :diamond
                x := 1:length(obs_y)
                y := obs_y
                ()
            end
        end
        
        if mean_pp
            pp_mean = vec(mean(pp_data, dims=1))
            mean_y = pp_mean .+ jitter_vals[1:length(pp_mean)]
            @series begin
                seriestype := :scatter
                label := mean_label
                color := colors[3]
                markersize := 4
                markershape := :square
                x := 1:length(mean_y)
                y := mean_y
                ()
            end
        end
    end
end

@recipe function f(chains::Chains, parameters::AbstractVector{Symbol}; colordim = :chain)
    colordim != :chain && error(
        "Symbol names are interpreted as parameter names, only compatible with ",
        "`colordim = :chain`",
    )

    ret = indexin(parameters, names(chains))
    any(y === nothing for y in ret) && error("Parameter not found")

    return chains, Int.(ret)
end

@recipe function f(
    chains::Chains,
    parameters::AbstractVector{<:Integer} = Int[];
    sections = _default_sections(chains),
    width = 500,
    height = 250,
    colordim = :chain,
    append_chains = false,
)
    _chains =
        isempty(parameters) ? Chains(chains, _clean_sections(chains, sections)) : chains
    c = append_chains ? pool_chain(_chains) : _chains
    ptypes = get(plotattributes, :seriestype, (:traceplot, :mixeddensity))
    ptypes = ptypes isa Symbol ? (ptypes,) : ptypes
    @assert all(ptype -> ptype ∈ supportedplots, ptypes)
    ntypes = length(ptypes)
    nrows, nvars, nchains = size(c)
    isempty(parameters) && (parameters = colordim == :chain ? (1:nvars) : (1:nchains))
    N = length(parameters)

    if :corner ∉ ptypes
        size --> (ntypes * width, N * height)
        legend --> false

        multiple_plots = N * ntypes > 1
        if multiple_plots
            layout := (N, ntypes)
        end

        i = 0
        for par in parameters
            for ptype in ptypes
                i += 1

                @series begin
                    if multiple_plots
                        subplot := i
                    end
                    colordim := colordim
                    seriestype := ptype
                    c, par
                end
            end
        end
    else
        ntypes > 1 && error(":corner is not compatible with multiple seriestypes")
        Corner(c, names(c)[parameters])
    end
end

struct Corner
    c::Any
    parameters::Any
end

@recipe function f(corner::Corner)
    # Convert labels to string because `Symbol` is not supported generally supported.
    label --> permutedims(map(string, corner.parameters))
    compact --> true
    size --> (600, 600)
    # NOTE: Don't use the indices from `chains(chains)`.
    # See https://github.com/TuringLang/MCMCChains.jl/issues/413.
    ar = collect(
        Array(corner.c.value[:, corner.parameters, i]) for i = 1:length(chains(corner.c))
    )
    RecipesBase.recipetype(:cornerplot, vcat(ar...))
end

function _compute_plot_data(
    i::Integer,
    chains::Chains,
    par_names::AbstractVector{Symbol};
    hpd_val = [0.05, 0.2],
    q = [0.1, 0.9],
    spacer = 0.4,
    _riser = 0.2,
    barbounds = (-Inf, Inf),
    show_mean = true,
    show_median = true,
    show_qi = false,
    show_hpdi = true,
    fill_q = true,
    fill_hpd = false,
    ordered = false,
)

    chain_dic = Dict(zip(quantile(chains)[:, 1], quantile(chains)[:, 4]))
    sorted_chain = sort(collect(zip(values(chain_dic), keys(chain_dic))))
    sorted_par = [sorted_chain[i][2] for i = 1:length(par_names)]
    par = (ordered ? sorted_par : par_names)
    hpdi = sort(hpd_val)

    chain_sections = MCMCChains.group(chains, Symbol(par[i]))
    chain_vec = vec(chain_sections.value.data)
    lower_hpd =
        [MCMCChains.hpd(chain_sections, alpha = hpdi[j]).nt.lower for j = 1:length(hpdi)]
    upper_hpd =
        [MCMCChains.hpd(chain_sections, alpha = hpdi[j]).nt.upper for j = 1:length(hpdi)]
    h = _riser + spacer * (i - 1)
    qs = quantile(chain_vec, q)
    k_density = kde(chain_vec)
    if fill_hpd
        x_int = filter(x -> lower_hpd[1][1] <= x <= upper_hpd[1][1], k_density.x)
        val = pdf(k_density, x_int) .+ h
    elseif fill_q
        x_int = filter(x -> qs[1] <= x <= qs[2], k_density.x)
        val = pdf(k_density, x_int) .+ h
    else
        x_int = k_density.x
        val = k_density.density .+ h
    end
    chain_med = median(chain_vec)
    chain_mean = mean(chain_vec)
    min = minimum(k_density.density .+ h)
    q_int = (show_qi ? [qs[1], chain_med, qs[2]] : [chain_med])

    return par,
    hpdi,
    lower_hpd,
    upper_hpd,
    h,
    qs,
    k_density,
    x_int,
    val,
    chain_med,
    chain_mean,
    min,
    q_int
end

@recipe function f(
    p::RidgelinePlot;
    hpd_val = [0.05, 0.2],
    q = [0.1, 0.9],
    spacer = 0.5,
    _riser = 0.2,
    show_mean = true,
    show_median = true,
    show_qi = false,
    show_hpdi = true,
    fill_q = true,
    fill_hpd = false,
    ordered = false,
)

    chn = p.args[1]
    par_names = p.args[2]

    for i = 1:length(par_names)
        par,
        hpdi,
        lower_hpd,
        upper_hpd,
        h,
        qs,
        k_density,
        x_int,
        val,
        chain_med,
        chain_mean,
        min,
        q_int = _compute_plot_data(
            i,
            chn,
            par_names;
            hpd_val = hpd_val,
            q = q,
            spacer = spacer,
            _riser = _riser,
            show_mean = show_mean,
            show_median = show_median,
            show_qi = show_qi,
            show_hpdi = show_hpdi,
            fill_q = fill_q,
            fill_hpd = fill_hpd,
            ordered = ordered,
        )

        yticks --> (
            length(par_names) > 1 ?
            (_riser .+ ((1:length(par_names)) .- 1) .* spacer, string.(par)) : :default
        )
        yaxis --> (length(par_names) > 1 ? "Parameters" : "Density")
        @series begin
            seriestype := :hline
            label := nothing
            linecolor := "#BBBBBB"
            linewidth --> 1.2
            [h]
        end
        @series begin
            seriestype := :path
            label := nothing
            fillrange --> min
            fillalpha --> 0.8
            x_int, val
        end
        @series begin
            seriestype := :path
            label := nothing
            linecolor --> "#000000"
            k_density.x, k_density.density .+ h
        end
        @series begin
            seriestype := :path
            label --> (show_mean ? (i == 1 ? "Mean" : nothing) : nothing)
            linecolor --> "dark red"
            linewidth --> (show_mean ? 1.2 : 0)
            [chain_mean, chain_mean], [min, min + pdf(k_density, chain_mean)]
        end
        @series begin
            seriestype := :path
            label --> (show_median ? (i == 1 ? "Median" : nothing) : nothing)
            linecolor --> "#000000"
            linewidth --> (show_median ? 1.2 : 0)
            [chain_med, chain_med], [min, min + pdf(k_density, chain_med)]
        end
        @series begin
            seriestype := :scatter
            label := (show_qi ? (i == 1 ? "Q$(q[1]), Q$(q[2])" : nothing) : nothing)
            markershape --> (show_qi ? :diamond : :circle)
            markercolor --> "#000000"
            markersize --> (show_qi ? 2 : 0)
            q_int, [h]
        end
        @series begin
            seriestype := :path
            label := nothing
            linecolor := "#000000"
            linewidth --> (show_qi ? 1.2 : 0)
            [qs[1], qs[2]], [h, h]
        end
        @series begin
            seriestype := :path
            label := (
                show_hpdi ? (i == 1 ? "$(Integer((1-hpdi[1])*100))% HPDI" : nothing) :
                nothing
            )
            linewidth --> (show_hpdi ? 2 : 0)
            seriesalpha --> 0.80
            linecolor --> :darkblue
            [lower_hpd[1][1], upper_hpd[1][1]], [h, h]
        end
    end
end

@recipe function f(
    p::ForestPlot;
    hpd_val = [0.05, 0.2],
    q = [0.1, 0.9],
    spacer = 0.5,
    _riser = 0.2,
    show_mean = true,
    show_median = true,
    show_qi = false,
    show_hpdi = true,
    fill_q = true,
    fill_hpd = false,
    ordered = false,
)

    chn = p.args[1]
    par_names = p.args[2]

    for i = 1:length(par_names)
        par,
        hpdi,
        lower_hpd,
        upper_hpd,
        h,
        qs,
        k_density,
        x_int,
        val,
        chain_med,
        chain_mean,
        min,
        q_int = _compute_plot_data(
            i,
            chn,
            par_names;
            hpd_val = hpd_val,
            q = q,
            spacer = spacer,
            _riser = _riser,
            show_mean = show_mean,
            show_median = show_median,
            show_qi = show_qi,
            show_hpdi = show_hpdi,
            fill_q = fill_q,
            fill_hpd = fill_hpd,
            ordered = ordered,
        )

        yticks --> (
            length(par_names) > 1 ?
            (_riser .+ ((1:length(par_names)) .- 1) .* spacer, string.(par)) : :default
        )
        yaxis --> (length(par_names) > 1 ? "Parameters" : "Density")

        for j = 1:length(hpdi)
            @series begin
                seriestype := :path
                label := (
                    show_hpdi ?
                    (i == 1 ? "$(Integer((1-hpdi[j])*100))% HPDI" : nothing) : nothing
                )
                linecolor --> j
                linewidth --> (show_hpdi ? 1.5 * j : 0)
                seriesalpha --> 0.80
                [lower_hpd[j][1], upper_hpd[j][1]], [h, h]
            end
        end
        @series begin
            seriestype := :scatter
            label := (show_median ? (i == 1 ? "Median" : nothing) : nothing)
            markershape --> :diamond
            markercolor --> "#000000"
            markersize --> (show_median ? length(hpdi) : 0)
            [chain_med], [h]
        end
        @series begin
            seriestype := :scatter
            label := (show_mean ? (i == 1 ? "Mean" : nothing) : nothing)
            markershape --> :circle
            markercolor --> :gray
            markersize --> (show_mean ? length(hpdi) : 0)
            [chain_mean], [h]
        end
        @series begin
            seriestype := :scatter
            label := (show_qi ? (i == 1 ? "Q1 = $(q[1]), Q3 = $(q[2])" : nothing) : nothing)
            markershape --> (show_qi ? :diamond : :circle)
            markercolor --> "#000000"
            markersize --> (show_qi ? 2 : 0)
            q_int, [h]
        end
        @series begin
            seriestype := :path
            label := nothing
            linecolor := "#000000"
            linewidth --> (show_qi ? 1.2 : 0.0)
            [qs[1], qs[2]], [h, h]
        end
    end
end
