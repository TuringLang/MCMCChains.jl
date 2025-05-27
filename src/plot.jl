@shorthands meanplot
@shorthands autocorplot
@shorthands mixeddensity
@shorthands pooleddensity
@shorthands traceplot
@shorthands corner
@shorthands violinplot

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
    # Passes `_actual_indices_for_violin`` to `_ViolinPlot`` to ensure correct x-axis labels (parameter names or chain numbers).
    elseif st ∈ (:violinplot, :violin)
        show_boxplot_kw = get(plotattributes, :show_boxplot, true)
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
