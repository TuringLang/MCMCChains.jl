@shorthands meanplot
@shorthands autocorplot
@shorthands mixeddensity
@shorthands pooleddensity
@shorthands traceplot
@shorthands corner
@userplot RidgelinePlot
@userplot ForestPlot

struct _TracePlot; c; val; end
struct _MeanPlot; c; val;  end
struct _DensityPlot; c; val;  end
struct _HistogramPlot; c; val;  end
struct _AutocorPlot; lags; val;  end

# define alias functions for old syntax
const translationdict = Dict(
                        :traceplot => _TracePlot,
                        :meanplot => _MeanPlot,
                        :density => _DensityPlot,
                        :histogram => _HistogramPlot,
                        :autocorplot => _AutocorPlot,
                        :pooleddensity => _DensityPlot
                      )

const supportedplots = push!(collect(keys(translationdict)), :mixeddensity, :corner)

@recipe f(c::Chains, s::Symbol) = c, [s]

@recipe function f(
    chains::Chains, i::Int;
    colordim = :chain,
    barbounds = (-Inf, Inf),
    maxlag = nothing,
    append_chains = false
)
    st = get(plotattributes, :seriestype, :traceplot)
    c = append_chains || st == :pooleddensity ? pool_chain(chains) : chains

    if colordim == :parameter
        title --> "Chain $(MCMCChains.chains(c)[i])"
        label --> permutedims(map(string, names(c)))
        val = c.value[:, :, i]
    elseif colordim == :chain
        title --> string(names(c)[i])
        label --> permutedims(map(x -> "Chain $x", MCMCChains.chains(c)))
        val = c.value[:, i, :]
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
        ac = autocor(c; sections = nothing, lags = lags, append_chains=false)
        ac_mat = cat(reduce.(hcat, Iterators.drop.(ac, 1))...; dims=3)
        val = colordim == :parameter ? ac_mat[:, :, i]' : ac_mat[i, :, :]
        _AutocorPlot(lags, val)
    elseif st ∈ supportedplots
        translationdict[st](c, val)
    else
        range(c), val
    end
end

@recipe function f(p::_DensityPlot)
    xaxis --> "Sample value"
    yaxis --> "Density"
    trim --> true
    [collect(skipmissing(p.val[:,k])) for k in 1:size(p.val, 2)]
end

@recipe function f(p::_HistogramPlot)
    xaxis --> "Sample value"
    yaxis --> "Frequency"
    fillalpha --> 0.7
    bins --> 25
    trim --> true
    [collect(skipmissing(p.val[:,k])) for k in 1:size(p.val, 2)]
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

@recipe function f(
    chains::Chains,
    parameters::AbstractVector{Symbol};
    colordim = :chain
)
    colordim != :chain &&
        error("Symbol names are interpreted as parameter names, only compatible with ",
              "`colordim = :chain`")

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
    append_chains = false
)
    _chains = isempty(parameters) ? Chains(chains, _clean_sections(chains, sections)) : chains
    c = append_chains ? pool_chain(_chains) : _chains
    ptypes = get(plotattributes, :seriestype, (:traceplot, :mixeddensity))
    ptypes = ptypes isa Symbol ? (ptypes,) : ptypes
    @assert all(ptype -> ptype ∈ supportedplots, ptypes)
    ntypes = length(ptypes)
    nrows, nvars, nchains = size(c)
    isempty(parameters) && (parameters = colordim == :chain ? (1:nvars) : (1:nchains))
    N = length(parameters)

    if :corner ∉ ptypes
        size --> (ntypes*width, N*height)
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
    c
    parameters
end

@recipe function f(corner::Corner)
    # Convert labels to string because `Symbol` is not supported generally supported.
    label --> permutedims(map(string, corner.parameters))
    compact --> true
    size --> (600, 600)
    # NOTE: Don't use the indices from `chains(chains)`.
    # See https://github.com/TuringLang/MCMCChains.jl/issues/413.
    ar = collect(Array(corner.c.value[:, corner.parameters, i]) for i in 1:length(chains(corner.c)))
    RecipesBase.recipetype(:cornerplot, vcat(ar...))
end

function _compute_plot_data(
    i::Integer,
    chains::Chains,
    par_names::AbstractVector{Symbol};
    hdi_prob = [0.94, 0.8],
    q = [0.1, 0.9],
    spacer = 0.4,
    _riser = 0.2,
    barbounds = (-Inf, Inf),
    show_mean = true,
    show_median = true,
    show_qi = false,
    show_hdii = true,
    fill_q = true,
    fill_hdi = false,
    ordered = false
)

    chain_dic = Dict(zip(quantile(chains)[2], quantile(chains)[5]))
    sorted_chain = sort(collect(zip(values(chain_dic), keys(chain_dic))))
    sorted_par = [sorted_chain[i][2] for i in 1:length(par_names)]
    par = (ordered ? sorted_par : par_names)
    hdii = sort(hdi_prob)

    chain_sections = MCMCChains.group(chains, Symbol(par[i]))
    chain_vec = vec(chain_sections.value.data)
    lower_hdi = [MCMCChains.hdi(chain_sections, prob = hdii[j])[:lower]
        for j in 1:length(hdii)]
    upper_hdi = [MCMCChains.hdi(chain_sections, prob = hdii[j])[:upper]
        for j in 1:length(hdii)]
    h = _riser + spacer*(i-1)
    qs = quantile(chain_vec, q)
    k_density = kde(chain_vec)
    if fill_hdi
        x_int = filter(x -> lower_hdi[1][1] <= x <= upper_hdi[1][1], k_density.x)
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

    return par, hdii, lower_hdi, upper_hdi, h, qs, k_density, x_int, val, chain_med,
        chain_mean, min, q_int
end

@recipe function f(
    p::RidgelinePlot;
    hdi_prob = [0.94, 0.8],
    q = [0.1, 0.9],
    spacer = 0.5,
    _riser = 0.2,
    show_mean = true,
    show_median = true,
    show_qi = false,
    show_hdii = true,
    fill_q = true,
    fill_hdi = false,
    ordered = false
)

    chn = p.args[1]
    par_names = p.args[2]

    for i in 1:length(par_names)
        par, hdii, lower_hdi, upper_hdi, h, qs, k_density, x_int, val, chain_med, chain_mean,
            min, q_int = _compute_plot_data(i, chn, par_names; hdi_prob = hdi_prob, q = q,
            spacer = spacer, _riser = _riser, show_mean = show_mean, show_median = show_median,
            show_qi = show_qi, show_hdii = show_hdii, fill_q = fill_q, fill_hdi = fill_hdi,
            ordered = ordered)

        yticks --> (length(par_names) > 1 ?
            (_riser .+ ((1:length(par_names)) .- 1) .* spacer, string.(par)) : :default)
        yaxis --> (length(par_names) > 1 ? "Parameters" : "Density" )
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
            label := (show_hdii ? (i == 1 ? "$(round(Int, hdii[i]*100))% HDI" : nothing)
                : nothing)
            linewidth --> (show_hdii ? 2 : 0)
            seriesalpha --> 0.80
            linecolor --> :darkblue
            [lower_hdi[1][1], upper_hdi[1][1]], [h, h]
        end
    end
end

@recipe function f(
    p::ForestPlot;
    hdi_prob = [0.94, 0.8],
    q = [0.1, 0.9],
    spacer = 0.5,
    _riser = 0.2,
    show_mean = true,
    show_median = true,
    show_qi = false,
    show_hdii = true,
    fill_q = true,
    fill_hdi = false,
    ordered = false
)

    chn = p.args[1]
    par_names = p.args[2]

    for i in 1:length(par_names)
        par, hdii, lower_hdi, upper_hdi, h, qs, k_density, x_int, val, chain_med, chain_mean,
            min, q_int = _compute_plot_data(i, chn, par_names; hdi_prob = hdi_prob, q = q,
            spacer = spacer, _riser = _riser, show_mean = show_mean, show_median = show_median,
            show_qi = show_qi, show_hdii = show_hdii, fill_q = fill_q, fill_hdi = fill_hdi,
            ordered = ordered)

        yticks --> (length(par_names) > 1 ?
            (_riser .+ ((1:length(par_names)) .- 1) .* spacer, string.(par)) : :default)
        yaxis --> (length(par_names) > 1 ? "Parameters" : "Density" )

        for j in 1:length(hdii)
            @series begin
                seriestype := :path
                label := (show_hdii ?
                    (i == 1 ? "$(round(Int, hdii[j]*100))% HDI" : nothing) : nothing)
                linecolor --> j
                linewidth --> (show_hdii ? 1.5*j : 0)
                seriesalpha --> 0.80
                [lower_hdi[j][1], upper_hdi[j][1]], [h, h]
            end
        end
        @series begin
            seriestype := :scatter
            label := (show_median ? (i == 1 ? "Median" : nothing) : nothing)
            markershape --> :diamond
            markercolor --> "#000000"
            markersize --> (show_median ? length(hdii) : 0)
            [chain_med], [h]
        end
        @series begin
            seriestype := :scatter
            label := (show_mean ? (i == 1 ? "Mean" : nothing) : nothing)
            markershape --> :circle
            markercolor --> :gray
            markersize --> (show_mean ? length(hdii) : 0)
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
