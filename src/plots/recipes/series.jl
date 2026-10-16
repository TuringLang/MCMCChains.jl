@recipe function f(p::_DensityPlot)
    xaxis --> "Sample value"
    yaxis --> "Density"
    trim --> true
    [collect(skipmissing(p.val[:, k])) for k = 1:size(p.val, 2)]
end

@recipe function f(p::_HistogramPlot)
    xaxis --> "Sample value"
    yaxis --> "Frequency"
    fillalpha --> FILL_ALPHA
    linealpha --> 0.8
    # Read off the draws, rather than a fixed count that is too coarse on a long chain and
    # too fine on a short one.
    bins --> :auto
    trim --> true

    series = [collect(skipmissing(p.val[:, k])) for k = 1:size(p.val, 2)]
    # Chains of one parameter would otherwise each be binned over their own range, which
    # gives bars of different widths and counts that cannot be compared. Parameters are a
    # different matter: two of them may be on scales that share no useful edges at all.
    if p.colordim === :chain
        bins := _histogram_edges(series, get(plotattributes, :bins, :auto))
    end
    series
end

@recipe function f(p::_EcdfPlot)
    seriestype := :steppost
    xaxis --> "Sample value"
    yaxis --> "ECDF"
    xs, ys = _ecdf_series(p.val)
    xs, ys
end

@recipe function f(p::_CountPlot; relative = false)
    levels, counts = _value_counts(p.val)
    heights = if relative
        totals = sum(counts; dims = 1)
        all(>(0), totals) ||
            throw(ArgumentError("countplot proportions require draws in every series"))
        counts ./ totals
    else
        counts
    end

    xguide --> "Sample value (categorical)"
    yguide --> (relative ? "Proportion" : "Count")
    nchains = size(heights, 2)
    width = something(get(plotattributes, :bar_width, nothing), 0.8)
    width isa Real && isfinite(width) && width > 0 ||
        throw(ArgumentError("countplot requires a finite, positive bar_width"))
    width /= nchains
    for c = 1:nchains
        @series begin
            seriestype := :_countbar
            bar_width := width
            countplot_values := _CountValue.(levels)
            countplot_offset := (c - (nchains + 1) / 2) .* width
            collect(eachindex(levels)), heights[:, c]
        end
    end
    nothing
end

function RecipesBase.apply_recipe(
    attrs::AbstractDict{Symbol,Any},
    ::Type{Val{:_countbar}},
    x,
    y,
    z,
)
    values, offset = attrs[:countplot_values], attrs[:countplot_offset]
    delete!(attrs, :countplot_values)
    delete!(attrs, :countplot_offset)
    axis = get(attrs, :orientation, :vertical) in (:vertical, :v, :vert) ? :x : :y
    # The bar recipe registers exact values on the subplot's categorical axis.
    series = RecipesBase.apply_recipe(attrs, Val{:bar}, values, y, z)
    for data in series
        data.plotattributes[axis] = data.plotattributes[axis] .+ offset
    end
    return series
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
    yaxis --> "Sample value"
    fillalpha --> FILL_ALPHA
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

@recipe function f(p::_RankPlot; nbins = 20)
    edges, counts = _rank_bin_counts(p.val, nbins)
    centres = [(edges[i] + edges[i+1]) / 2 for i = 1:(length(edges)-1)]
    nchains = size(counts, 2)
    colours = _chain_colours(nchains)
    labels = get(plotattributes, :label, nothing)

    xaxis --> "Rank (pooled over chains)"
    yaxis --> "Count"

    # Each chain is named and coloured on its own series. Plots hands a row vector of
    # colours or labels out one series at a time, so the reference line would take the
    # first of each and leave every chain wearing the one beside it.
    for c = 1:nchains
        @series begin
            seriestype := :step
            linecolor := colours[c]
            label := labels isa AbstractArray ? labels[min(c, length(labels))] : nothing
            centres, counts[:, c]
        end
    end

    # A flat line at the expected count is the reference the eye compares against. Returned
    # rather than emitted, because a recipe that ends without data still adds a series, and
    # that series takes the next label off the row vector.
    seriestype := :hline
    label := nothing
    linecolor := AXIS_COLOUR
    linestyle := :dash
    linewidth := 1
    [size(p.val, 1) / nbins]
end
