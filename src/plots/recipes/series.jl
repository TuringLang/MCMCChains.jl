@recipe function f(p::_DensityPlot)
    xaxis --> "Sample value"
    yaxis --> "Density"
    color_palette --> CHAIN_PALETTE
    trim --> true
    [collect(skipmissing(p.val[:, k])) for k = 1:size(p.val, 2)]
end

@recipe function f(p::_HistogramPlot)
    xaxis --> "Sample value"
    yaxis --> "Frequency"
    color_palette --> CHAIN_PALETTE
    fillalpha --> FILL_ALPHA
    linealpha --> 0.8
    # Freedman and Diaconis against Sturges, rather than a fixed count that is too coarse on
    # a long chain and too fine on a short one.
    bins --> :auto
    trim --> true

    series = [collect(skipmissing(p.val[:, k])) for k = 1:size(p.val, 2)]
    # Each chain would otherwise be binned over its own range, which gives bars of different
    # widths between chains and counts that cannot be compared.
    bins := _histogram_edges(series, get(plotattributes, :bins, :auto))
    series
end

@recipe function f(p::_MeanPlot)
    seriestype := :path
    color_palette --> CHAIN_PALETTE
    xaxis --> "Iteration"
    yaxis --> "Mean"
    range(p.c), cummean(p.val)
end

@recipe function f(p::_AutocorPlot)
    seriestype := :path
    color_palette --> CHAIN_PALETTE
    xaxis --> "Lag"
    yaxis --> "Autocorrelation"
    p.lags, p.val
end

@recipe function f(p::_TracePlot)
    seriestype := :path
    color_palette --> CHAIN_PALETTE
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
    color_palette --> CHAIN_PALETTE
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

    seriestype := :step
    xaxis --> "Rank (pooled over chains)"
    yaxis --> "Count"
    # Named rather than taken from the palette in turn, so that the reference line below
    # does not shift chain 1 off the colour it has in every other plot.
    linecolor --> _chain_colours(size(counts, 2))
    # A flat line at the expected count is the reference the eye compares against.
    @series begin
        seriestype := :hline
        label := nothing
        linecolor := AXIS_COLOUR
        linestyle := :dash
        linewidth := 1
        [size(p.val, 1) / nbins]
    end
    centres, counts
end
