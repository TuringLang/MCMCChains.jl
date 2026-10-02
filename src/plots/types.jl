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
