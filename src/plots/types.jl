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
    # colordim: the series are chains of one parameter, or parameters of one chain. Only
    # chains of one parameter can share bin edges.
    colordim::Any
end
_HistogramPlot(c, val) = _HistogramPlot(c, val, :chain)
struct _RankPlot
    c::Any
    val::Any
end
struct _EcdfPlot
    c::Any
    val::Any
end
struct _CountPlot{C,V}
    c::C
    val::V
end

# Plots treats unwrapped numbers as coordinates.
struct _CountValue{T}
    value::T
end
Base.isequal(a::_CountValue, b::_CountValue) = isequal(a.value, b.value)
Base.hash(v::_CountValue, h::UInt) = hash(v.value, h)
Base.show(io::IO, v::_CountValue) = print(io, v.value)
struct _AutocorPlot
    lags::Any
    val::Any
end
# One dot per parameter, used by every per-parameter diagnostic.
struct _ParameterDots
    par_names::Any
    val::Any
    # references: values to draw a dashed line at, such as the R-hat threshold.
    references::Any
    label::Any
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
    :countplot => _CountPlot,
    :autocorplot => _AutocorPlot,
    :rankplot => _RankPlot,
    :ecdfplot => _EcdfPlot,
    :pooleddensity => _DensityPlot,
    :violinplot => _ViolinPlot,
    :violin => _ViolinPlot,
)

const supportedplots = push!(collect(keys(translationdict)), :mixeddensity)
