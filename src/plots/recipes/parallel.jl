"""
    _polyline(val, rows, npar)

One `x`, `y` pair tracing every row of `val` listed in `rows`, separated by `NaN`.

Handing Plots a matrix would make one series per draw, and a few thousand series is both
slow to draw and a few thousand legend entries. A break in the line costs nothing.
"""
function _polyline(val, rows, npar)
    n = length(rows)
    xs = Vector{Float64}(undef, n * (npar + 1))
    ys = Vector{Float64}(undef, n * (npar + 1))
    for (i, r) in enumerate(rows)
        at = (i - 1) * (npar + 1)
        for j = 1:npar
            xs[at+j] = j
            ys[at+j] = val[r, j]
        end
        xs[at+npar+1] = NaN
        ys[at+npar+1] = NaN
    end
    return xs, ys
end

@recipe function f(
    p::ParallelPlot;
    standardise = true,
    num_draws = nothing,
    random_seed = nothing,
)
    chn = _chain_arg(p, "parallelplot")
    _apply_chrome!(plotattributes)

    par_names, draws = _draws_by_parameter(chn)
    npar = length(par_names)
    npar > 1 || throw(ArgumentError("parallelplot needs at least two parameters"))

    val = standardise ? _standardise(draws) : float(draws)
    divergent = _divergences(chn)
    rng =
        random_seed === nothing ? Random.default_rng() : Random.MersenneTwister(random_seed)
    rows = _thin_draws(size(val, 1), num_draws, rng)

    xticks --> (1:npar, string.(par_names))
    xlims --> (0.5, npar + 0.5)
    xaxis --> "Parameters"
    yaxis --> (standardise ? "Standardised value" : "Sample value")
    legend --> (divergent === nothing ? false : :outertopright)

    ordinary = divergent === nothing ? rows : filter(r -> !divergent[r], rows)
    flagged = divergent === nothing ? Int[] : filter(r -> divergent[r], rows)
    palette = plot_style().palette

    if !isempty(ordinary)
        @series begin
            seriestype := :path
            label := divergent === nothing ? nothing : "Draws"
            linecolor := palette[1]
            # Every line at full strength is a solid block, so the bundle fades as it grows.
            linealpha --> clamp(20 / length(ordinary), 0.02, 1.0)
            linewidth --> 1
            _polyline(val, ordinary, npar)
        end
    end

    if !isempty(flagged)
        @series begin
            seriestype := :path
            label := "Divergent"
            linecolor := palette[min(2, length(palette))]
            linealpha --> 0.7
            linewidth --> 1.2
            _polyline(val, flagged, npar)
        end
    end
end
