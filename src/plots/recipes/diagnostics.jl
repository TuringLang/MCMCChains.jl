@recipe function f(p::_ParameterDots)
    n = length(p.par_names)

    yticks --> (1:n, string.(p.par_names))
    ylims --> (0.5, n + 0.5)
    yaxis --> "Parameters"
    xaxis --> p.label
    legend --> false

    for r in p.references
        @series begin
            seriestype := :vline
            label := nothing
            linecolor := AXIS_COLOUR
            linestyle := :dash
            linewidth := 1
            [r]
        end
    end

    seriestype := :scatter
    markercolor --> first(plot_style().palette)
    markerstrokewidth --> 0
    markersize --> 5
    p.val, collect(1:n)
end

@recipe function f(p::EssPlot; kind = :bulk, relative = false)
    chn = _chain_arg(p, "essplot")
    _apply_chrome!(plotattributes)

    df = ess(chn; kind = kind, relative = relative)
    niter, _, nchains = size(Chains(chn, :parameters))
    # Vehtari et al. (2021) ask for an effective sample size of 400 before reading R-hat.
    threshold = relative ? 400 / (niter * nchains) : 400.0
    label =
        relative ? "Effective sample size / draws ($kind)" : "Effective sample size ($kind)"

    _ParameterDots(df.nt.parameters, collect(df.nt.ess), [threshold], label)
end

@recipe function f(p::RhatPlot; kind = :rank)
    chn = _chain_arg(p, "rhatplot")
    _apply_chrome!(plotattributes)

    df = rhat(chn; kind = kind)
    # Vehtari et al. (2021) recommend sampling until R-hat is below 1.01.
    _ParameterDots(df.nt.parameters, collect(df.nt.rhat), [1.01], "R-hat ($kind)")
end

# Does not apply, bad, worth a look, fine. Vermillion and bluish green carry the meaning
# without leaning on the red and green that colour blindness confuses.
const _SCORE_COLOURS = Dict(-1 => "#9A9A9A", 0 => "#D55E00", 1 => "#E69F00", 2 => "#009E73")

@recipe function f(p::DiagnosticsPlot)
    chn = _chain_arg(p, "diagnosticsplot")
    _apply_chrome!(plotattributes)

    par_names, column_names, values, scores = _diagnostics_table(chn)
    nrows, ncols = size(values)

    legend := false
    grid := false
    # The cells are the frame, so the only thing the axes still owe us is the labels.
    framestyle := :grid
    foreground_color_axis := :transparent
    foreground_color_border := :transparent
    xticks := (1:ncols, column_names)
    yticks := (1:nrows, string.(par_names))
    xlims --> (0.4, ncols + 0.6)
    ylims --> (0.4, nrows + 0.6)
    # Rows read top to bottom, the way a table does.
    yflip := true
    xmirror := true
    size --> (150 * ncols + 80, 34 * nrows + 90)

    for level = -1:2
        cells = [(i, j) for i = 1:nrows, j = 1:ncols if scores[i, j] == level]
        isempty(cells) && continue
        xs = Float64[]
        ys = Float64[]
        for (i, j) in cells
            append!(xs, [j - 0.47, j + 0.47, j + 0.47, j - 0.47, NaN])
            append!(ys, [i - 0.42, i - 0.42, i + 0.42, i + 0.42, NaN])
        end
        @series begin
            seriestype := :shape
            label := nothing
            linewidth := 0
            fillcolor := _SCORE_COLOURS[level]
            xs, ys
        end
    end

    annotations := [
        (
            j,
            i,
            (_format_cell(values[i, j], _DIAGNOSTIC_COLUMNS[j].digits), 9, :white, :center),
        ) for i = 1:nrows for j = 1:ncols
    ]

    ()
end

@recipe function f(p::EvolutionPlot; diagnostic = :ess, npoints = 20)
    chn = _chain_arg(p, "evolutionplot")
    _apply_chrome!(plotattributes)

    points, par_names, values = _evolution(chn, diagnostic, npoints)

    xaxis --> "Draws per chain"
    yaxis --> (diagnostic === :ess ? "Effective sample size" : "R-hat")
    label --> permutedims(string.(par_names))
    linewidth --> 1.5
    markershape --> :circle
    markersize --> 3
    markerstrokewidth --> 0

    reference = diagnostic === :ess ? 400.0 : 1.01
    @series begin
        seriestype := :hline
        label := nothing
        linecolor := AXIS_COLOUR
        linestyle := :dash
        linewidth := 1
        [reference]
    end

    points, values
end

@recipe function f(p::McsePlot; relative = false)
    chn = _chain_arg(p, "mcseplot")
    _apply_chrome!(plotattributes)

    df = mcse(chn)
    val = collect(df.nt.mcse)
    label = "Monte Carlo standard error"
    if relative
        _, draws = _draws_by_parameter(chn)
        val = val ./ vec(std(draws; dims = 1))
        label = "Monte Carlo standard error / SD"
    end

    _ParameterDots(df.nt.parameters, val, Float64[], label)
end
