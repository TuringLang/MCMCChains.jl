"""
    CHAIN_PALETTE

Colours used to tell chains apart.

The Okabe-Ito palette, designed to stay distinguishable under protanopia, deuteranopia and
tritanopia. Blue, vermillion and green lead because the common case is a handful of chains,
and those three separate most strongly from each other.
"""
const CHAIN_PALETTE = [
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000",  # black
]

# Filled series stack on top of each other when several chains are drawn, so they need to be
# see-through to stay readable.
const FILL_ALPHA = 0.45

# Grey rather than black for anything that is not data, so the data carries the contrast.
const TEXT_COLOUR = "#262626"
const AXIS_COLOUR = "#545454"

const _DEFAULT_STYLE = (
    palette = CHAIN_PALETTE,
    grid = :none,
    background = nothing,
    framestyle = :axes,
    fontfamily = nothing,
    legend_panel = :first,
)

const _STYLE = Ref{NamedTuple}(_DEFAULT_STYLE)

const _GRID_KINDS = (:none, :lines, :dots)
const _LEGEND_PANELS = (:first, :all, :none)

"""
    plot_style()

The look MCMCChains gives its plots.

See [`plot_style!`](@ref) to change it.
"""
plot_style() = _STYLE[]

"""
    plot_style!(; kwargs...)

Change the look MCMCChains gives its plots, and return the new style.

Everything set here is a default, so a keyword passed to a single plot still wins.

# Keywords

- `palette`: the colours chains are drawn in, in order. Defaults to [`MCMCChains.CHAIN_PALETTE`](@ref).
- `grid`: `:none`, `:lines` or `:dots`.
- `background`: the colour behind the plot, `:transparent` for none, or `nothing` to leave it to Plots.
- `framestyle`: how much of the frame to draw, such as `:axes`, `:box`, `:origin`, `:grid` or `:none`.
- `fontfamily`: the font text is drawn in, or `nothing` for the Plots default.
- `legend_panel`: which panel of a multi-panel plot names the chains, one of `:first`, `:all` or `:none`.

# Examples

```julia
MCMCChains.plot_style!(grid = :dots, background = :transparent)
MCMCChains.plot_style!(palette = ["#332288", "#117733", "#DDCC77", "#CC6677"])
MCMCChains.reset_plot_style!()
```
"""
function plot_style!(; kwargs...)
    style = _STYLE[]
    for key in keys(kwargs)
        haskey(style, key) || throw(
            ArgumentError(
                "`$key` is not a style, expected one of $(join(keys(style), ", "))",
            ),
        )
    end

    grid = get(kwargs, :grid, style.grid)
    grid in _GRID_KINDS ||
        throw(ArgumentError("`grid` must be one of $(join(_GRID_KINDS, ", "))"))

    legend_panel = get(kwargs, :legend_panel, style.legend_panel)
    legend_panel in _LEGEND_PANELS ||
        throw(ArgumentError("`legend_panel` must be one of $(join(_LEGEND_PANELS, ", "))"))

    palette = get(kwargs, :palette, style.palette)
    isempty(palette) && throw(ArgumentError("`palette` needs at least one colour"))

    _STYLE[] = merge(style, NamedTuple(kwargs))
    return plot_style()
end

"""
    reset_plot_style!()

Put the look MCMCChains gives its plots back to the shipped default.
"""
function reset_plot_style!()
    _STYLE[] = _DEFAULT_STYLE
    return plot_style()
end

"""
    _chain_colours(nchains)

The first `nchains` palette entries as a row vector, which is how Plots reads one colour per
column of a matrix of series.
"""
function _chain_colours(nchains)
    palette = plot_style().palette
    return permutedims([palette[mod1(i, length(palette))] for i = 1:nchains])
end
