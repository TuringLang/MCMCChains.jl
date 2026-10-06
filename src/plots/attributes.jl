"""
    _interval_ylims(riser, spacer, nparams, top)

Vertical limits for ridgeline and forest plots.

Plots sizes the axis to the data, so the outermost rows end up flush against the frame and
their tick labels are clipped. `top` is the highest drawn value, which for a forest plot is
the last baseline and for a ridgeline plot is the tallest density.
"""
function _interval_ylims(riser, spacer, nparams, top)
    pad = spacer / 2
    return (riser - pad, max(top, riser + (nparams - 1) * spacer) + pad)
end

"""
    CHAIN_PALETTE

Colours used to distinguish chains.

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

"""
    _axis_default!(plotattributes, key, value)

Set an axis attribute, unless the caller already set it.

Plots expands an axis attribute such as `grid` into `xgrid`, `ygrid` and `zgrid` before a
recipe runs, and does not always leave the original behind, so a plain `get!` on `grid`
finds nothing and overrules the caller.
"""
function _axis_default!(plotattributes, key::Symbol, value)
    for prefix in ("", "x", "y", "z")
        haskey(plotattributes, Symbol(prefix, key)) && return nothing
    end
    plotattributes[key] = value
    return nothing
end

"""
    _apply_chrome!(plotattributes)

Apply the shared look to everything that is not data.

Gridlines, a full frame and a legend box are ink that encodes nothing, and on a diagnostic
plot they compete with the thing being judged, so they are off by default. Text and axes are
grey rather than black for the same reason. Each setting is a default, so anything the
caller passes wins.
"""
function _apply_chrome!(plotattributes)
    _axis_default!(plotattributes, :grid, false)
    _axis_default!(plotattributes, :framestyle, :axes)
    get!(plotattributes, :foreground_color_legend, nothing)
    _axis_default!(plotattributes, :foreground_color_axis, AXIS_COLOUR)
    _axis_default!(plotattributes, :foreground_color_border, AXIS_COLOUR)
    _axis_default!(plotattributes, :foreground_color_guide, TEXT_COLOUR)
    get!(plotattributes, :foreground_color_text, TEXT_COLOUR)
    get!(plotattributes, :titlefontsize, 12)
    _axis_default!(plotattributes, :guidefontsize, 10)
    _axis_default!(plotattributes, :tickfontsize, 8)
    get!(plotattributes, :legendfontsize, 8)
    get!(plotattributes, :left_margin, (8, :mm))
    get!(plotattributes, :bottom_margin, (3, :mm))
    return nothing
end
