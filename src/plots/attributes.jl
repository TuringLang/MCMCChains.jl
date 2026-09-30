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

Taken from ArviZ's default style, which picked them for distinguishability under the common
forms of colour blindness.
"""
const CHAIN_PALETTE = [
    "#36ACC6",  # cyan
    "#F66D7F",  # rose
    "#FAC364",  # amber
    "#7C2695",  # purple
    "#228306",  # green
    "#A252F4",  # violet
    "#63F0EA",  # turquoise
    "#A77E4F",  # brown
]

# Filled series stack on top of each other when several chains are drawn, so they need to be
# see-through to stay readable.
const FILL_ALPHA = 0.45

# Grey rather than black for anything that is not data, so the data carries the contrast.
const TEXT_COLOUR = "#262626"
const AXIS_COLOUR = "#545454"

"""
    _apply_chrome!(plotattributes)

Apply the shared look to everything that is not data.

Follows ArviZ's default style: no grid, only the left and bottom spines, an unboxed legend,
and grey rather than black for text and axes, so the data carries the contrast. Each
setting is a default, so anything the caller passes wins.
"""
function _apply_chrome!(plotattributes)
    get!(plotattributes, :grid, false)
    get!(plotattributes, :framestyle, :axes)
    get!(plotattributes, :foreground_color_legend, nothing)
    get!(plotattributes, :foreground_color_axis, AXIS_COLOUR)
    get!(plotattributes, :foreground_color_border, AXIS_COLOUR)
    get!(plotattributes, :foreground_color_text, TEXT_COLOUR)
    get!(plotattributes, :foreground_color_guide, TEXT_COLOUR)
    get!(plotattributes, :titlefontsize, 12)
    get!(plotattributes, :guidefontsize, 10)
    get!(plotattributes, :tickfontsize, 8)
    get!(plotattributes, :legendfontsize, 8)
    get!(plotattributes, :left_margin, (8, :mm))
    get!(plotattributes, :bottom_margin, (3, :mm))
    return nothing
end
