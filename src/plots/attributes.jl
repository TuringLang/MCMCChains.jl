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
tritanopia.
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

# Grey rather than black for anything that is not data, so the data carries the contrast.
const AXIS_COLOUR = "#545454"
