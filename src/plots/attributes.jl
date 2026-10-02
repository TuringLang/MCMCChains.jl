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
