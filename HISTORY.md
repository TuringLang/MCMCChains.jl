# MCMCChains Changelog

## 7.8.0

Added `rankplot`.
Added `diagnosticsplot`, which puts every convergence diagnostic for every parameter in one coloured grid.
Added `essplot`, `rhatplot` and `mcseplot`, which draw the convergence diagnostics MCMCChains already reports.
Added `evolutionplot`, which recomputes a diagnostic on growing prefixes of the chain.
Added `parallelplot`, which draws every draw across all parameters and picks out the divergent ones.
Added `nutsplot` for the acceptance rate, step size, tree depth and divergence count an HMC sampler records.
Added the `:ecdfplot` series type.
Added `plot_style!` and `reset_plot_style!` to set the palette, the grid, the background, the frame, the font and where a multi-panel plot names its chains.
Multi-panel plots name the chains in one panel, instead of leaving the lines unlabelled.
A histogram bins every chain of a parameter over the same edges, so the bars line up between chains.
Histograms choose their bin count from the draws, by Freedman and Diaconis against Sturges, rather than always using 25.
Plots drop the grid, the top and right spines, and the legend box, and draw chains from the Okabe-Ito palette.
`energyplot` and `ppcplot` take that look too, which they were missing.
`rankplot` draws chain 1 in the same colour as every other plot does.
An attribute passed to a plot is no longer overruled by that look, which was dropping `grid`, `framestyle` and the font sizes.
`violinplot` gained a y axis label.
`ridgelineplot` and `forestplot` now plot every parameter when no parameter list is given.
Multi-panel plots reserve room for the y axis label, which was being clipped.

## 7.7.0

Remove support for PrettyTables.jl versions prior to 3.0.

## 7.6.0

Compatibility for PrettyTables@3.

Minimum Julia version bumped to 1.10.
