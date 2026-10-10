# MCMCChains Changelog

## 7.8.0

Added `rankplot`.
Plots drop the grid, the top and right spines, and the legend box, and draw chains from the Okabe-Ito palette.
An attribute passed to a plot is no longer overruled by that look, which was dropping `grid`, `framestyle` and the font sizes.
`rankplot` draws chain 1 in the same colour as every other plot does.
A histogram bins every chain of a parameter over the same edges, so the bars line up between chains.
Histograms choose their bin count from the draws, by Freedman and Diaconis against Sturges, rather than always using 25.
`violinplot` gained a y axis label.
`ridgelineplot` and `forestplot` now plot every parameter when no parameter list is given.
Multi-panel plots reserve room for the y axis label, which was being clipped.

## 7.7.0

Remove support for PrettyTables.jl versions prior to 3.0.

## 7.6.0

Compatibility for PrettyTables@3.

Minimum Julia version bumped to 1.10.

## 7.5.0

Add a method for `MCMCDiagnosticTools.bfmi(::Chains)`. This computes the Bayesian Fraction of Missing Information for a chain or set of chains. Previously one had to extract a raw `Array` from the `Chains` object and pass that to `bfmi`.
