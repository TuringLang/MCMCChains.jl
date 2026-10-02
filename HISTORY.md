# MCMCChains Changelog

## 7.8.0

Added `rankplot`.
`ridgelineplot` and `forestplot` now plot every parameter when no parameter list is given.
`violinplot` gained a y axis label.
Multi-panel plots reserve room for the y axis label, which was being clipped.
Plots drop the grid, the top and right spines, and the legend box, and draw chains from the Okabe-Ito palette.

## 7.7.0

Remove support for PrettyTables.jl versions prior to 3.0.

## 7.6.0

Compatibility for PrettyTables@3.

Minimum Julia version bumped to 1.10.

## 7.5.0

Add a method for `MCMCDiagnosticTools.bfmi(::Chains)`. This computes the Bayesian Fraction of Missing Information for a chain or set of chains. Previously one had to extract a raw `Array` from the `Chains` object and pass that to `bfmi`.
