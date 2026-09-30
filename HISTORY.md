# MCMCChains Changelog

## 7.8.0

`ridgelineplot` and `forestplot` now plot every parameter when no parameter list is given, which their docstrings already described.
Both now leave room for the outermost rows and draw the legend outside the axes, so a row is no longer drawn on the frame or hidden behind the legend.
Multi-panel plots such as `meanplot` and `mixeddensity` now reserve room for the y axis label, which Plots otherwise draws off the left edge of a tall figure.

## 7.7.0

Remove support for PrettyTables.jl versions prior to 3.0.

## 7.6.0

Compatibility for PrettyTables@3.

Minimum Julia version bumped to 1.10.

## 7.5.0

Add a method for `MCMCDiagnosticTools.bfmi(::Chains)`. This computes the Bayesian Fraction of Missing Information for a chain or set of chains. Previously one had to extract a raw `Array` from the `Chains` object and pass that to `bfmi`.
