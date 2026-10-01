# MCMCChains Changelog

## 7.8.0

`ridgelineplot` and `forestplot` now plot every parameter when no parameter list is given, which their docstrings already described.
Both now leave room for the outermost rows and draw the legend outside the axes, so a row is no longer drawn on the frame or hidden behind the legend.
Multi-panel plots such as `meanplot` and `mixeddensity` now reserve room for the y axis label, which Plots otherwise draws off the left edge of a tall figure.
Plots now share a single look: no grid, only the left and bottom spines, an unboxed legend, and grey rather than black chrome, so the ink that remains is the data.
Chains are drawn from the Okabe-Ito palette, which stays distinguishable under the common forms of colour blindness and keeps a chain the same colour in every panel.
`violinplot` gained the y axis label it was missing.
Added `rankplot`, which ranks every draw against every other draw and shows each chain's share of the rank range against the uniform count it should hit.
Vehtari et al. (2021) recommend it in place of the trace plot, which loses its diagnostic value once chains are long.

## 7.7.0

Remove support for PrettyTables.jl versions prior to 3.0.

## 7.6.0

Compatibility for PrettyTables@3.

Minimum Julia version bumped to 1.10.

## 7.5.0

Add a method for `MCMCDiagnosticTools.bfmi(::Chains)`. This computes the Bayesian Fraction of Missing Information for a chain or set of chains. Previously one had to extract a raw `Array` from the `Chains` object and pass that to `bfmi`.
