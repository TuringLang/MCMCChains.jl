# MCMCChains Changelog

## 8.0.0

### Added

- `describe` accepts `digits` to set the decimal places in its summary and quantile tables.
- `hcat` now adds scalar or vector values by name, including matrices of named tuples, while preserving iteration labels and chain metadata.
- `var` and `std`, with one estimate per parameter.
- `cov`, with parameter names on both axes of the covariance matrix.

### Changed

- `corner`, `corner!` and `seriestype=:corner` were replaced by `StatsPlots.cornerplot` and `StatsPlots.cornerplot!`.

### Fixed

- `StatsPlots.cornerplot` and `StatsPlots.cornerplot!` threw when given a `Chains` object.
- `hcat` failed to add named-tuple draws stored in an `Any`-typed matrix and could return a plain matrix for one draw.
- `ChainDataFrame` summaries could not be converted with `Array`, and `convert(Array, summary)` failed for a single statistic.
- `get` returned numeric array parameters in stored column order instead of array index order.
- `Chains` selection by a vector of parameter names returned columns in stored order instead of requested order.
- `Chains` assignment by a vector of parameter names placed values in stored order instead of requested order.
- `get` and `get_params` with `flatten=true` lost draws, while `get` also threw for section selections or included names from other sections.

## 7.8.1

- `countplot`, which counts distinct sampled values on a categorical axis.
- `histogram` placed empty bars between a parameter's few sampled values and now bins on those values under `bins = :auto` when the edges fit.

## 7.8.0

### Added

- `diagnosticsplot`, which puts every convergence diagnostic for every parameter in one coloured grid.
- `essplot`, `rhatplot` and `mcseplot`, one dot per parameter for the diagnostics the package already reports.
- `evolutionplot`, which recomputes a diagnostic on growing prefixes of the chain.
- `parallelplot`, which draws every draw across all parameters and picks out the divergent ones.
- `nutsplot`, for the acceptance rate, step size, tree depth and divergence count an HMC sampler records.
- `rankplot`.
- The `:ecdfplot` series type.
- `plot_style!` and `reset_plot_style!`, which set the palette, the grid, the background, the frame, the font, and where a multi-panel plot names its chains.

### Changed

- Plots drop the grid, the top and right spines, and the legend box, and draw chains from the Okabe-Ito palette.
- Multi-panel plots name the chains in one panel, instead of leaving the lines unlabelled.
- Histograms choose their bin count from the draws rather than always using 25.

### Fixed

- `ridgelineplot` and `forestplot` threw when given no parameter list, though the docstring said they would plot every parameter.
- `energyplot` and `ppcplot` were missing the shared look.
- `rankplot` drew chain 1 in a different colour from every other plot.
- Histograms binned each chain over its own range, so bars differed in width between chains.
- An attribute passed to a plot was overruled by the shared look, which dropped `grid`, `framestyle` and the font sizes.
- Multi-panel plots clipped the y axis label.
- `violinplot` had no y axis label.

## 7.7.0

Remove support for PrettyTables.jl versions prior to 3.0.

## 7.6.0

Compatibility for PrettyTables@3.

Minimum Julia version bumped to 1.10.

## 7.5.0

Add a method for `MCMCDiagnosticTools.bfmi(::Chains)`. This computes the Bayesian Fraction of Missing Information for a chain or set of chains. Previously one had to extract a raw `Array` from the `Chains` object and pass that to `bfmi`.
