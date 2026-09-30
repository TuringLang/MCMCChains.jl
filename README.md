# MCMCChains.jl

[![CI](https://github.com/TuringLang/MCMCChains.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/TuringLang/MCMCChains.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/TuringLang/MCMCChains.jl/branch/main/graph/badge.svg?token=TFxRFbKONS)](https://codecov.io/gh/TuringLang/MCMCChains.jl)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://TuringLang.github.io/MCMCChains.jl/stable/)
[![Dev Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://TuringLang.github.io/MCMCChains.jl/dev/)

Implementation of Julia types for summarizing MCMC simulations and utility functions for diagnostics and visualizations.

## Example

```julia
using MCMCChains
using StatsPlots

val = randn(100, 3, 2) .+ [1, 2, 3]'
val = hcat(val, rand(1:2, 100, 1, 2))
chn = Chains(val, [:A, :B, :C, :D])

plot(chn; size=(840, 600))
```

![Basic plot for Chains](https://turinglang.github.io/MCMCChains.jl/dev/default_plot.svg)

See the [docs](https://TuringLang.github.io/MCMCChains.jl/dev/) for more information.

## License Notice

Note that this package heavily uses and adapts code from the Mamba.jl package licensed under MIT License, see `LICENSE.md`.
