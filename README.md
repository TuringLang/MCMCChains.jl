# MCMCChains.jl

![CI](https://github.com/TuringLang/MCMCChains.jl/workflows/CI/badge.svg?branch=main)
[![codecov](https://codecov.io/gh/TuringLang/MCMCChains.jl/branch/main/graph/badge.svg?token=TFxRFbKONS)](https://codecov.io/gh/TuringLang/MCMCChains.jl)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://TuringLang.github.io/MCMCChains.jl/stable/)
[![Dev Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://TuringLang.github.io/MCMCChains.jl/dev/)

Implementation of Julia types for summarizing MCMC simulations and utility functions for diagnostics and visualizations.

```
using MCMCChains
using StatsPlots

plot(chn)
```
![Basic plot for Chains](https://turinglang.github.io/MCMCChains.jl/dev/default_plot.svg)

See the [docs](https://TuringLang.github.io/MCMCChains.jl/dev/) for more information.

## License Notice

Note that this package heavily uses and adapts code from the Mamba.jl package licensed under MIT License, see `LICENSE.md`.
