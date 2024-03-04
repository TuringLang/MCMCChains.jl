# MCMCChains

Implementation of Julia types for summarizing MCMC simulations and utility functions for [diagnostics](@ref Diagnostics) and [visualizations](@ref StatsPlots.jl).

```@docs
ridgelineplot(chains::Chains, par_names::Vector{Symbol}; hpd_val = [0.05, 0.2],
q = [0.1, 0.9], spacer = 0.5, _riser = 0.2, show_mean = true, show_median = true,
show_qi = false, show_hpdi = true, fill_q = true, fill_hpd = false, ordered = false)
```

```@docs
forestplot(chains::Chains, par_names::Vector{Symbol}; hpd_val = [0.05, 0.2], q = [0.1, 0.9],
spacer = 0.5, _riser = 0.2, show_mean = true, show_median = true, show_qi = false,
show_hpdi = true, fill_q = true, fill_hpd = false, ordered = false)
```
