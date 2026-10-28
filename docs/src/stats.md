# Posterior statistics

## Display precision

Use `digits` to choose the decimal places shown in `describe`'s summary and quantile tables:

```@example precision
using MCMCChains
chn = Chains(reshape(1:400, 100, 2, 2) ./ 3, ["mu", "theta[1]"])
describe(chn; digits = 2)
```

When `digits` is omitted, `describe` reads `:digits` from its `IOContext`. The same setting applies when displaying a summary table directly:

```@example precision
show(IOContext(stdout, :digits => 2), MIME"text/plain"(), mean(chn))
```

## Selecting parameters

Use sections to choose which parameter groups appear in the summary and quantile tables:

```@example precision
grouped = set_section(chn, Dict(:parameters => ["mu"], :vectors => ["theta[1]"]))
describe(grouped; digits = 2)
```

Use `sections = :vectors` for the vector entries or `sections = nothing` for all entries.

## API

The methods listed below are defined in `src/stats.jl`.

```@docs
autocor
describe
mean
MCMCChains.var
MCMCChains.std
MCMCChains.cov
summarystats
quantile
hpd
```
