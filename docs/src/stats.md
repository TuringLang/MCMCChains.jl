# Posterior statistics

## Display precision

Set `:digits` in an `IOContext` to choose the number of decimal places printed by `describe`:

```@example precision
using MCMCChains
chn = Chains(reshape(1:400, 100, 2, 2) ./ 3, [:a, :b])
describe(IOContext(stdout, :digits => 2), chn)
```

The same setting applies when displaying a summary table:

```@example precision
show(IOContext(stdout, :digits => 2), MIME"text/plain"(), mean(chn))
```

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
