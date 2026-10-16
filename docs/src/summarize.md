# Summarize

## Extract summary values

`mean(chn)` and `summarystats(chn)` return a `ChainDataFrame` with parameter names in the `parameters` column:

```@example summary
using MCMCChains
chn = Chains(reshape(1:24, 4, 3, 2), [:a, :b, :c])
summary = mean(chn)
```

Index by parameter and statistic to get a value or vector:

```@example summary
summary[:a, :mean]
```

```@example summary
summary[[:a, :c], :mean]
```

Use `:` to keep all columns in a smaller table:

```@example summary
summary[[:a, :c], :]
```

Use Tables.jl to extract named columns, or convert to a DataFrame:

```@example summary
using Tables
Tables.columntable(summary)
```

```@example summary
using DataFrames
DataFrame(summary)
```

## API

The methods listed below are defined in `src/summarize.jl`.

```@autodocs
Modules = [MCMCChains]
Pages = ["summarize.jl"]
```
