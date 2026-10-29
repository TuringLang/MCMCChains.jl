# Summarize

## Extract summary values

`mean(chn)` and `summarystats(chn)` return a `ChainDataFrame` with parameter names in the `parameters` column:

```@example summary
using MCMCChains
chn = Chains(reshape(1:120, 20, 3, 2), [:a, :b, :c])
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
Tables.columntable(summary).mean
```

```@example summary
using DataFrames
DataFrame(summary)
```

`Array(summary)` returns a vector for one statistic or a matrix for several statistics. Matrix rows follow `summary[:, :parameters]`, and columns follow `names(summary)[2:end]`.

`convert(Array, summary)` keeps the two-dimensional table layout even for one statistic.

The vector follows the parameter order in the chain. Sort the grouped chain first if you need natural index order:

```@example summary
grouped_chn = Chains(reshape(1:120, 20, 3, 2), ["x[2]", "x[1]", "y"])
grouped_mean = mean(sort(group(grouped_chn, :x)))
Array(grouped_mean)
```

The conversion interface preserves the statistic dimension:

```@example summary
size(convert(Array, grouped_mean))
```

## API

The methods listed below are defined in `src/summarize.jl`.

```@autodocs
Modules = [MCMCChains]
Pages = ["summarize.jl"]
```
