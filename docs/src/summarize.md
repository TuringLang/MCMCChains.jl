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

`Array(summary)` returns the statistic columns as a matrix. Rows follow `summary[:, :parameters]`, and columns follow `names(summary)[2:end]`.

For a vector of means from a grouped parameter, select the named column:

```@example summary
grouped_chn = Chains(reshape(1:120, 20, 3, 2), ["x[1]", "x[2]", "y"])
grouped_mean = mean(group(grouped_chn, :x))
Tables.getcolumn(grouped_mean, :mean)
```

The matrix form keeps its statistic dimension, including when it contains only one column:

```@example summary
size(Array(grouped_mean))
```

## API

The methods listed below are defined in `src/summarize.jl`.

```@autodocs
Modules = [MCMCChains]
Pages = ["summarize.jl"]
```
