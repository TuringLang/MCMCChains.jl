# Gadfly.jl plots

To plot the Chains via [Gadfly.jl](https://github.com/GiovineItalia/Gadfly.jl), use the DataFrames constructor:

```@example gadfly
using DataFrames
using CategoricalArrays
using Gadfly
write_svg(path, p; w=6inch, h=4inch) = Gadfly.draw(Gadfly.SVG(path, w, h), p) # hide
using MCMCChains

chn = Chains(randn(100, 2, 3), [:A, :B])
df = DataFrame(chn)
df[!, :chain] = categorical(df.chain)

plot(df, x=:A, color=:chain, Geom.density, Guide.ylabel("Density"))
```

## Multiple parameters

Or, to show multiple parameters in one plot, use `DataFrames.stack`

```@example gadfly
sdf = stack(df, names(chn), variable_name=:parameter)
first(sdf, 5)
```

and `Gadfly.Geom.subplot_grid`

```@example gadfly
plot(sdf, ygroup=:parameter, x=:value, color=:chain,
    Geom.subplot_grid(Geom.density), Guide.ylabel("Density"))
```

This is very flexible.
For example, we can look at the first two chains only by using `DataFrames.filter`

```@example gadfly
first_chain = filter([:chain] => c -> c == 1 || c == 2, sdf)

plot(first_chain, xgroup=:parameter, ygroup=:chain, x=:value,
    Geom.subplot_grid(Geom.density, Guide.xlabel(orientation=:horizontal)),
    Guide.xlabel("Parameter"), Guide.ylabel("Chain"))
```

## Trace

```@example gadfly
plot(first_chain, ygroup=:parameter, x=:iteration, y=:value, color=:chain,
    Geom.subplot_grid(Geom.point), Guide.ylabel("Sample value"))
```
