# Makie.jl plots

This page shows an example of plotting MCMCChains.jl with Makie.jl.
The example is meant to provide an useful basis to build upon.
Let's define some random chain and load the required packages:

```@example makie
using CairoMakie
using DataFrames
using MCMCChains

chns = Chains(randn(300, 5, 3), [:A, :B, :C, :D, :E])
```

A basic way to visualize the chains is to show the drawn samples at each iteration.
Colors depict different chains:

```@example makie
params = names(chns, :parameters)

df = DataFrame(chns)
n_chains = length(unique(df.chain))
n_samples = nrow(df) / n_chains

# Alternatively, use `CategoricalArrays.categorical`.
df[!, :chain] = string.(df.chain)

fig = Figure(; resolution=(1_000, 800))

# Create and store separate axes for showing iterations.
values_axs = [Axis(fig[i, 1]; ylabel=string(c)) for (i, c) in enumerate(params)]
for (ax, col) in zip(values_axs, params)
    for i in 1:n_chains
        chain = string(i)
        values = filter(:chain => ==(chain), df)[:, col]
        lines!(ax, 1:n_samples, values; label=chain)
    end
end

# Thanks to having stored the axes before, we can apply some extra tweaks.
# These tweaks usually depend on the kind of model being fitted.
values_axs[end].xlabel = "Iteration"
hidexdecorations!.(values_axs[1:end-1]; grid=false)

fig
```

Next, we can add a second row of plots next to it which show the density estimate for these samples:

```@example makie
density_axs = [Axis(fig[i, 2]; ylabel=string(c)) for (i, c) in enumerate(params)]
for (ax, col) in zip(density_axs, params)
    for i in 1:n_chains
        chain = string(i)
        values = filter(:chain => ==(chain), df)[:, col]
        density!(ax, values; label=chain)
    end
end

# Just like above, we add some extra tweaks.
density_axs[end].xlabel = "Parameter estimate"
linkxaxes!(density_axs...)
hideydecorations!.(density_axs)
hidexdecorations!.(density_axs[1:end-1]; grid=false)

# And we can also drop the y ticks for the plots on the left.
hideydecorations!.(values_axs; label=false)

fig
```

Finally, let's add a simple legend.
Thanks to setting `label` above, this legend will have the right entries:

```@example makie
axislegend(first(density_axs))

fig
```
