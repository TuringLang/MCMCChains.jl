# Makie.jl plots

This page shows an example of plotting MCMCChains.jl with Makie.jl.
The example is meant to provide an useful basis to build upon.
Let's define some random chain and load the required packages:

```@example makie
using MCMCChains

chns = Chains(randn(300, 5, 3), [:A, :B, :C, :D, :E])
```

A basic way to visualize the chains is to show the drawn samples at each iteration.
Colors depict different chains.

```@example makie
using CairoMakie
CairoMakie.activate!(type = "svg")

params = names(chns, :parameters)

n_chains = length(chains(chns))
n_samples = length(chns)

fig = Figure(; resolution=(1_000, 800))

for (i, param) in enumerate(params)
    ax = Axis(fig[i, 1]; ylabel=string(param))
    for chain in 1:n_chains
        values = chns[:, param, chain]
        lines!(ax, 1:n_samples, values; label=string(chain))
    end

    hideydecorations!(ax; label=false)
    if i < length(params)
        hidexdecorations!(ax; grid=false)
    else
        ax.xlabel = "Iteration"
    end
end

fig
```

Next, we can add a second row of plots next to it which show the density estimate for these samples:

```@example makie
for (i, param) in enumerate(params)
    ax = Axis(fig[i, 2]; ylabel=string(param))
    for chain in 1:n_chains
        values = chns[:, param, chain]
        density!(ax, values; label=string(chain))
    end

    hideydecorations!(ax)
    if i < length(params)
        hidexdecorations!(ax; grid=false)
    else
        ax.xlabel = "Parameter estimate"
    end
end

axes = [only(contents(fig[i, 2])) for i in 1:length(params)]
linkxaxes!(axes...)

fig
```

Finally, let's add a simple legend.
Thanks to setting `label` above, this legend will have the right labels:

```@example makie
axislegend(first(axes))

fig
```
