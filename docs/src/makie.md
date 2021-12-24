# Makie.jl plots

Compared to Gadfly.jl and StatsPlots.jl, Makie.jl is the most flexible plotting library that you can use.
On this page is an example function for plotting with Makie.jl which you can use directly or further tweak to your needs.

```@example makie
using CairoMakie
using DataFrames
using MCMCChains

chns = Chains(randn(300, 5, 3), [:A, :B, :C, :D, :E])
```

```@example makie
function plot_chains(chns; density_func=density!)
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
    hideydecorations!.(values_axs; label=false)
    hidexdecorations!.(values_axs[1:end-1]; grid=false)

    # Create and store separate axes for showing the parameter density estimate.
    density_axs = [Axis(fig[i, 2]; ylabel=string(c)) for (i, c) in enumerate(params)]
    for (ax, col) in zip(density_axs, params)
        for i in 1:n_chains
            chain = string(i)
            values = filter(:chain => ==(chain), df)[:, col]
            density_func(ax, values; label=chain)
        end
    end

    # Just like above, we add some extra tweaks.
    density_axs[end].xlabel = "Parameter estimate"
    linkxaxes!(density_axs...)
    hideydecorations!.(density_axs)
    hidexdecorations!.(density_axs[1:end-1]; grid=false)

    return fig
end
nothing # hide
```

```@example makie
plot_chains(chns)
```
