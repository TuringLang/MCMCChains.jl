module MCMCChainsMakieExt

import MCMCChains
import Makie


function MCMCChains.trace(chns::T) where {T<:MCMCChains.Chains}
    params = MCMCChains.names(chns, :parameters)

    n_chains = length(MCMCChains.chains(chns))
    n_samples = length(chns)
    n_params = length(params)

    colors = Makie.to_colormap(:tol_vibrant)
    width = 600
    height = max(400, 80 * n_params)

    fig = Makie.Figure(; size=(width, height))

    for (i, param) in enumerate(params)
        ax = Makie.Axis(fig[i+1, 1]; ylabel=string(param))
        for chain in 1:n_chains
            values = chns[:, param, chain]
            Makie.lines!(
                ax,
                1:n_samples,
                values;
                label=string(chain),
                color=(colors[chain], 0.7),
                linewidth=0.7
            )
        end

        Makie.hideydecorations!(ax; label=false)
        if i < n_params
            Makie.hidexdecorations!(ax; grid=false)
        else
            ax.xlabel = "Iteration"
        end
    end

    for (i, param) in enumerate(params)
        ax = Makie.Axis(fig[i+1, 2]; ylabel=string(param))
        for chain in 1:n_chains
            values = chns[:, param, chain]
            Makie.density!(
                ax,
                values;
                label=string(chain),
                color=(colors[chain], 0.1),
                strokewidth=1,
                strokecolor=(colors[chain], 0.7)
            )
        end

        Makie.hideydecorations!(ax)
        if i < n_params
            Makie.hidexdecorations!(ax; grid=false)
        else
            ax.xlabel = "Parameter estimate"
        end
    end

    axes = [only(Makie.contents(fig[i+1, 2])) for i in 1:n_params]
    Makie.linkxaxes!(axes...)

    Makie.Legend(fig[1, 1:2], first(axes), "Chain", orientation=:horizontal, titlehalign=:left, halign=:left, titleposition=:left)

    Makie.rowgap!(fig.layout, 10)
    Makie.colgap!(fig.layout, 10)

    return fig
end
end