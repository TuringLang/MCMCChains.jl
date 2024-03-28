module MCMCChainsMakieExt

using MCMCChains
import Makie


@recipe(MCMCChains.Chains) do scene
    Theme()
end

function Makie.plot!(chns::MCMCChains.Chains)
    params = names(chns, :parameters)

    n_chains = length(chains(chns))
    n_samples = length(chns)

    fig = Figure()

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

    axislegend(first(axes))

    return fig
end

end