module MCMCChainsMakieExt

import MCMCChains
import Makie


# @recipe(MCMCChains.Chains) do scene
#     Theme()
# end

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


function MCMCChains.myplot(chns::T) where {T<:MCMCChains.Chains}
    params = names(chns, :parameters)

    n_chains = length(MCMCChains.chains(chns))
    n_samples = length(chns)

    colors = Makie.cgrad(:roma100, n_chains, categorical=true)
    fig = Makie.Figure()

    for (i, param) in enumerate(params)
        ax = Makie.Axis(fig[i+1, 1]; ylabel=string(param))
        for chain in 1:n_chains
            values = chns[:, param, chain]
            Makie.lines!(ax, 1:n_samples, values; label=string(chain), color=(colors[chain], 0.7), linewidth=0.7)
        end

        Makie.hideydecorations!(ax; label=false)
        if i < length(params)
            Makie.hidexdecorations!(ax; grid=false)
        else
            ax.xlabel = "Iteration"
        end
    end

    for (i, param) in enumerate(params)
        ax = Makie.Axis(fig[i+1, 2]; ylabel=string(param))
        for chain in 1:n_chains
            values = chns[:, param, chain]
            Makie.density!(ax, values; label=string(chain), color=(colors[chain], 0.7))
        end

        Makie.hideydecorations!(ax)
        if i < length(params)
            Makie.hidexdecorations!(ax; grid=false)
        else
            ax.xlabel = "Parameter estimate"
        end
    end

    axes = [only(Makie.contents(fig[i+1, 2])) for i in 1:length(params)]
    Makie.linkxaxes!(axes...)

    Makie.Legend(fig[1, 1:2], first(axes), "Chain", orientation=:horizontal)

    return fig
end
end