@recipe function f(p::NutsPlot; kind = :acceptance)
    chn = _chain_arg(p, "nutsplot")
    _apply_chrome!(plotattributes)

    val, label = _nuts_statistic(chn, kind)
    chain_ids = MCMCChains.chains(chn)

    if kind === :divergence
        divergences = vec(sum(val .> 0; dims = 1))
        title --> "Divergent transitions"
        xticks --> (1:length(chain_ids), string.(chain_ids))
        xaxis --> "Chain"
        yaxis --> "Draws"
        legend --> false
        seriestype := :bar
        bar_width --> 0.6
        linewidth --> 0
        fillcolor --> _chain_colours(length(chain_ids))
        collect(1:length(chain_ids)), divergences
    else
        title --> label
        xaxis --> label
        yaxis --> (kind === :treedepth ? "Draws" : "Density")
        label --> permutedims(["Chain $c" for c in chain_ids])
        seriestype := (kind === :treedepth ? :histogram : :density)
        fillalpha --> FILL_ALPHA
        val
    end
end
