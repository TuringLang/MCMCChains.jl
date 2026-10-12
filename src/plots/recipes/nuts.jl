@recipe function f(p::NutsPlot; kind = :acceptance)
    chn = _chain_arg(p, "nutsplot")
    _apply_chrome!(plotattributes)

    val, label = _nuts_statistic(chn, kind)
    chain_ids = MCMCChains.chains(chn)

    if kind === :divergence
        divergences = vec(sum(val .> 0; dims = 1))
        colours = _chain_colours(length(chain_ids))
        title --> "Divergent transitions"
        xticks --> (1:length(chain_ids), string.(chain_ids))
        xaxis --> "Chain"
        yaxis --> "Draws"
        legend --> false
        # One series per bar. Plots hands a row vector of colours out one series at a time,
        # not one bar at a time, so a single series would paint every chain alike. The last
        # bar is returned rather than emitted, because a recipe that ends without data
        # still adds a series.
        for i = 1:(length(chain_ids)-1)
            @series begin
                seriestype := :bar
                bar_width --> 0.6
                linewidth --> 0
                label := nothing
                fillcolor := colours[i]
                [i], [divergences[i]]
            end
        end
        seriestype := :bar
        bar_width --> 0.6
        linewidth --> 0
        label := nothing
        fillcolor := colours[end]
        [length(chain_ids)], [divergences[end]]
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
