@recipe function f(chains::Chains, parameters::AbstractVector{Symbol}; colordim = :chain)
    colordim != :chain && error(
        "Symbol names are interpreted as parameter names, only compatible with ",
        "`colordim = :chain`",
    )

    ret = indexin(parameters, names(chains))
    any(y === nothing for y in ret) && error("Parameter not found")

    return chains, Int.(ret)
end

@recipe function f(
    chains::Chains,
    parameters::AbstractVector{<:Integer} = Int[];
    sections = _default_sections(chains),
    width = 500,
    height = 250,
    colordim = :chain,
    append_chains = false,
)
    _chains =
        isempty(parameters) ? Chains(chains, _clean_sections(chains, sections)) : chains
    c = append_chains ? pool_chain(_chains) : _chains
    ptypes = get(plotattributes, :seriestype, (:traceplot, :mixeddensity))
    ptypes = ptypes isa Symbol ? (ptypes,) : ptypes
    @assert all(ptype -> ptype ∈ supportedplots, ptypes)
    ntypes = length(ptypes)
    nrows, nvars, nchains = size(c)
    isempty(parameters) && (parameters = colordim == :chain ? (1:nvars) : (1:nchains))
    N = length(parameters)

    if :corner ∉ ptypes
        size --> (ntypes * width, N * height)
        # Read before the default below, so that a caller who asks for a legend keeps it on
        # every panel and a caller who asks for none gets none.
        asked_for_legend = haskey(plotattributes, :legend)
        legend --> false
        _apply_chrome!(plotattributes)

        multiple_plots = N * ntypes > 1
        if multiple_plots
            layout := (N, ntypes)
        end

        # Without this the panels are a set of unnamed coloured lines. One panel carries the
        # names by default, because repeating the same legend on every panel is ink that
        # says nothing new and lands on top of the data.
        legend_panel = asked_for_legend ? :none : plot_style().legend_panel

        i = 0
        for par in parameters
            for ptype in ptypes
                i += 1

                @series begin
                    if multiple_plots
                        subplot := i
                    end
                    if legend_panel === :all || (legend_panel === :first && i == 1)
                        legend := :best
                    end
                    colordim := colordim
                    seriestype := ptype
                    c, par
                end
            end
        end
    else
        ntypes > 1 && error(":corner is not compatible with multiple seriestypes")
        Corner(c, names(c)[parameters])
    end
end
