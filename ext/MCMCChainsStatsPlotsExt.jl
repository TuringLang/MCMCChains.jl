module MCMCChainsStatsPlotsExt

using MCMCChains
import StatsPlots

function _corner_input(c::Chains, parameters, sections)
    selected =
        if parameters === nothing || (parameters isa AbstractVector && isempty(parameters))
            Chains(c, MCMCChains._clean_sections(c, sections))
        else
            c[:, parameters isa Symbol ? [parameters] : parameters, :]
        end
    data = reduce(vcat, (Array(selected.value[:, :, i]) for i = 1:size(selected, 3)))
    return data, permutedims(string.(MCMCChains.names(selected)))
end

function _corner_kwargs(labels, kwargs)
    attrs = Dict{Symbol,Any}(kwargs)
    MCMCChains._apply_chrome!(attrs)
    get!(attrs, :label, labels)
    get!(attrs, :compact, length(labels) > 1)
    get!(attrs, :size, (600, 600))
    attrs[:left_margin] = get(kwargs, :left_margin, 8 * StatsPlots.Plots.mm)
    attrs[:bottom_margin] = get(kwargs, :bottom_margin, 3 * StatsPlots.Plots.mm)
    return attrs
end

function StatsPlots.cornerplot(
    c::Chains,
    parameters = nothing;
    sections = MCMCChains._default_sections(c),
    kwargs...,
)
    data, labels = _corner_input(c, parameters, sections)
    return StatsPlots.cornerplot(data; _corner_kwargs(labels, kwargs)...)
end

function StatsPlots.cornerplot!(
    c::Chains,
    parameters = nothing;
    sections = MCMCChains._default_sections(c),
    kwargs...,
)
    data, labels = _corner_input(c, parameters, sections)
    return StatsPlots.cornerplot!(data; _corner_kwargs(labels, kwargs)...)
end

function StatsPlots.cornerplot!(
    p::StatsPlots.Plots.Plot,
    c::Chains,
    parameters = nothing;
    sections = MCMCChains._default_sections(c),
    kwargs...,
)
    data, labels = _corner_input(c, parameters, sections)
    return StatsPlots.cornerplot!(p, data; _corner_kwargs(labels, kwargs)...)
end

end
