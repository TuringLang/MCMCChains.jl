module MCMCChainsStatsPlotsExt

using MCMCChains
import StatsPlots

const _CORNER_MARGINS =
    (left_margin = 8 * StatsPlots.Plots.mm, bottom_margin = 3 * StatsPlots.Plots.mm)

_corner_kwargs(kwargs) = merge(_CORNER_MARGINS, (; kwargs...))

StatsPlots.cornerplot(c::Chains, args...; kwargs...) =
    MCMCChains.corner(c, args...; _corner_kwargs(kwargs)...)

StatsPlots.cornerplot!(c::Chains, args...; kwargs...) =
    MCMCChains.corner!(c, args...; _corner_kwargs(kwargs)...)

StatsPlots.cornerplot!(p::StatsPlots.Plots.Plot, c::Chains, args...; kwargs...) =
    MCMCChains.corner!(p, c, args...; _corner_kwargs(kwargs)...)

end
