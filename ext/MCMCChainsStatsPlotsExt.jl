module MCMCChainsStatsPlotsExt

using MCMCChains
import StatsPlots

StatsPlots.cornerplot(c::Chains, args...; kwargs...) =
    MCMCChains.corner(c, args...; kwargs...)

StatsPlots.cornerplot!(c::Chains, args...; kwargs...) =
    MCMCChains.corner!(c, args...; kwargs...)

StatsPlots.cornerplot!(p::StatsPlots.Plots.Plot, c::Chains, args...; kwargs...) =
    MCMCChains.corner!(p, c, args...; kwargs...)

end
