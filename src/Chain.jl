module Chain

import Showoff: showoff
import StatsBase: autocor, autocov, countmap, counts, describe, predict,
       quantile, sample, sem, summarystats
# import Gadfly: draw, Geom, Guide, Layer, layer, PDF, PGF, Plot, plot, PNG, PS,
#       render, Scale, SVG, Theme
import Plots: plot, histogram, contour

export Chains, plot

#using Plots: Plot, unicodeplots, pyplot, gr, histogram, contour

abstract type AbstractChains end

struct Chains <: AbstractChains
	value::Array{Float64, 3}
	range::AbstractRange{Int}
	names::Vector{AbstractString}
	chains::Vector{Int}
end

# imports
include("utils.jl")

include("chains.jl")
include("chainsummary.jl")
include("discretediag.jl")
include("fileio.jl")
include("gelmandiag.jl")
include("gewekediag.jl")
include("heideldiag.jl")
include("mcse.jl")
#include("modelchains.jl")
#include("modelstats.jl")
include("rafterydiag.jl")
include("stats.jl")
include("plot.jl")

end # module
