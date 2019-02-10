module MCMCChain

import Showoff: showoff
import StatsBase: autocor, autocov, countmap, counts, describe, predict,
       quantile, sample, sem, summarystats
import LinearAlgebra: diag
import Serialization: serialize, deserialize

using RecipesBase
import RecipesBase: plot

using Distributions
using SpecialFunctions

export Chains, getindex, setindex!
export plot, traceplot, meanplot, densityplot, histogramplot, mixeddensityplot, autcorplot
export describe

# export diagnostics functions
export discretediag, gelmandiag, gewekediag, heideldiag, rafterydiag

abstract type AbstractChains end

"""
    Chains type

Parameters:

- `value`: `iterations × variables × chains` Data array
- `range`: Range describing the iterations (considering thinning)
- `names`: List of variable names (strings)
- `chains`: List of chain ids
"""
struct Chains{T<:Real} <: AbstractChains
    logevidence::T
    value::Array{Union{Missing, T}, 3}
    range::AbstractRange{Int}
    names::Vector
    uniquenames::Dict{Symbol, Int}
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
#include("plot2.jl")

end # module
