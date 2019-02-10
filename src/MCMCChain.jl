module MCMCChain

import Showoff: showoff
import StatsBase: autocor, autocov, countmap, counts, describe, predict,
       quantile, sample, sem, summarystats
import LinearAlgebra: diag
import Base: sort

using RecipesBase
import RecipesBase: plot

using Serialization
using Distributions
using SpecialFunctions
using AxisArrays

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
    value::AxisArray
    logevidence::T
    name_map::Dict{Any, Vector}
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
