module MCMCChains

import Showoff: showoff
import StatsBase: autocor, autocov, countmap, counts, describe, predict,
       quantile, sample, sem, summarystats
import LinearAlgebra: diag
import Serialization: serialize, deserialize
import Base: sort, range, names
import Statistics: cor

using RecipesBase
import RecipesBase: plot

using Serialization
using Distributions
using SpecialFunctions
using AxisArrays
const axes = Base.axes

export Chains, getindex, setindex!, chains
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
struct Chains{A, T} <: AbstractChains
    value::AxisArray{Union{Missing,A},3}
    logevidence::T
    name_map::Dict{Any, Vector}
    info::Dict{Symbol, Any}
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
include("rafterydiag.jl")
include("stats.jl")
include("plot.jl")

end # module
