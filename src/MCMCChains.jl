module MCMCChains

import Showoff: showoff
import StatsBase: autocor, autocov, countmap, counts, describe, predict,
       quantile, sample, sem, summarystats, sample, AbstractWeights
import LinearAlgebra: diag
import Serialization: serialize, deserialize
import Base: sort, range, names, get, hash, convert, show, display
import Statistics: cor
import Core.Array
import DataFrames: DataFrame, names

using RecipesBase
import RecipesBase: plot

using Serialization
using Distributions
using KernelDensity
using SpecialFunctions
using AxisArrays
const axes = Base.axes

export Chains, getindex, setindex!, chains, setinfo, chainscat
export describe, set_section, get_params, sections
export sample, AbstractWeights
export Array, DataFrame, sort_sections, convert
export summarize, summarystats, ChainDataFrame

# export diagnostics functions
export discretediag, gelmandiag, gewekediag, heideldiag, rafterydiag
export autocor

abstract type AbstractChains end

"""
    Chains type

Parameters:

- `value`: An `AxisArray` object with axes `iter` × `var` × `chains`
- `logevidence` : A field containing the logevidence.
- `name_map` : A `NamedTuple` mapping each variable to a section.
- `info` : A `NamedTuple` containing miscellaneous information relevant to the chain.
The `info` field can be set using `setinfo(c::Chains, n::NamedTuple)`.
"""
struct Chains{A, T, K<:NamedTuple, L<:NamedTuple} <: AbstractChains
    value::AxisArray{A,3}
    logevidence::T
    name_map::K
    info::L
end

# imports
include("utils.jl")

include("chains.jl")
include("chainsummary.jl")
include("constructors.jl")
include("summarize.jl")
include("discretediag.jl")
include("fileio.jl")
include("gelmandiag.jl")
include("gewekediag.jl")
include("heideldiag.jl")
include("mcse.jl")
#include("modelchains.jl")
#include("modelstats.jl")
include("rafterydiag.jl")
include("sampling.jl")
include("stats.jl")
include("plot.jl")
#include("plot2.jl")

end # module
