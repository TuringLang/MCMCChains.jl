module MCMCChains

using AxisArrays
const axes = Base.axes

using AbstractMCMC
import DataFrames: DataFrame
using Distributions
using RecipesBase
using SpecialFunctions
using StatsBase: autocov, counts, sem, AbstractWeights
import StatsBase: autocor, describe, quantile, sample, summarystats

using LinearAlgebra: diag
import Serialization: serialize, deserialize
import Random
import Statistics: std, cor, mean

export Chains, chains, chainscat
export set_section, get_params, sections, sort_sections, setinfo, set_names
export mean
export autocor, describe, sample, summarystats, AbstractWeights
export ChainDataFrame, DataFrame
export summarize

# export diagnostics functions
export discretediag, gelmandiag, gewekediag, heideldiag, rafterydiag
export hpd, ess

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
include("constructors.jl")
include("summarize.jl")
include("discretediag.jl")
include("fileio.jl")
include("gelmandiag.jl")
include("gewekediag.jl")
include("heideldiag.jl")
include("mcse.jl")
include("rafterydiag.jl")
include("sampling.jl")
include("stats.jl")
include("modelstats.jl")
include("plot.jl")

end # module
