module MCMCChains

using AxisArrays
const axes = Base.axes

import AbstractMCMC
import AbstractMCMC: chainscat
using Distributions
using RecipesBase
using SpecialFunctions
using Formatting
import StatsBase: autocov, counts, sem, AbstractWeights,
    autocor, describe, quantile, sample, summarystats, cov
using Requires
import PrettyTables
import Tables
import TableTraits
import IteratorInterfaceExtensions

using LinearAlgebra: diag
import Serialization: serialize, deserialize
import Random
import Statistics: std, cor, mean, var

export Chains, chains, chainscat
export set_section, get_params, sections, sort_sections, setinfo
export replacenames, namesingroup, group
export autocor, describe, sample, summarystats, AbstractWeights, mean, quantile
export ChainDataFrame
export summarize

# Export diagnostics functions
export discretediag, gelmandiag, gewekediag, heideldiag, rafterydiag
export hpd, ess

"""
    Chains

Parameters:

- `value`: An `AxisArray` object with axes `iter` × `var` × `chains`
- `logevidence` : A field containing the logevidence.
- `name_map` : A `NamedTuple` mapping each variable to a section.
- `info` : A `NamedTuple` containing miscellaneous information relevant to the chain.
The `info` field can be set using `setinfo(c::Chains, n::NamedTuple)`.
"""
struct Chains{T,A<:AxisArray{T,3},L,K<:NamedTuple,I<:NamedTuple} <: AbstractMCMC.AbstractChains
    value::A
    logevidence::L
    name_map::K
    info::I
end

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
include("tables.jl")

end # module
