module MCMCChains

using AxisArrays
const axes = Base.axes

import AbstractMCMC
import AbstractMCMC: chainscat
# import DataFrames
# import DataFrames: eachcol, DataFrame
using Distributions
using RecipesBase
using SpecialFunctions
using Formatting
import StatsBase: autocov, counts, sem, AbstractWeights,
    autocor, describe, quantile, sample, summarystats, cov
using Requires

using LinearAlgebra: diag
import Serialization: serialize, deserialize
import Random
import Statistics: std, cor, mean, var

export Chains, chains, chainscat
export set_section, get_params, sections, sort_sections, setinfo, set_names
export autocor, describe, sample, summarystats, AbstractWeights, mean, quantile
export ChainDataFrame, DataFrame
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
struct Chains{A, T, K<:NamedTuple, L<:NamedTuple} <: AbstractMCMC.AbstractChains
    value::AxisArray{A,3}
    logevidence::T
    name_map::K
    info::L
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
# include("modelstats.jl")
include("plot.jl")

function __init__()
    @require DataFrames="a93c6f00-e57d-5684-b7b6-d8193f3e46c0" include("dataframes-compat.jl")
end

end # module
