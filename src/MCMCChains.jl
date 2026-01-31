module MCMCChains

using AxisArrays
const axes = Base.axes
import AbstractMCMC
import AbstractMCMC: chainscat, chainsstack
using Distributions
using RecipesBase
using Dates
using KernelDensity: kde, pdf
import DataAPI
import StatsBase:
    autocov,
    counts,
    sem,
    AbstractWeights,
    autocor,
    describe,
    quantile,
    sample,
    summarystats,
    cov,
    ecdf

import MCMCDiagnosticTools
import MLJModelInterface
import NaturalSort
import OrderedCollections
import PrettyTables
import StatsFuns
import Tables
import TableTraits
import IteratorInterfaceExtensions

import LinearAlgebra
import Random
import Statistics: std, cor, mean, var, mean!

export Chains, chains, chainscat, chainsstack
export setrange, resetrange
export set_section, get_params, sections, sort_sections, setinfo
export replacenames, namesingroup, group
export autocor, describe, sample, summarystats, AbstractWeights, mean, quantile
export ChainDataFrame
export summarize
export write_stancsv, read_stancsv

# Reexport diagnostics functions
using MCMCDiagnosticTools:
    discretediag,
    ess,
    ess_rhat,
    AutocovMethod,
    FFTAutocovMethod,
    BDAAutocovMethod,
    gelmandiag,
    gelmandiag_multivariate,
    gewekediag,
    heideldiag,
    mcse,
    rafterydiag,
    rhat,
    rstar
export discretediag
export ess, ess_rhat, rhat, AutocovMethod, FFTAutocovMethod, BDAAutocovMethod
export gelmandiag, gelmandiag_multivariate
export gewekediag
export heideldiag
export mcse
export rafterydiag
export rstar

export hpd

"""
    Chains

Parameters:

- `value`: An `AxisArray` object with axes `iter` × `var` × `chains`
- `logevidence` : A field containing the logevidence.
- `name_map` : A `NamedTuple` mapping each variable to a section.
- `info` : A `NamedTuple` containing miscellaneous information relevant to the chain.
The `info` field can be set using `setinfo(c::Chains, n::NamedTuple)`.
"""
struct Chains{T,A<:AxisArray{T,3},L,K<:NamedTuple,I<:NamedTuple} <:
       AbstractMCMC.AbstractChains
    value::A
    logevidence::L
    name_map::K
    info::I
end

include("utils.jl")
include("chains.jl")
include("constructors.jl")
include("ess_rhat.jl")
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
include("rstar.jl")

function write_stancsv end
function read_stancsv end

const _CSV_REQUIRED_MSG = """
StanCSV functions require the CSV.jl package.
Please load it first: `using CSV`
"""

write_stancsv(file, chn::Chains; kwargs...) = error(_CSV_REQUIRED_MSG)
write_stancsv(file, chn::Chains, ::Val{:all}; kwargs...) = error(_CSV_REQUIRED_MSG)
read_stancsv(file) = error(_CSV_REQUIRED_MSG)
read_stancsv(files::AbstractVector) = error(_CSV_REQUIRED_MSG)

end # module
