module Chain

export Chains

abstract type AbstractChains end

struct Chains <: AbstractChains
	value::Array{Float64, 3}
	range::AbstractRange{Int}
	names::Vector{AbstractString}
	chains::Vector{Int}
end

# imports
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
#include("plot.jl")


end # module
