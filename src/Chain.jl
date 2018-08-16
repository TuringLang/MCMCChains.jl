module Chain

export Chains

abstract type AbstractChains end

struct Chains <: AbstractChains
	value::Array{Float64, 3}
	range::AbstractRange{Int}
	names::Vector{AbstractString}
	chains::Vector{Int}
end

end # module
