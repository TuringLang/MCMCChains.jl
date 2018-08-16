module Chain


abstract type AbstractChains end

struct Chains <: AbstractChains
	value::Array{Float64, 3}
  range::Range{Int}
  names::Vector{AbstractString}
  chains::Vector{Int}
end

struct ModelChains <: AbstractChains
	value::Array{Float64, 3}
  range::Range{Int}
  names::Vector{AbstractString}
  chains::Vector{Int}
  #model::Model
end

end # module
