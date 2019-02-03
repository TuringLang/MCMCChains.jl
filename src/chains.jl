#################### Chains ####################

## Constructors ##

# Default name map.
const DEFAULT_MAP = Dict{Symbol, Vector{Symbol}}(:parameters => [])

# Set default parameter names if not given.
function Chains(val::AbstractArray{T,3}) where T <: Real
    parameter_names = map(i->Symbol("Param$i"), 1:size(val, 2))
    return Chains(val, parameter_names)
end

# Generic chain constructor.
function Chains(val::AbstractArray{T,3},
    parameter_names::Vector,
    name_map = DEFAULT_MAP) where T <: Real

    # If we received an array of pairs, convert it to a dictionary.
    if typeof(name_map) <: Array
        name_map = Dict(name_map)
    end

    names = [:iter, :var, :chain]
    axvals = [
              1:size(val, 1),
              parameter_names,
              map(i->Symbol("Chain$i"), 1:size(val, 3)),
    ]

    if length(keys(name_map)) == 1
        name_map[first(keys(name_map))] = parameter_names
    else
        # Check that all parameters are assigned.
        for param in parameter_names
            found = false
            for (_, assigned) in name_map
                if param in assigned
                    found = true
                    break
                end
            end

            # Assign to :parameters by default.
            if !found
                push!(name_map[:parameters], param)
            end
        end
    end

    axs = ntuple(i -> Axis{names[i]}(axvals[i]), 3)
    A = AxisArray(convert(Array{Union{Missing,T},3}, val), axs...)
    return Chains{T}(A, zero(T), name_map)
end


#################### Indexing ####################

# If only a single argument is given, index by parameter.
Base.getindex(c::Chains, i) = getindex(c.value, :, i, :)
Base.setindex!(c::Chains, v, i) = setindex!(c.value, v, :, i, :)

# Otherwise, forward the indexing to AxisArrays.
Base.getindex(c::Chains, i...) = getindex(c.value, i...)
Base.setindex!(c::Chains, v, i...) = setindex!(c.value, v, i...)

#################### Concatenation ####################

function Base.cat(dim::Integer, c1::AbstractChains, args::AbstractChains...)
  dim == 1 ? cat1(c1, args...) :
  dim == 2 ? cat2(c1, args...) :
  dim == 3 ? cat3(c1, args...) :
    throw(ArgumentError("cannot concatenate along dimension $dim"))
end

function cat1(c1::AbstractChains, args::AbstractChains...)
  range = c1.range
  for c in args
    last(range) + step(range) == first(c) ||
      throw(ArgumentError("noncontiguous chain iterations"))
    step(range) == step(c) ||
      throw(ArgumentError("chain thinning differs"))
    range = first(range):step(range):last(c)
  end

  names = c1.names
  all(c -> c.names == names, args) ||
    throw(ArgumentError("chain names differ"))

  chains = c1.chains
  all(c -> c.chains == chains, args) ||
    throw(ArgumentError("sets of chains differ"))

  value = cat(1, c1.value, map(c -> c.value, args)...)
  Chains(value, start=first(range), thin=step(range), names=names,
         chains=chains)
end

function cat2(c1::AbstractChains, args::AbstractChains...)
  range = c1.range
  all(c -> c.range == range, args) ||
    throw(ArgumentError("chain ranges differ"))

  names = c1.names
  n = length(names)
  for c in args
    names = union(names, c.names)
    n += length(c.names)
    n == length(names) ||
      throw(ArgumentError("non-unique chain names"))
  end

  chains = c1.chains
  all(c -> c.chains == chains, args) ||
    throw(ArgumentError("sets of chains differ"))

  value = cat(2, c1.value, map(c -> c.value, args)...)
  Chains(value, start=first(range), thin=step(range), names=names,
         chains=chains)
end

function cat3(c1::AbstractChains, args::AbstractChains...)
  range = c1.range
  all(c -> c.range == range, args) ||
    throw(ArgumentError("chain ranges differ"))

  names = c1.names
  all(c -> c.names == names, args) ||
    throw(ArgumentError("chain names differ"))

  value = cat(3, c1.value, map(c -> c.value, args)...)
  Chains(value, start=first(range), thin=step(range), names=names)
end

Base.hcat(c1::AbstractChains, args::AbstractChains...) = cat(2, c1, args...)

Base.vcat(c1::AbstractChains, args::AbstractChains...) = cat(1, c1, args...)


#################### Base Methods ####################

function Base.keys(c::AbstractChains)
    names(c)
end

function Base.show(io::IO, c::Chains)
    print(io, "Object of type \"$(summary(c))\"\n\n")
    println(io, header(c))
    show(io, summarystats(c))
end

function Base.size(c::AbstractChains)
    dim = size(c.value)
    last(c), dim[2], dim[3]
end

function Base.size(c::AbstractChains, ind)
    size(c)[ind]
end

Base.first(c::AbstractChains) = first(c.value[Axis{:iter}].val)
Base.step(c::AbstractChains) = step(c.value[Axis{:iter}].val)
Base.last(c::AbstractChains) = last(c.value[Axis{:iter}].val)

#################### Auxilliary Functions ####################

function combine(c::AbstractChains)
  n, p, m = size(c.value)
  value = Array{Float64}(n * m, p)
  for j in 1:p
    idx = 1
    for i in 1:n, k in 1:m
      value[idx, j] = c.value[i, j, k]
      idx += 1
    end
  end
  value
end

function range(c::AbstractChains)
    return c.value[Axis{:iter}].val
end

function chains(c::AbstractChains)
    return c.value[Axis{:chain}].val
end

function names(c::AbstractChains)
    return c.value[Axis{:var}].val
end

function header(c::AbstractChains)
    rng = range(c)
    string(
        # "Log model evidence = $(c.logevidence)\n", #FIXME: Uncomment.
        "Iterations        = $(first(c)):$(last(c))\n",
        "Thinning interval = $(step(c))\n",
        "Chains            = $(join(map(string, chains(c)), ", "))\n",
        "Samples per chain = $(length(range(c)))\n",
        "Parameters        = $(join(map(string, names(c)), ", "))\n"
    )
end

function indiscretesupport(c::AbstractChains,
                           bounds::Tuple{Real, Real}=(0, Inf))
  nrows, nvars, nchains = size(c.value)
  result = Array{Bool}(undef, nvars * (nrows > 0))
  for i in 1:nvars
    result[i] = true
    for j in 1:nrows, k in 1:nchains
      x = c.value[j, i, k]
      if !isinteger(x) || x < bounds[1] || x > bounds[2]
        result[i] = false
        break
      end
    end
  end
  result
end

function link(c::AbstractChains)
  cc = copy(c.value)
  for j in 1:length(c.names)
    x = cc[:, j, :]
    if minimum(x) > 0.0
      cc[:, j, :] = maximum(x) < 1.0 ? logit.(x) : log.(x)
    end
  end
  cc
end
