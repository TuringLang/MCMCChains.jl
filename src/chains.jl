#################### Chains ####################

## Constructors ##

# Default name map.
const DEFAULT_MAP = Dict{Symbol, Vector{Any}}(:parameters => [])

# Set default parameter names if not given.
function Chains(val::AbstractArray{A,3};
        start::Int=1,
        thin::Int=1,
        evidence = 0.0) where {A<:Union{Real, Union{Missing, Real}}}
    parameter_names = map(i->"Param$i", 1:size(val, 2))
    return Chains(val, parameter_names, start=start, thin=thin)
end

# Generic chain constructor.
function Chains(val::AbstractArray{A,3},
        parameter_names::Vector,
        name_map = copy(DEFAULT_MAP);
        start::Int=1,
        thin::Int=1,
        evidence = 0.0) where {A<:Union{Real, Union{Missing, Real}}}

    # If we received an array of pairs, convert it to a dictionary.
    if typeof(name_map) <: Array
        name_map = Dict(name_map)
    end

    # Make sure that we have a :parameters index.
    if !in(:parameters, keys(name_map))
        name_map[:parameters] = []
    end

    # Preclean the name_map of names that aren't in the
    # parameter_names vector.
    for section in keys(name_map)
        filter!(x -> x ∈ parameter_names, name_map[section])
    end

    # Provide an alias for start, to avoid clashing with the range call.
    names = [:iter, :var, :chain]
    axvals = [
        Base.range(start, step=thin, length=size(val, 1)),
        parameter_names,
        map(i->Symbol("Chain$i"), 1:size(val, 3)),
    ]

    # Check if we should stick anything in the :internals field.
    removed_names = []
    for param in parameter_names
        if endswith(string(param), "_") || startswith(string(param), "_")
            # Add the internals key if it doesn't already exist.
            in(:internals, keys(name_map)) || (name_map[:internals] = [])
            push!(name_map[:internals], param)
        end
    end

    # Remove the names we put in the internals field.
    filter!(x -> x ∉ removed_names, parameter_names)

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

            # Assign to :parameters by default, or :internals if it starts
            # or ends with an underscore.
            if !found
                if endswith(string(param), "_") || startswith(string(param), "_")
                    # Add the internals key if it doesn't already exist.
                    in(:internals, keys(name_map)) || (name_map[:internals] = [])
                    push!(name_map[:internals], param)
                else
                    push!(name_map[:parameters], param)
                end
            end
        end
    end

    axs = ntuple(i -> Axis{names[i]}(axvals[i]), 3)
    arr = AxisArray(convert(Array{Union{Missing,A},3}, val), axs...)
    return sort(Chains{typeof(evidence)}(arr, evidence, name_map, Dict{Symbol, Any}()))
end

# Retrieve a new chain with only a specific section pulled out.
function Chains(c::Chains{T}, section::Union{Vector, Any};
                sorted=false) where {T<:Real}
    section = typeof(section) <: AbstractArray ? section : [section]

    # If we received an empty list, return the original chain.
    if isempty(section) return c end

    # Gather the names from the relevant sections.
    names = []
    for s in section
        # Make sure the section exists first.
        in(s, keys(c.name_map)) ||
            throw(ArgumentError("$section not found in Chains name map."))

        names = vcat(names, c.name_map[s])
    end

    # Extract wanted values.
    new_vals = c.value[:, names, :]

    # Create the new chain.
    new_chn = Chains{T}(new_vals,
        c.logevidence,
        Dict([s => c.name_map[s] for s in section]),
        c.info)

    if sorted
        return sort(new_chn)
    else
        return new_chn
    end
end


#################### Indexing ####################

Base.getindex(c::Chains, i1::T) where T<:Union{AbstractUnitRange, StepRange} = c[i1, :, :]
Base.getindex(c::Chains, i1::Integer) = c[i1:i1, :, :]
Base.getindex(c::Chains, v::Symbol) = c[[v]]
Base.getindex(c::Chains, v::String) = c[:, [v], :]
Base.getindex(c::Chains, v::Vector{String}) = c.value[:, v, :]

function Base.getindex(c::Chains, v::Vector{Symbol})
    v_str = string.(v)
    idx = indexin(v_str, names(c))
    syms = []
    for i in eachindex(idx)
        value = v_str[i]
        if idx[i] == nothing
            append!(syms,
                collect(Iterators.filter(k -> occursin(value*"[", string(k)), names(c))))
        else
            push!(syms, value)
        end
    end

    sort!(syms, lt=MCMCChain.natural)
    return c.value[:, syms, :]
end

function Base.getindex(c::Chains{T}, i...) where T
    # Make sure things are in array form to preserve the axes.
    ind = (i[1],
           typeof(i[2]) <: Union{AbstractArray, Colon} ?  i[2] : [i[2]],
           typeof(i[3]) <: Union{AbstractArray, Colon} ?  i[3] : [i[3]]
    )

    newval = getindex(c.value, ind...)
    names = newval.axes[2].val
    return Chains{T}(newval,
        c.logevidence,
        _trim_name_map(names, c.name_map),
        c.info)
end

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
    # show(io, summarystats(c))
end

function Base.size(c::AbstractChains)
    dim = size(c.value)
    last(c), dim[2], dim[3]
end

function Base.size(c::AbstractChains, ind)
    size(c)[ind]
end

Base.length(c::AbstractChains) = length(range(c))
Base.first(c::AbstractChains) = first(c.value[Axis{:iter}].val)
Base.step(c::AbstractChains) = step(c.value[Axis{:iter}].val)
Base.last(c::AbstractChains) = last(c.value[Axis{:iter}].val)

#################### Auxilliary Functions ####################

function combine(c::AbstractChains)
  n, p, m = size(c.value)
  value = Array{Float64}(undef, n * m, p)
  for j in 1:p
    idx = 1
    for i in 1:n, k in 1:m
      value[idx, j] = c.value[i, j, k]
      idx += 1
    end
  end
  value
end

function Base.range(c::AbstractChains)
    return c.value[Axis{:iter}].val
end

function chains(c::AbstractChains)
    return c.value[Axis{:chain}].val
end

function Base.names(c::AbstractChains)
    return c.value[Axis{:var}].val
end

# Return a new chain for each section.
function get_sections(c::AbstractChains, sections::Vector = [])
    sections = length(sections) == 0 ? collect(keys(c.name_map)) : sections
    return [Chains(c, section) for section in sections]
end

# Return a new chain for each section.
function get_sections(c::AbstractChains, section::Union{Symbol, String})
    return get_sections(c, [section])
end

function header(c::AbstractChains)
    rng = range(c)
    sections = []
    for (section, nms) in c.name_map
        push!(sections, string("$section",
            repeat(" ", 18 - length(string(section))),
            "= $(join(map(string, nms), ", "))\n"
        ))
    end
    string(
        # "Log model evidence = $(c.logevidence)\n", #FIXME: Uncomment.
        "Iterations        = $(first(c)):$(last(c))\n",
        "Thinning interval = $(step(c))\n",
        "Chains            = $(join(map(string, chains(c)), ", "))\n",
        "Samples per chain = $(length(range(c)))\n",
        sections...
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

function _trim_name_map(names::Vector, name_map::Dict)
    n = copy(name_map)
    for (key, values) in name_map
        intersection = values ∩ names
        if length(intersection) > 0
            n[key] = intersection
        else
            delete!(n, key)
        end
    end
    return n
end

### Chains specific functions ###
function sort(c::Chains{T}) where T<:Real
    v = c.value
    x, y, z = size(v)
    unsorted = collect(zip(1:y, v.axes[2].val))
    sorted = sort(unsorted, by = x -> string(x[2]), lt=MCMCChain.natural)
    new_axes = (v.axes[1], Axis{:var}([n for (_, n) in sorted]), v.axes[3])
    new_v = copy(v.data)
    for i in eachindex(sorted)
        new_v[:, i, :] = v[:, sorted[i][1], :]
    end
    aa = AxisArray(new_v, new_axes...)
    return MCMCChain.Chains{T}(aa, c.logevidence, c.name_map, c.info)
end
