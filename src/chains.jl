#################### Chains ####################

## Constructors ##

# Default name map.
const DEFAULT_MAP = Dict{Symbol, Vector{Any}}(:parameters => [])

# Set default parameter names if not given.
function Chains(val::AbstractArray{A,3};
        start::Int=1,
        thin::Int=1,
        evidence = 0.0,
        info=NamedTuple()) where {A<:Union{Real, Union{Missing, Real}}}
    parameter_names = map(i->"Param$i", 1:size(val, 2))
    return Chains(val, parameter_names, start=start, thin=thin)
end

# Generic chain constructor.
function Chains(val::AbstractArray{A,3},
        parameter_names::Vector,
        name_map = copy(DEFAULT_MAP);
        start::Int=1,
        thin::Int=1,
        evidence = 0.0,
        info=NamedTuple()) where {A<:Union{Real, Union{Missing, Real}}}

    # If we received an array of pairs, convert it to a dictionary.
    if typeof(name_map) <: Array
        name_map = Dict(name_map)
    elseif typeof(name_map) <: NamedTuple
        name_map = _namedtuple2dict(name_map)
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

    # Construct axis names and ranges.
    names = [:iter, :var, :chain]
    axvals = [
        Base.range(start, step=thin, length=size(val, 1)),
        parameter_names,
        map(i->Symbol("Chain$i"), 1:size(val, 3)),
    ]

    if length(keys(name_map)) == 1
        name_map[first(keys(name_map))] = parameter_names
    else
        # Store unassigned variables.
        unassigned = Set([])

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
                push!(unassigned, param)
            end
        end

        # Assign all unassigned parameter names.
        for param in unassigned
            push!(name_map[:parameters], param)
        end
    end

    # Make the name_map NamedTuple.
    name_map_tupl = _dict2namedtuple(name_map)

    # Ensure that we have a hashedsummary key in info.
    if !in(:hashedsummary, keys(info))
        s = (hash(0), ChainSummaries("", []))
        info = merge(info, (hashedsummary = Ref(s),))
    end

    # Construct the AxisArray.
    axs = ntuple(i -> Axis{names[i]}(axvals[i]), 3)
    arr = AxisArray(convert(Array{Union{Missing,A},3}, val), axs...)
    return sort(
        Chains{A, typeof(evidence), typeof(name_map_tupl), typeof(info)}(
            arr, evidence, name_map_tupl, info)
    )
end

# Retrieve a new chain with only a specific section pulled out.
function Chains(c::Chains{A, T, K, L}, section::Union{Vector, Any};
                sorted=false) where {A, T<:Real, K, L}
    section = typeof(section) <: AbstractArray ? section : [section]

    # If we received an empty list, return the original chain.
    if isempty(section) return c end

    # Gather the names from the relevant sections.
    names = []
    if isa(section, Vector)
        for s in section
            # Make sure the section exists first.
            in(s, keys(c.name_map)) ||
                throw(ArgumentError("$section not found in Chains name map."))

            names = vcat(names, c.name_map[s])
        end
    else
        in(s, keys(c.name_map)) ||
            throw(ArgumentError("$section not found in Chains name map."))

        names = c.name_map[s]
    end

    # Extract wanted values.
    new_vals = c.value[:, names, :]

    # Create the new section map.
    new_section_map = _dict2namedtuple(Dict([s => c.name_map[s] for s in section]))

    # Create the new chain.
    new_chn = Chains{A, T, typeof(new_section_map), L}(new_vals,
        c.logevidence,
        new_section_map,
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
Base.getindex(c::Chains, v::String) = Array(c.value[:, [v], :])
Base.getindex(c::Chains, v::Vector{String}) = Array(c.value[:, v, :])

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

    sort!(syms, lt=MCMCChains.natural)
    return Array(c.value[:, syms, :])
end

function Base.getindex(c::Chains{A, T, K, L}, i...) where {A, T, K, L}
    # Make sure things are in array form to preserve the axes.
    ind = (i[1],
           typeof(i[2]) <: Union{AbstractArray, Colon} ?  i[2] : [i[2]],
           typeof(i[3]) <: Union{AbstractArray, Colon} ?  i[3] : [i[3]]
    )

    newval = getindex(c.value, ind...)
    names = newval.axes[2].val
    new_name_map = _trim_name_map(names, c.name_map)
    return Chains{A, T, typeof(new_name_map), L}(newval,
        c.logevidence,
        new_name_map,
        c.info)
end

Base.setindex!(c::Chains, v, i...) = setindex!(c.value, v, i...)


#################### Base Methods ####################

function Base.show(io::IO, c::Chains)
    print(io, "Object of type Chains, with data of type $(summary(c.value.data))\n\n")
    println(io, header(c))

    # Grab the value hash.
    h = hash(c.value)

    if :hashedsummary in keys(c.info)
        s = c.info.hashedsummary.x
        if s[1] == h
            show(io, s[2])
        else
            new_summary = summarystats(c, suppress_header=true)
            c.info.hashedsummary.x = (h, new_summary)
            show(io, new_summary)
        end
    else
        show(io, summarystats(c, suppress_header=true))
    end
end

function Base.size(c::AbstractChains)
  dim = size(c.value)
  last(c), dim[2], dim[3]
end

Base.keys(c::AbstractChains) = names(c)
Base.size(c::AbstractChains, ind) = size(c)[ind]
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

"""
    range(c::AbstractChains)

Returns the range used in a `Chains` object.
"""
function range(c::AbstractChains)
    return c.value[Axis{:iter}].val
end

"""
    chains(c::AbstractChains)

Returns the names or symbols of each chain in an `AbstractChains` object.
"""
function chains(c::AbstractChains)
    return c.value[Axis{:chain}].val
end

"""
    names(c::AbstractChains)

Return the parameter names in a `Chains` object.
"""
function names(c::AbstractChains)
    return c.value[Axis{:var}].val
end

"""
    get_sections(c::AbstractChains, sections::Vector = [])

Returns multiple `Chains` objects, each containing only a single section.
"""
function get_sections(c::AbstractChains, sections::Vector = [])
    sections = length(sections) == 0 ? collect(keys(c.name_map)) : sections
    return [Chains(c, section) for section in sections]
end

# Return a new chain for each section.
function get_sections(c::AbstractChains, section::Union{Symbol, String})
    return get_sections(c, [section])
end

"""
    header(c::Chains; section=missing)

Returns a string containing summary information for a `Chains` object.
If the `section` keyword is used, this function prints only the relevant section
header.

Example:

```julia
# Printing the whole header.
header(chn)

# Print only one section's header.
header(chn, section = :parameter)
```
"""
function header(c::AbstractChains; section=missing)
    rng = range(c)

    # Function to make section strings.
    section_str(sec, arr) = string(
        "$sec",
        repeat(" ", 18 - length(string(sec))),
        "= $(join(map(string, arr), ", "))\n"
    )

    # Set up string array.
    section_strings = String[]

    # Get section lines.
    if section isa Missing
        for (sec, nms) in pairs(c.name_map)
            section_string = section_str(sec, nms)
            push!(section_strings, section_string)
        end
    else
        section in keys(c.name_map) ||
            throw(ArgumentError("$section not found in name map."))
        section_string = section_str(section, c.name_map[section])
        push!(section_strings, section_string)
    end

    # Return header.
    return string(
        # "Log model evidence = $(c.logevidence)\n", #FIXME: Uncomment.
        "Iterations        = $(first(c)):$(last(c))\n",
        "Thinning interval = $(step(c))\n",
        "Chains            = $(join(map(string, chains(c)), ", "))\n",
        "Samples per chain = $(length(range(c)))\n",
        section_strings...
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

"""
    _trim_name_map(names::Vector, name_map::NamedTuple)

This is an internal function used to remove values from a name map
and return a new name_map.
"""
function _trim_name_map(names::Vector, name_map::NamedTuple)
    n = _namedtuple2dict(name_map)
    for (key, values) in name_map
        intersection = values ∩ names
        if length(intersection) > 0
            n[key] = intersection
        else
            delete!(n, key)
        end
    end
    return _dict2namedtuple(n)
end

### Chains specific functions ###
"""
    sort(c::Chains)

Returns a new column-sorted version of `c`, using natural sort order.
"""
function sort(c::Chains{A, T, K, L}) where {A, T, K, L}
    v = c.value
    x, y, z = size(v)
    unsorted = collect(zip(1:y, v.axes[2].val))
    sorted = sort(unsorted, by = x -> string(x[2]), lt=MCMCChains.natural)
    new_axes = (v.axes[1], Axis{:var}([n for (_, n) in sorted]), v.axes[3])
    new_v = copy(v.data)
    for i in eachindex(sorted)
        new_v[:, i, :] = v[:, sorted[i][1], :]
    end
    aa = AxisArray(new_v, new_axes...)
    return MCMCChains.Chains{A, T, K, L}(aa, c.logevidence, c.name_map, c.info)
end

"""
    setinfo(c::Chains, n::NamedTuple)

Returns a new `Chains` object with a `NamedTuple` type `n` placed in the `info` field.

Example:

```julia
new_chn = setinfo(chn, NamedTuple{(:a, :b)}((1, 2)))
```
"""
function setinfo(c::Chains{A, T, K}, n::NamedTuple) where {A, T, K}
    return Chains{A, T, K, typeof(n)}(
        c.value,
        c.logevidence,
        c.name_map,
        n
    )
end


#################### Concatenation ####################

function Base.cat(c1::AbstractChains, args::AbstractChains...; dims::Integer = 3)
  dims == 1 ? cat1(c1, args...) :
  dims == 2 ? cat2(c1, args...) :
  dims == 3 ? cat3(c1, args...) :
    throw(ArgumentError("cannot concatenate along dimension $dim"))
end

function cat1(c1::AbstractChains, args::AbstractChains...)
    rng = range(c1)
    for c in args
        last(rng) + step(rng) == first(c) ||
            throw(ArgumentError("noncontiguous chain iterations"))
        step(rng) == step(c) ||
            throw(ArgumentError("chain thinning differs"))
        rng = first(rng):step(rng):last(c)
    end

    nms = names(c1)
    all(c -> names(c) == nms, args) ||
        throw(ArgumentError("chain names differ"))

    chns = chains(c1)
    all(c -> chains(c1) == chns, args) ||
        throw(ArgumentError("sets of chains differ"))

    name_map = _ntdictmerge(c1.name_map, map(c -> c.name_map, args)...)

    value = cat(c1.value, map(c -> c.value, args)..., dims=1)
    Chains(value, nms, name_map, start=first(rng), thin=step(rng),
        info = c1.info)
end

function cat2(c1::AbstractChains, args::AbstractChains...)
  rng = range(c1)
  all(c -> range(c) == rng, args) ||
    throw(ArgumentError("chain ranges differ"))

  nms = names(c1)
  n = length(nms)
  for c in args
    nms = union(nms, names(c))
    n += length(names(c))
    n == length(nms) ||
      throw(ArgumentError("non-unique parameter names"))
  end

  name_map = _ntdictmerge(c1.name_map, map(c -> c.name_map, args)...)

  chns = chains(c1)
  all(c -> chains(c) == chns, args) ||
    throw(ArgumentError("sets of chains differ"))

  value = cat(c1.value, map(c -> c.value, args)..., dims=2)
  Chains(value, nms, name_map, start=first(rng), thin=step(rng),
      info = c1.info)
end

function cat3(c1::AbstractChains, args::AbstractChains...)
  rng = range(c1)
  all(c -> range(c) == rng, args) ||
    throw(ArgumentError("chain ranges differ"))

  nms = names(c1)

  name_map = _ntdictmerge(c1.name_map, map(c -> c.name_map, args)...)

  value = cat(c1.value.data, map(c -> c.value.data, args)..., dims=3)
  Chains(value, nms, name_map, start=first(rng), thin=step(rng),
      info = c1.info)
end

Base.hcat(c1::AbstractChains, args::AbstractChains...) = cat(2, c1, args...)
# Base.vcat(c1::AbstractChains, args::AbstractChains...) = cat(1, c1, args...)

function Base.vcat(cs::Chains...)
    c1 = cs[1]

    evidence = c1.logevidence
    nms = names(c1)
    rng = range(c1)
    chns = chains(c1)

    value = c1.value
    info = c1.info
    name_map = c1.name_map

    all(c -> names(c) == nms, cs) ||
        throw(ArgumentError("parameter names differ"))

    all(c -> chains(c) == chns, cs) ||
        throw(ArgumentError("sets of chains differ"))

    for c in cs[2:end]
        @assert evidence == c.logevidence
        @assert rng == range(c)

        value = cat(value, c.value, dims=1)
        info = _ntdictmerge(info, map(c -> c.info, args)...)
        name_map = _ntdictmerge(name_map, map(c -> c.name_map, args)...)
    end

    chn = Chains(value,
                nms,
                name_map,
                info=info,
                evidence=evidence)
    return chn
end
