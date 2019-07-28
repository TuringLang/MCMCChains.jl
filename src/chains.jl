#################### Chains ####################

## Constructors ##

# Default name map.
const DEFAULT_MAP = Dict{Symbol, Vector{Any}}(:parameters => [])

# Set default parameter names if not given.
# function Chains(val::AbstractArray{A,3};
#         start::Int=1,
#         thin::Int=1,
#         evidence = missing,
#         info=NamedTuple()) where {A<:Union{Real, Union{Missing, Real}}}
#
#     return Chains(val, parameter_names, start=start, thin=thin)
# end
#
# Constructor to handle a vector of vectors.
function Chains(val::Vector{Vector{A}},
				parameter_names::Vector{String} = map(i->"Param$i", 1:length(first(val))),
        name_map = copy(DEFAULT_MAP);
        start::Int=1,
        thin::Int=1,
        evidence = missing,
        info::NamedTuple=NamedTuple()
) where {A<:Union{Real, Union{Missing, Real}}}
	println("It's working")
	return Chains(Array(hcat(val...)'), parameter_names, name_map, start=start,
           thin=thin, evidence=evidence, info=info)
end

# Constructor to handle a 1D array.
function Chains(val::AbstractArray{A,1},
        parameter_names::Vector{String} = map(i->"Param$i", 1:size(val, 2)),
        name_map = copy(DEFAULT_MAP);
        start::Int=1,
        thin::Int=1,
        evidence = missing,
        info::NamedTuple=NamedTuple()
) where {A<:Union{Real, Union{Missing, Real}}}
	return Chains(val[:,:,:], parameter_names, name_map, start=start,
           thin=thin, evidence=evidence, info=info)
end

# Set default parameter names if not given.
function Chains(val::AbstractArray{A,2},
        parameter_names::Vector{String} = map(i->"Param$i", 1:size(val, 2)),
        name_map = copy(DEFAULT_MAP);
        start::Int=1,
        thin::Int=1,
        evidence = missing,
        info::NamedTuple=NamedTuple()
) where {A<:Union{Real, Union{Missing, Real}}}
    Chains(val[:,:,:], parameter_names, name_map, start=start,
           thin=thin, evidence=evidence, info=info)
end

# Generic chain constructor.
function Chains(val::AbstractArray{A,3},
        parameter_names::Vector{String} = map(i->"Param$i", 1:size(val, 2)),
        name_map = copy(DEFAULT_MAP);
        start::Int=1,
        thin::Int=1,
        evidence = missing,
        info::NamedTuple=NamedTuple()) where {A<:Union{Real, Union{Missing, Real}}}

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
        collect(1:size(val, 3)),
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
        empty_df_vec = [ChainDataFrame("", DataFrame())]
        s = (hash(0), empty_df_vec)
        info = merge(info, (hashedsummary = Ref(s),))
    end

    # Construct the AxisArray.
    axs = ntuple(i -> Axis{names[i]}(axvals[i]), 3)
    arr = AxisArray(val, axs...)
    return sort(
        Chains{A, typeof(evidence), typeof(name_map_tupl), typeof(info)}(
            arr, evidence, name_map_tupl, info)
    )
end

# Retrieve a new chain with only a specific section pulled out.
function Chains(c::Chains{A, T, K, L}, section::Union{Vector, Any};
                sorted=false) where {A, T, K, L}
    section = typeof(section) <: AbstractArray ? section : [section]

    # If we received an empty list, return the original chain.
    if isempty(section)
        if sorted
            return sort(c)
        else
            return c
        end
    end

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
function _sym2index(c::Chains, v::Union{Vector{Symbol}, Vector{String}}; sort = true)
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
    return sort ? sort!(syms, lt=MCMCChains.natural) : syms
end

Base.getindex(c::Chains, i1::T) where T<:Union{AbstractUnitRange, StepRange} = c[i1, :, :]
Base.getindex(c::Chains, i1::Integer) = c[i1:i1, :, :]
Base.getindex(c::Chains, v::Symbol) = c[[v]]
Base.getindex(c::Chains, v::String) = c[:, [v], :]
Base.getindex(c::Chains, v::Vector{String}) = c[:, v, :]

function Base.getindex(c::Chains, v::Vector{Symbol})
    syms = _sym2index(c, v)
    return c[:, syms, :]
end

function Base.getindex(c::Chains{A, T, K, L}, i...) where {A, T, K, L}
    # Make sure things are in array form to preserve the axes.
    ind = [typeof(i[1]) <: Integer ? (i[1]:i[1]) : i[1],
           typeof(i[2]) <: Union{AbstractArray, Colon} ?  i[2] : [i[2]],
           typeof(i[3]) <: Union{AbstractArray, Colon} ?  i[3] : [i[3]]
    ]

    # Check to see if we received a symbol in i[2].
    if ind[2] != Colon() && typeof(first(ind[2])) <: Symbol
        ind[2] = _sym2index(c, ind[2])
    end

    newval = getindex(c.value, ind...)
    names = newval.axes[2].val
    new_name_map = _trim_name_map(names, c.name_map)
    return Chains{A, T, typeof(new_name_map), L}(newval,
        c.logevidence,
        new_name_map,
        c.info)
end

Base.setindex!(c::Chains, v, i...) = setindex!(c.value, v, i...)
Base.lastindex(c::Chains) = lastindex(c.value, 1)
Base.lastindex(c::Chains, d::Integer) = lastindex(c.value, d)

"""
    Base.get(c::Chains, v::Symbol; flatten=false)
    Base.get(c::Chains, vs::Vector{Symbol}; flatten=false)

Returns a `NamedTuple` with `v` as the key, and matching paramter
names as the values.

Passing `flatten=true` will return a `NamedTuple` with keys ungrouped.

Example:

```julia
x = get(c, :param1)
x = get(c, [:param1, :param2])
```
"""
Base.get(c::Chains, v::Symbol; flatten = false) = get(c, [v], flatten=flatten)
function Base.get(c::Chains, vs::Vector{Symbol}; flatten = false)
    pairs = Dict()
    for v in vs
        syms = _sym2index(c, [v])
        len = length(syms)
        val = ()
        if len > 1
            val = ntuple(i -> c.value[:,syms[i],:], length(syms))
        elseif len == 1
            val = c.value[:,syms[1],:]
        else
            continue
        end

        if flatten
            for i in eachindex(syms)
                pairs[syms[i]] = val[i]
            end
        else
            pairs[v] = val
        end
    end
    return _dict2namedtuple(pairs)
end

"""
    get(c::Chains; section::Union{Vector{Symbol}, Symbol; flatten=false}

Returns all parameters in a given section(s) as a `NamedTuple`.

Passing `flatten=true` will return a `NamedTuple` with keys ungrouped.

Example:

```julia
x = get(chn, section = :parameters)
x = get(chn, section = [:internals, :parameters])
```
"""
function Base.get(c::Chains;
        section::Union{Vector{Symbol}, Symbol},
        flatten = false)
    section = section isa Symbol ? [section] : section
    not_found = Symbol[]
    names = Set(String[])
    for v in section
        if v in keys(c.name_map)
            # If the name contains a bracket,
            # split it so get can group them correctly.
            nms = flatten ?
                c.name_map[v] :
                map(n -> String(split(n, "[")[1]), c.name_map[v])
            push!(names, nms...)
        else
            push!(not_found, v)
        end
    end

    if length(not_found) > 0
        throw(ArgumentError("$not_found not found in chains name map."))
    end

    return get(c, Symbol.(names), flatten = flatten)
end

"""
    get_params(c::Chains; flatten=false)

Returns all parameters packaged as a `NamedTuple`. Variables with a bracket
in their name (as in "P[1]") will be grouped into the returned value as P.

Passing `flatten=true` will return a `NamedTuple` with keys ungrouped.

Example:

```julia
x = get_params(chn)
x.P
```
"""
get_params(c::Chains; flatten = false) = get(c, section = sections(c), flatten=flatten)

#################### Base Methods ####################

function Base.show(io::IO, c::Chains)
    print(io, "Object of type Chains, with data of type $(summary(c.value.data))\n\n")
    println(io, header(c))

    # Grab the value hash.
    h = hash(c)

    if :hashedsummary in keys(c.info)
        s = c.info.hashedsummary.x
        if s[1] == h
            show(io, s[2])
        else
            new_summary = describe(c)
            c.info.hashedsummary.x = (h, new_summary)
            show(io, new_summary)
        end
    else
        show(io, describe(c, suppress_header=true))
    end
end

Base.keys(c::AbstractChains) = names(c)
Base.size(c::AbstractChains) = size(c.value)
Base.size(c::AbstractChains, ind) = size(c)[ind]
Base.length(c::AbstractChains) = length(range(c))
Base.first(c::AbstractChains) = first(c.value[Axis{:iter}].val)
Base.step(c::AbstractChains) = step(c.value[Axis{:iter}].val)
Base.last(c::AbstractChains) = last(c.value[Axis{:iter}].val)

Base.convert(::Type{Array}, chn::MCMCChains.Chains) = convert(Array, chn.value)

#################### Auxilliary Functions ####################

function Base.hash(c::Chains)
    val = hash(c.value) + hash(c.info) + hash(c.name_map) + hash(c.logevidence)
    return hash(val)
end

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
    names(c::AbstractChains, sections)

Return the parameter names in a `Chains` object.
"""
function names(c::AbstractChains)
    return c.value[Axis{:var}].val
end

"""
    names(c::AbstractChains, sections::Union{Symbol, Vector{Symbol}})

Return the parameter names in a `Chains` object, given an array of sections.
"""
function names(c::AbstractChains,
    sections::Union{Symbol, Vector{Symbol}})
    # Check that sections is an array.
    sections = typeof(sections) <: AbstractArray ?
        sections :
        [sections]

    nms = []

    for i in sections
        push!(nms, c.name_map[i]...)
    end
    return nms
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
    sections(c::AbstractChains)

Retrieve a list of the sections in a chain.
"""
sections(c::AbstractChains) = collect(keys(c.name_map))

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

    #

    # Return header.
    return string(
        ismissing(c.logevidence) ? "" : "Log evidence      = $(c.logevidence)\n",
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
  cc = copy(c.value.data)
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
    for (key, values) in n
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

set_section(c::Chains, nt::NamedTuple) = set_section(c, _namedtuple2dict(nt))

"""
    set_section(c::Chains, nt::Dict)

Changes a chains name mapping to a provided dictionary. This also supports a NamedTuple.
Any parameters in the chain that are unassigned will be placed into
the :parameters section.
"""
function set_section(c::Chains{A, T, K, L}, d::Dict) where {A,T,K,L}
    # Add :parameters if it's not there.
    if !(:parameters in keys(d))
        d[:parameters] = []
    end

    # Make sure all the names are in the new name map.
    nms = Set([])
    for values in values(d)
        for val in values
            push!(nms, val)
        end
    end
    missing_names = setdiff(names(c), nms)

    # Assign everything to :parameters if anything's missing.
    if length(missing_names) > 0
        @warn "Section mapping does not contain all parameter names, " *
            "$missing_names assigned to :parameters."
        push!(d[:parameters], missing_names...)
    end

    nt = _dict2namedtuple(d)
    return Chains{A, T, typeof(nt), L}(
        c.value,
        c.logevidence,
        nt,
        c.info
    )
end

function _use_showall(c::AbstractChains, section::Symbol)
    if section == :parameters && !in(:parameters, keys(c.name_map))
        return true
    end
    return false
end

function _clean_sections(c::AbstractChains, sections::Union{Vector{Symbol}, Symbol})
    sections = sections isa AbstractArray ? sections : [sections]
    ks = collect(keys(c.name_map))
    return ks ∩ sections
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

chainscat(c1::AbstractChains, args::AbstractChains...) = cat(c1, args..., dims=3)
Base.hcat(c1::AbstractChains, args::AbstractChains...) = cat(c1, args..., dims=2)
Base.vcat(c1::AbstractChains, args::AbstractChains...) = cat(c1, args..., dims=1)

function pool_chain(c::Chains{A, T, K, L}) where {A, T, K, L}
    val = c.value.data
    concat = vcat([val[:,:,j] for j in 1:size(val,3)]...)
    return Chains(cat(concat, dims=3), names(c), c.name_map; info=c.info)
end
