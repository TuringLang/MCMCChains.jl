#################### Chains ####################

## Constructors ##

# Constructor to handle a vector of vectors.
Chains(val::AbstractVector{<:AbstractVector{<:Union{Missing, Real}}}, args...; kwargs...) =
	Chains(copy(reduce(hcat, val)'), args...; kwargs...)

# Constructor to handle a 1D array.
Chains(val::AbstractVector{<:Union{Missing, Real}}, args...; kwargs...) =
	Chains(reshape(val, :, 1, 1), args...; kwargs...)

# Constructor to handle a 2D array
Chains(val::AbstractMatrix{<:Union{Missing, Real}}, args...; kwargs...) =
    Chains(reshape(val, size(val, 1), size(val, 2), 1), args...; kwargs...)

# Constructor to handle parameter names that are not Symbols.
function Chains(
    val::AbstractArray{<:Union{Missing,Real},3},
    parameter_names::AbstractVector,
    args...;
    kwargs...
)
    return Chains(val, Symbol.(parameter_names), args...; kwargs...)
end

# Generic chain constructor.
function Chains(
    val::AbstractArray{<:Union{Missing, Real},3},
    parameter_names::AbstractVector{Symbol} = Symbol.(:param_, 1:size(val, 2)),
    name_map = (parameters = parameter_names,);
    start::Int = 1,
    thin::Int = 1,
    evidence = missing,
    info::NamedTuple = NamedTuple()
)
    # Make sure that we have a `:parameters` index and # Copying can avoid state mutation.
    _name_map = initnamemap(name_map)

    # Preclean the name_map of names that aren't in the
    # parameter_names vector.
    for names in _name_map
        filter!(x -> x ∈ parameter_names, names)
    end

    # Store unassigned variables.
    unassigned = Set(Symbol[])

    # Check that all parameters are assigned.
    for param in parameter_names
        if all(param ∉ names for names in _name_map)
            push!(unassigned, param)
        end
    end

    # Assign all unassigned parameter names.
    append!(_name_map[:parameters], unassigned)

    # Construct the AxisArray.
    arr = AxisArray(val;
                    iter = range(start, step=thin, length=size(val, 1)),
                    var = parameter_names,
                    chain = 1:size(val, 3))

    # Create the new chain.
    return Chains(arr, evidence, _name_map, info)
end

# Retrieve a new chain with only a specific section pulled out.
Chains(c::Chains, section::Union{Symbol,String}) = Chains(c, (section,))
function Chains(chn::Chains, sections)
    # Make sure the sections exist first.
	  all(haskey(chn.name_map, Symbol(x)) for x in sections) ||
		    error("some sections are not present in the chain")

	  # Create the new section map.
    name_map = (; (Symbol(k) => chn.name_map[Symbol(k)] for k in sections)...)

    # Extract wanted values.
    value = chn.value[:, reduce(vcat, name_map), :]

    # Create the new chain.
    return Chains(value, chn.logevidence, name_map, chn.info)
end
Chains(chain::Chains, ::Nothing) = chain

# Groups of parameters

"""
    namesingroup(chains::Chains, sym::Union{String,Symbol})

Return the names of all parameters in a chain that belong to the group `:name`.

This is based on the MCMCChains convention that parameters with names of the form `:name[index]` 
belong to one group of parameters called `:name`.

If the chain contains a parameter of name `:name` it will be returned as well.
"""
namesingroup(chains::Chains, sym::String) = namesingroup(chains, Symbol(sym))
function namesingroup(chains::Chains, sym::Symbol)
    # Start by looking up the symbols in the list of parameter names.
    names_of_params = names(chains)
    regex = Regex("^$sym\$|^$sym\\[")
    indices = findall(x -> match(regex, string(x)) !== nothing, names(chains))
    return names_of_params[indices]
end

"""
    group(chains::Chains, name::Union{String,Symbol})

Return a subset of the chain chain with all parameters in the group `Symbol(name)`.
"""
function group(chains::Chains, name::Union{String,Symbol})
    return chains[:, namesingroup(chains, name), :]
end

#################### Indexing ####################

Base.getindex(c::Chains, i::Integer) = c[i, :, :]
Base.getindex(c::Chains, i::AbstractVector{<:Integer}) = c[i, :, :]

Base.getindex(c::Chains, v::String) = c[:, Symbol(v), :]
Base.getindex(c::Chains, v::AbstractVector{String}) = c[:, Symbol.(v), :]

Base.getindex(c::Chains, v::Symbol) = c[:, v, :]
Base.getindex(c::Chains, v::AbstractVector{Symbol}) = c[:, v, :]

Base.getindex(chn::Chains, i, j, k) = _getindex(chn, chn.value[_toindex(i, j, k)...])

_getindex(::Chains, data) = data
function _getindex(chains::Chains, data::AxisArray{<:Any,3})
    names = data.axes[2].val
    namemap = namemap_intersect(chains.name_map, names)
    return Chains(data, chains.logevidence, namemap, chains.info)
end

# convert strings to symbols but try to keep all dimensions for multiple parameters
_toindex(i, j, k) = (i, string2symbol(j), k)
_toindex(i::Integer, j, k) = (i:i, string2symbol(j), k)
_toindex(i, j, k::Integer) = (i, string2symbol(j), k:k)
_toindex(i::Integer, j, k::Integer) = (i:i, string2symbol(j), k:k)

# return an array or a number if a single parameter is specified
const SingleIndex = Union{Symbol,String,Integer}
_toindex(i, j::SingleIndex, k) = (i, string2symbol(j), k)
_toindex(i::Integer, j::SingleIndex, k) = (i, string2symbol(j), k)
_toindex(i, j::SingleIndex, k::Integer) = (i, string2symbol(j), k)
_toindex(i::Integer, j::SingleIndex, k::Integer) = (i, string2symbol(j), k)

Base.setindex!(c::Chains, v, i...) = setindex!(c.value, v, i...)
Base.lastindex(c::Chains) = lastindex(c.value, 1)
Base.lastindex(c::Chains, d::Integer) = lastindex(c.value, d)

"""
    Base.get(c::Chains, v::Symbol; flatten=false)
    Base.get(c::Chains, vs::Vector{Symbol}; flatten=false)

Return a `NamedTuple` with `v` as the key, and matching paramter
names as the values.

Passing `flatten=true` will return a `NamedTuple` with keys ungrouped.

Example:

```julia
x = get(c, :param_1)
x = get(c, [:param_1, :param_2])
```
"""
Base.get(c::Chains, v::Symbol; flatten = false) = get(c, [v], flatten=flatten)
function Base.get(c::Chains, vs::Vector{Symbol}; flatten = false)
    pairs = Dict()
    for v in vs
        syms = namesingroup(c, v)
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

Return all parameters in a given section(s) as a `NamedTuple`.

Passing `flatten=true` will return a `NamedTuple` with keys ungrouped.

Example:

```julia
x = get(chn, section = :parameters)
x = get(chn, section = [:internals, :parameters])
```
"""
function Base.get(
    c::Chains;
    section::Union{Symbol,AbstractVector{Symbol}},
    flatten = false
)
    names = Set(Symbol[])
    regex = r"[^\[]*"
    _section = section isa Symbol ? (section,) : section
    for v in _section
        v in keys(c.name_map) || error("section $v does not exist")

        # If the name contains a bracket,
        # split it so get can group them correctly.
        if flatten
            append!(names, c.name_map[v])
        else
            for name in c.name_map[v]
                m = match(regex, string(name))
                push!(names, Symbol(m.match))
            end
        end
    end

    return get(c, collect(names); flatten = flatten)
end

"""
    get_params(c::Chains; flatten=false)

Return all parameters packaged as a `NamedTuple`. Variables with a bracket
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

function Base.show(io::IO, chains::Chains)
    print(io, "MCMC chain (", summary(chains.value.data), ")")
end

function Base.show(io::IO, mime::MIME"text/plain", chains::Chains)
    print(io, "Chains ", chains, ":\n\n", header(chains))

    # Show summary stats.
    summaries = describe(chains)
    for summary in summaries
        println(io)
        show(io, mime, summary)
    end
end

Base.keys(c::Chains) = names(c)
Base.size(c::Chains) = size(c.value)
Base.size(c::Chains, ind) = size(c)[ind]
Base.length(c::Chains) = length(range(c))
Base.first(c::Chains) = first(c.value[Axis{:iter}].val)
Base.step(c::Chains) = step(c.value[Axis{:iter}].val)
Base.last(c::Chains) = last(c.value[Axis{:iter}].val)

Base.convert(::Type{Array}, chn::Chains) = convert(Array, chn.value)

#################### Auxilliary Functions ####################

"""
    range(chains::Chains)

Return the range of iteration indices of the `chains`.
"""
Base.range(chains::Chains) = chains.value[Axis{:iter}].val

"""
    setrange(chains::Chains, range)

Generate a new chain from `chains` with iterations indexed by `range`.

The new chain and `chains` share the same data in memory.
"""
function setrange(chains::Chains, range::AbstractRange{<:Integer})
    if length(chains) != length(range)
        error("length of `range` (", length(range),
              ") is not equal to the number of iterations (", length(chains), ")")
    end

    value = AxisArray(chains.value.data;
                      iter = range, var = names(chains), chain = MCMCChains.chains(chains))

    return Chains(value, chains.logevidence, chains.name_map, chains.info)
end

"""
    resetrange(chains::Chains)

Generate a new chain from `chains` with iterations indexed by `1:n`, where `n` is the number
of samples per chain.

The new chain and `chains` share the same data in memory.
"""
resetrange(chains::Chains) = setrange(chains, 1:size(chains, 1))

"""
    chains(c::Chains)

Return the names or symbols of each chain in a `Chains` object.
"""
function chains(c::Chains)
    return c.value[Axis{:chain}].val
end

"""
    names(chains::Chains)

Return the parameter names in the `chains`.
"""
Base.names(chains::Chains) = chains.value[Axis{:var}].val

"""
    names(chains::Chains, section::Symbol)

Return the parameter names of a `section` in the `chains`.
"""
Base.names(chains::Chains, section::Symbol) = convert(Vector{Symbol}, chains.name_map[section])

"""
    names(chains::Chains, sections)

Return the parameter names of the `sections` in the `chains`.
"""
function Base.names(c::Chains, sections)
    names = Symbol[]
    for section in sections
        append!(names, c.name_map[section])
    end
    return names
end

"""
    get_sections(chains[, sections])

Return multiple `Chains` objects, each containing only a single section.
"""
function get_sections(chains::Chains, sections = keys(chains.name_map))
    return [Chains(chains, section) for section in sections]
end
get_sections(chains::Chains, section::Union{Symbol, String}) = Chains(chains, section)

"""
    sections(c::Chains)

Retrieve a list of the sections in a chain.
"""
sections(c::Chains) = collect(keys(c.name_map))

"""
    header(c::Chains; section=missing)

Return a string containing summary information for a `Chains` object.
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
function header(c::Chains; section=missing)
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

function indiscretesupport(c::Chains,
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

function link(c::Chains)
  cc = copy(c.value.data)
  for j in axes(cc, 2)
    x = cc[:, j, :]
    if minimum(x) > 0.0
      cc[:, j, :] = maximum(x) < 1.0 ? logit.(x) : log.(x)
    end
  end
  cc
end

### Chains specific functions ###

"""
    sort(c::Chains[; lt=NaturalSort.natural])

Return a new column-sorted version of `c`.

By default the columns are sorted in natural sort order.
"""
function Base.sort(c::Chains; lt = NaturalSort.natural)
    v = c.value
    x, y, z = size(v)
    unsorted = collect(zip(1:y, v.axes[2].val))
    sorted = sort(unsorted, by = x -> string(x[2]), lt=lt)
    new_axes = (v.axes[1], Axis{:var}([n for (_, n) in sorted]), v.axes[3])
    new_v = copy(v.data)
    for i in eachindex(sorted)
        new_v[:, i, :] = v[:, sorted[i][1], :]
    end

    aa = AxisArray(new_v, new_axes...)

    # Sort the name map too:
    namemap = deepcopy(c.name_map)
    for names in namemap
        sort!(names, by=string, lt=lt)
    end

    return Chains(aa, c.logevidence, namemap, c.info)
end

"""
    setinfo(c::Chains, n::NamedTuple)

Return a new `Chains` object with a `NamedTuple` type `n` placed in the `info` field.

Example:

```julia
new_chn = setinfo(chn, NamedTuple{(:a, :b)}((1, 2)))
```
"""
function setinfo(c::Chains, n::NamedTuple)
    return Chains(c.value, c.logevidence, c.name_map, n)
end

"""
    set_section(chains::Chains, namemap)

Create a new `Chains` object from `chains` with the provided `namemap` mapping of parameter
names.

Both chains share the same underlying data. Any parameters in the chain that are unassigned
will be placed into the `:parameters` section.
"""
function set_section(chains::Chains, namemap)
    # Initialize the name map.
    _namemap = initnamemap(namemap)

    # Make sure all the names are in the new name map.
    newnames = Set(Symbol[])
    names_of_params = names(chains)
    for names in _namemap
        filter!(x -> x ∈ names_of_params, names)
        for name in names
            push!(newnames, name)
        end
    end
    missingnames = setdiff(names_of_params, newnames)

    # Assign everything that is missing to :parameters.
    if !isempty(missingnames)
        @warn "Section mapping does not contain all parameter names, " *
            "$missingnames assigned to :parameters."
        for name in missingnames
            push!(_namemap.parameters, name)
        end
    end

    return Chains(chains.value, chains.logevidence, _namemap, chains.info)
end

_default_sections(c::Chains) = haskey(c.name_map, :parameters) ? :parameters : nothing

function _clean_sections(chains::Chains, sections)
    return filter(sections) do section
        haskey(chains.name_map, Symbol(section))
    end
end
function _clean_sections(chains::Chains, section::Union{String,Symbol})
    return haskey(chains.name_map, Symbol(section)) ? section : ()
end
_clean_sections(::Chains, ::Nothing) = nothing


#################### Concatenation ####################

Base.cat(c::Chains, cs::Chains...; dims = Val(1)) = _cat(dims, c, cs...)
Base.cat(c::T, cs::T...; dims = Val(1)) where T<:Chains = _cat(dims, c, cs...)

Base.vcat(c::Chains, cs::Chains...) = _cat(Val(1), c, cs...)
Base.vcat(c::T, cs::T...) where T<:Chains = _cat(Val(1), c, cs...)

Base.hcat(c::Chains, cs::Chains...) = _cat(Val(2), c, cs...)
Base.hcat(c::T, cs::T...) where T<:Chains = _cat(Val(2), c, cs...)

AbstractMCMC.chainscat(c::Chains, cs::Chains...) = _cat(Val(3), c, cs...)

_cat(dim::Int, cs::Chains...) = _cat(Val(dim), cs...)

function _cat(::Val{1}, c1::Chains, args::Chains...)
    # check inputs
    thin = step(c1)
    all(c -> step(c) == thin, args) || throw(ArgumentError("chain thinning differs"))
    nms = names(c1)
    all(c -> names(c) == nms, args) || throw(ArgumentError("chain names differ"))
    chns = chains(c1)
    all(c -> chains(c) == chns, args) || throw(ArgumentError("sets of chains differ"))

    # concatenate all chains
    data = mapreduce(c -> c.value.data, vcat, args; init = c1.value.data)
    value = AxisArray(data;
                      iter = range(first(c1); length = size(data, 1), step = thin),
                      var = nms,
                      chain = chns)

    return Chains(value, missing, c1.name_map, c1.info)
end

function _cat(::Val{2}, c1::Chains, args::Chains...)
    # check inputs
    rng = range(c1)
    all(c -> range(c) == rng, args) || throw(ArgumentError("chain ranges differ"))
    chns = chains(c1)
    all(c -> chains(c) == chns, args) || throw(ArgumentError("sets of chains differ"))

    # combine names and sections of parameters
    nms = names(c1)
    n = length(nms)
    for c in args
        nms = union(nms, names(c))
        n += length(names(c))
        n == length(nms) || throw(ArgumentError("non-unique parameter names"))
    end

    name_map = mapreduce(c -> c.name_map, merge_union, args; init = c1.name_map)

    # concatenate all chains
    data = mapreduce(c -> c.value.data, hcat, args; init = c1.value.data)
    value = AxisArray(data; iter = rng, var = nms, chain = chns)

    return Chains(value, missing, name_map, c1.info)
end

function _cat(::Val{3}, c1::Chains, args::Chains...)
    # check inputs
    rng = range(c1)
    all(c -> range(c) == rng, args) || throw(ArgumentError("chain ranges differ"))
    nms = names(c1)
    all(c -> names(c) == nms, args) || throw(ArgumentError("chain names differ"))

    # concatenate all chains
    data = mapreduce(c -> c.value.data, (x, y) -> cat(x, y; dims = 3), args;
                     init = c1.value.data)
    value = AxisArray(data; iter = rng, var = nms, chain = 1:size(data, 3))

    return Chains(value, missing, c1.name_map, c1.info)
end

function pool_chain(c::Chains)
    data = c.value.data
    pool_data = reshape(permutedims(data, [1, 3, 2]), :, size(data, 2), 1)
    return Chains(pool_data, names(c), c.name_map; info=c.info)
end

"""
    replacenames(chains::Chains, dict::AbstractDict) 

Replace parameter names by creating a new `Chains` object that shares the same underlying data.

# Examples

```jldoctest; setup = :(using MCMCChains)
julia> chn = Chains(rand(100, 2, 2), ["one", "two"]);

julia> chn2 = replacenames(chn, "one" => "A");

julia> names(chn2)
2-element Vector{Symbol}:
 :A
 :two

julia> chn3 = replacenames(chn2, Dict("A" => "one", "two" => "B"));

julia> names(chn3)
2-element Vector{Symbol}:
 :one
 :Z
```
"""
replacenames(chains::Chains, dict::AbstractDict) = replacenames(chains, pairs(dict)...)
function replacenames(chains::Chains, old_new::Pair...)
    isempty(old_new) && error("you have to specify at least one replacement")

    # Set new parameter names and a new name map.
    names_of_params = copy(names(chains))
    namemap = deepcopy(chains.name_map)
    for (old, new) in old_new
        symold_symnew = Symbol(old) => Symbol(new)

        replace!(names_of_params, symold_symnew)
        for names in namemap
            replace!(names, symold_symnew)
        end
    end

    value = AxisArray(
        chains.value.data;
        iter = range(chains), var = names_of_params, chain = 1:size(chains, 3)
    )

    return Chains(value, chains.logevidence, namemap, chains.info)
end
