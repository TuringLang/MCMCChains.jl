struct ChainDataFrame{NT<:NamedTuple}
    name::String
    nt::NT
    nrows::Int
    ncols::Int

    function ChainDataFrame(name::String, nt::NamedTuple)
        lengths = length(first(nt))
        all(x -> length(x) == lengths, nt) || error("Lengths must be equal.")

        return new{typeof(nt)}(name, nt, lengths, length(nt))
    end
end

ChainDataFrame(nt::NamedTuple) = ChainDataFrame("", nt)

Base.size(c::ChainDataFrame) = (c.nrows, c.ncols)
Base.names(c::ChainDataFrame) = collect(keys(c.nt))

# Display

function Base.show(io::IO, df::ChainDataFrame)
    print(io, df.name, " (", df.nrows, " x ", df.ncols, ")")
end

function Base.show(io::IO, ::MIME"text/plain", df::ChainDataFrame)
    digits = get(io, :digits, 4)
    formatter = PrettyTables.ft_printf("%.$(digits)f")

    println(io, df.name)
    PrettyTables.pretty_table(io, df.nt; formatter = formatter, tf = PrettyTables.borderless)
end

Base.isequal(c1::ChainDataFrame, c2::ChainDataFrame) = isequal(c1, c2)

# Index functions
function Base.getindex(c::ChainDataFrame, s::Union{Colon, Integer, UnitRange}, g::Union{Colon, Integer, UnitRange})
    convert(Array, getindex(c, c.nt[:parameters][s], collect(keys(c.nt))[g]))
end

function Base.getindex(c::ChainDataFrame, s::Union{Symbol, Vector{Symbol}}, g::Colon)
    getindex(c, c.nt[:parameters], collect(keys(c.nt)))
end

function Base.getindex(c::ChainDataFrame, s::Union{Symbol, Vector{Symbol}})
    getindex(c, s, collect(keys(c.nt)))
end

function Base.getindex(c::ChainDataFrame, s::Union{Colon, Integer, UnitRange}, ks::Union{Symbol, Vector{Symbol}})
    getindex(c, c.nt[:parameters][s], ks)
end

function Base.getindex(
    c::ChainDataFrame, 
    s::Union{String, Vector{String}, Symbol, Vector{Symbol}}, 
    ks::Union{Symbol, Vector{Symbol}}
)
    s = s isa AbstractArray ? s : [s]
    ks = ks isa AbstractArray ? ks : [ks]
    ind = indexin(Symbol.(s), Symbol.(c.nt[:parameters]))

    not_found = map(x -> x === nothing, ind)

    any(not_found) && error("Cannot find parameters $(s[not_found]) in chain")

    # If there are multiple columns, return a new CDF.
    if length(ks) > 1
        if !(:parameters in ks)
            ks = vcat(:parameters, ks)
        end
        nt = NamedTuple{tuple(ks...)}(tuple([c.nt[k][ind] for k in ks]...))
        return ChainDataFrame(c.name, nt)
    else
        # Otherwise, return a vector if there's multiple parameters
        # or just a scalar if there's one parameter.
        if length(s) == 1
            return c.nt[ks[1]][ind][1]
        else
            return c.nt[ks[1]][ind]
        end
    end
end

function Base.lastindex(c::ChainDataFrame, i::Integer)
    if i == 1
        return c.nrows
    elseif i ==2
        return c.ncols
    else
        error("No such dimension")
    end
end

function Base.convert(::Type{Array}, c::C) where C<:ChainDataFrame
    T = promote_eltype_namedtuple_tail(c.nt)
    arr = Array{T, 2}(undef, c.nrows, c.ncols - 1)
    
    for (i, k) in enumerate(Iterators.drop(keys(c.nt), 1))
        arr[:, i] = c.nt[k]
    end

    return arr
end

Base.convert(::Type{Array{ChainDataFrame,1}}, cs::Array{ChainDataFrame,1}) = cs
function Base.convert(::Type{Array}, cs::Array{C,1}) where C<:ChainDataFrame
    return mapreduce((x, y) -> cat(x, y; dims = Val(3)), cs) do c
        reshape(convert(Array, c), Val(3))
    end
end

"""

# Summarize a Chains object formatted as a ChainsDataFrame

Summarize method for a Chains object.

### Method
```julia
  summarize(
    chn::Chains,
    funs...;
    sections::Vector{Symbol}=[:parameters],
    func_names=[],
    etype=:bm
  )
```

### Required arguments
```julia
* `chn` : Chains object to convert to a ChainsDataFrame-formatted summary
```

### Optional arguments
```julia
* `funs...` : zero or more vector functions, e.g. mean, std, etc.
* `sections = [:parameters]` : Sections from the Chains object to be included
* `etype = :bm` : Default for df_mcse
```

### Examples
```julia
* `summarize(chns)` : Complete chain summary
* `summarize(chns[[:parm1, :parm2]])` : Chain summary of selected parameters
* `summarize(chns, sections=[:parameters])`  : Chain summary of :parameters section
* `summarize(chns, sections=[:parameters, :internals])` : Chain summary for multiple sections
```

"""
function summarize(chains::Chains, funs...;
        sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
        func_names::AbstractVector{Symbol} = Symbol[],
        append_chains::Bool=true,
        showall::Bool=false,
        name::String="",
        additional_df=nothing,
        sorted::Bool=false)
    # Check that we actually have :parameters.
    showall = showall || !in(:parameters, keys(chains.name_map))
    
    # If we weren't given any functions, fall back on summary stats.
    if length(funs) == 0
        return summarystats(chains,
            sections=sections,
            showall=showall)
    end
    
    # Generate a chain to work on.
    chn = Chains(chains, sections, sorted=sorted)

    # Obtain names of parameters.
    names_of_params = names(chn)

    # If no function names were given, make a new list.
    fnames = isempty(func_names) ? collect(nameof.(funs)) : func_names

    # Obtain the additional named tuple.
    additional_nt = additional_df === nothing ? NamedTuple() : additional_df.nt

    if append_chains
        # Evaluate the functions.
        data = to_matrix(chn)
        fvals = [[f(data[:, i]) for i in axes(data, 2)] for f in funs]

        # Build the ChainDataFrame.
        nt = merge((; parameters = names_of_params, zip(fnames, fvals)...), additional_nt)
        df = ChainDataFrame(name, nt)

        return df
    else
        # Evaluate the functions.
        data = to_vector_of_matrices(chn)
        vector_of_fvals = [[[f(x[:, i]) for i in axes(x, 2)] for f in funs] for x in data]

        # Build the ChainDataFrames.
        vector_of_nt = [
            merge((; parameters = names_of_params, zip(fnames, fvals)...), additional_nt)
            for fvals in vector_of_fvals
        ]
        vector_of_df = [
            ChainDataFrame(name * " (Chain $i)", nt)
            for (i, nt) in enumerate(vector_of_nt)
        ]

        return vector_of_df
    end
end