using DataFrames: colwise, names!
using StatsBase: mean, std, sem
import StatsBase: sem
import Base.size
import Base.names

struct ChainDataFrame
    name::String
    df::DataFrame
end

ChainDataFrame(df::DataFrame) = ChainDataFrame("", df)

Base.size(c::ChainDataFrame) = size(c.df)
Base.names(c::ChainDataFrame) = names(c.df)

function Base.show(io::IO, c::ChainDataFrame)
    println(io, c.name)
    show(io, c.df, summary = false, allrows=true)
end

Base.getindex(c::ChainDataFrame, args...) = getindex(c.df, args...)
Base.getindex(c::ChainDataFrame, s::Union{Symbol, Vector{Symbol}}) = c.df[s]

function Base.show(io::IO, cs::Vector{ChainDataFrame})
    println(io, summary(cs))
    for i in cs
        println(io)
        println(io, i.name)
        show(io, i.df, allrows=true, summary=false)
        println(io)
    end
end

# Allows overriding of `display`
function Base.show(io::IO, ::MIME"text/plain", cs::Vector{ChainDataFrame})
    show(io, cs)
end

function Base.getindex(c::ChainDataFrame, s::Union{Symbol, Vector{Symbol}}, m)
    s = s isa AbstractArray ? s : [s]
    s_ind = indexin(s, c.df[:, :parameters])
    return c.df[s_ind, m]
end

function Base.getindex(c::ChainDataFrame,
        s1::Vector{Symbol},
        s2::Union{Symbol, Vector{Symbol}})
    return c.df[map(x -> x in s1, c.df.parameters), s2]
end

function Base.getindex(c::ChainDataFrame,
        s1::Symbol,
        s2::Union{Symbol, Vector{Symbol}})
    return c.df[c.df.parameters .== s1, s2]
end


Base.lastindex(c::ChainDataFrame, i::Integer) = lastindex(c.df, i)

function Base.convert(::Type{T}, cs::Array{ChainDataFrame,1}) where T<:Array
    arrs = [convert(T, cs[j].df[:, 2:end]) for j in 1:length(cs)]
    return cat(arrs..., dims = 3)
end
"""

# Summarize a Chains object formatted as a DataFrame

Summarize method for an MCMCChains.Chains object.

### Method
```julia
  summarize(
    chn::MCMCChains.AbstractChains,
    funs...;
    sections::Vector{Symbol}=[:parameters],
    func_names=[],
    etype=:bm
  )
```

### Required arguments
```julia
* `chn` : Chains object to convert to a DataFrame-formatted summary
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
function summarize(chn::Chains, funs...;
        sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
        func_names=[],
        append_chains::Bool=true,
        showall::Bool=false,
        name::String="")
    # Check that we actually have :parameters.
    showall = showall ? true : !in(:parameters, keys(chn.name_map))

    # If we weren't given any functions, fall back on summary stats.
    if length(funs) == 0
        return summarystats(chn,
            sections=sections,
            showall=showall)
    end

    # Generate a dataframe to work on.
    df = DataFrame(chn, sections, showall=showall, append_chains=append_chains)

    # If no function names were given, make a new list.
    func_names = length(func_names) == 0 ?
        handle_funs(funs) : Symbol.(func_names)

    # Do all the math, make columns.
    columns = if append_chains
        vcat([names(df)], [colwise(f, df) for f in funs])
    else
        [vcat([names(df[1])], [colwise(f, i) for f in funs]) for i in df]
    end

    # Make a vector of column names.
    colnames = vcat([:parameters], func_names)

    # Build the dataframe.
    ret_df = if append_chains
        ChainDataFrame(name, DataFrame(columns, colnames))
    else
        [ChainDataFrame(name, DataFrame(i, colnames)) for i in columns]
    end

    return ret_df
end

function handle_funs(fns)
    tmp =  [string(f) for f in fns]
    Symbol.([split(tmp[i], ".")[end] for i in 1:length(tmp)])
end
