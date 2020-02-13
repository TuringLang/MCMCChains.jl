struct ChainDataFrame
    name::String
    df::DataFrame

    function ChainDataFrame(name::String, df::DataFrame; digits=missing)
        if !ismissing(digits)
            if isa(digits, Integer)
                nrow = size(df, 1)
                ncol = size(df, 2)
                for i in 1:nrow
                    for j in 1:ncol
                        if applicable(round, df[i,j])
                            df[i, j] = round(df[i,j], digits=digits)
                        end
                    end
                end
            elseif isa(digits, NamedTuple) || isa(digits, Dict)
                cols = keys(digits)
                nrow = size(df, 1)
                for c in cols
                    if c in names(df)
                        for r in 1:nrow
                            df[r,c] = round(df[r,c], digits=digits[c])
                        end
                    end
                end
            end
        end

        return new(name, df)
    end
end

ChainDataFrame(df::DataFrame) = ChainDataFrame("", df)

Base.size(c::ChainDataFrame) = size(c.df)
Base.names(c::ChainDataFrame) = names(c.df)
function Base.show(io::IO, c::ChainDataFrame)
    println(io, c.name)
    show(io, c.df, summary = false, allrows=true, allcols=true)
end

Base.getindex(c::ChainDataFrame, args...) = getindex(c.df, args...)
Base.getindex(c::ChainDataFrame, s::Union{Symbol, Vector{Symbol}}) = c.df[:, s]
Base.isequal(cs1::Vector{ChainDataFrame}, cs2::Vector{ChainDataFrame}) = isequal.(cs1, cs2)
Base.isequal(c1::ChainDataFrame, c2::ChainDataFrame) = isequal(c1, c2)

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
    return c.df[map(x -> x in s1, c.df[:, :parameters]), s2]
end

function Base.getindex(c::ChainDataFrame,
        s1::Symbol,
        s2::Union{Symbol, Vector{Symbol}})
    return c.df[c.df[:, :parameters] .== s1, s2]
end


Base.lastindex(c::ChainDataFrame, i::Integer) = lastindex(c.df, i)

Base.convert(::Type{Array{ChainDataFrame,1}}, cs::Array{ChainDataFrame,1}) = cs
function Base.convert(::Type{T}, cs::Array{ChainDataFrame,1}) where T<:Array
    arrs = [convert(T, cs[j].df[:, 2:end]) for j in 1:length(cs)]
    return cat(arrs..., dims = 3)
end
"""

# Summarize a Chains object formatted as a DataFrame

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
        name::String="",
        additional_df=nothing,
        digits=missing)
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
        vcat([names(df)],
             [[f(col) for col = eachcol(df, false)] for f in funs])
    else
        [vcat([names(df[1])],
              [[f(col) for col = eachcol(i, false)] for f in funs]) for i in df]
    end

    # Make a vector of column names.
    colnames = vcat([:parameters], func_names)

    # Build the dataframes.
    ret_df = if append_chains
        DataFrame(columns, colnames)
    else
        [DataFrame(i, colnames) for i in columns]
    end

    if additional_df != nothing
        if append_chains
            ret_df = join(ret_df, additional_df, on=:parameters)
        else
            ret_df = [join(r, additional_df, on=:parameters) for r in ret_df]
        end
    end

    ret_df = if append_chains
        ChainDataFrame(name, ret_df, digits=digits)
    else
        [ChainDataFrame(name, r, digits=digits) for r in ret_df]
    end

    return ret_df
end

function handle_funs(fns)
    tmp =  [string(f) for f in fns]
    Symbol.([split(tmp[i], ".")[end] for i in 1:length(tmp)])
end
