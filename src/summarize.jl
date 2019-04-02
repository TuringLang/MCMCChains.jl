using DataFrames: colwise, names!
using StatsBase: mean, std, sem
import StatsBase: sem
import Base.size
import Base.names

struct ChainDataFrame
    df::DataFrame
end

Base.size(c::ChainDataFrame) = size(c.df)
Base.names(c::ChainDataFrame) = names(c.df)
Base.show(io::IO, c::ChainDataFrame) = show(io, c.df)

Base.getindex(c::ChainDataFrame, args...) = getindex(c.df, args...)
Base.getindex(c::ChainDataFrame, s::Union{Symbol, Vector{Symbol}}) = return c.df[s]

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
    sections::Vector{Symbol}=Symbol[:parameters], 
    func_names=[], etype=:bm, args...)
    
    if length(funs) == 0
        return dfsummarystats(chn, sections)
    end

    # Generate a dataframe to work on.
    df = DataFrame(chn, sections)
    
    # If no function names were given, make a new list.
    func_names = length(func_names) == 0 ?
        handle_funs(funs) : Symbol.(func_names)

    # Do all the math, make columns.
    columns = vcat([names(df)], [colwise(f, df) for f in funs])

    # Make a vector of column names.
    colnames = vcat([:parameters], func_names)

    # Build the dataframe.
    df = DataFrame(columns)
    names!(df, colnames)
    return ChainDataFrame(df)
end

function handle_funs(fns)
  tmp =  [string(f) for f in fns] 
  Symbol.([split(tmp[i], ".")[end] for i in 1:length(tmp)])
end

function dfsummarystats(chn::MCMCChains.AbstractChains,
    sections::Vector{Symbol}=Symbol[:parameters]; etype=:bm, args...)
    sem(x) = sqrt(var(x) / length(x))
    df_mcse(x) = mcse(x, etype, args...)
    ess(x) = min((std(x) / df_mcse(x))^2, size(x, 1))
    funs = [mean, std, sem, df_mcse, ess]
    func_names = [:mean, :std, :naive_se, :mcse, :ess]
    return summarize(chn, funs...; sections=sections, func_names=func_names)
end
