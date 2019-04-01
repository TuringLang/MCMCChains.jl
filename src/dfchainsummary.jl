using DataFrames: colwise, names!
using StatsBase: mean, std, sem
import StatsBase: sem
#import MCMCChains: mcse

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

function summarize(chn::Chains, funs...;
    sections::Vector{Symbol}=Symbol[:parameters], func_names=[])
    if length(funs) == 0
        return dfsummarystats(chn, sections)
    end

    # Generate a dataframe to work on.
    df = DataFrame(chn, sections)

    # If no function names were given, make a new list.
    func_names = length(func_names) == 0 ?
        [Symbol(replace(split(string(f), ".", limit=2)[2], "\"" => "")) for f in funs] :
        Symbol.(func_names)

    # Do all the math, make columns.
    columns = vcat([names(df)], [colwise(f, df) for f in funs])

    # Make a vector of column names.
    colnames = vcat([:parameters], func_names)

    # Build the dataframe.
    df = DataFrame(columns)
    names!(df, colnames)
    return ChainDataFrame(df)
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

"""

# DataFrame Chain Summary

Array constructor from an MCMCChains.Chains object. Returns 3 dimensionsal
array or an Array of 2 dimensional Arrays. If only a single parameter is selected for
inclusion, a dimension is dropped in both cases, as is e.g. required by cde(), etc.

### Method
```julia
  dfchainsummary(
    chn::MCMCChains.AbstractChains,
    sections::Vector{Symbol};
    etype::Symbol,
    args
  )
```

### Required arguments
```julia
* `chn` : Chains object to convert to an Array
```

### Optional arguments
```julia
* `sections = Symbol[]` : Sections from the Chains object to be included
* `etype = :bm`  : Default mcse method
```

### Examples
```julia
* `dfchainsummary(chns)` : Complete chain summary
* `dfchainsummary(chns[:par])` : Chain summary of :par only
* `dfchainsummary(chns, [:parameters])`  : Chain summary of :parameter section
* `dfchainsummary(chns, [:parameters, :internals])`  : Chain summary includes multiple sections
```

"""
function dfchainsummary(chn::MCMCChains.AbstractChains,
  sections::Vector{Symbol}=Symbol[]; etype=:bm, args...)

  sem(x) = sqrt(var(x) / length(x))
  df_mcse(x) = mcse(x, etype, args...)

  df = DataFrame(chn, sections)

  # Add summary stats columns
  if size(chn.value, 1) > 200
    sum_df = DataFrame(
      :parameters => names(df),
      :mean => colwise(mean, df),
      :std => colwise(std, df),
      :naive_se => colwise(sem, df),
      :mcse => colwise(df_mcse, df),
      :ess => repeat([size(df, 1)], length(names(df)))
    )
  else
    sum_df = DataFrame(
      :parameters => names(df),
      :mean => colwise(mean, df),
      :std => colwise(std, df),
      :naive_se => colwise(sem, df)
    )
  end

    return ChainDataFrame(sum_df)
end
