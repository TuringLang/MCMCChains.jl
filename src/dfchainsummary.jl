using DataFrames: colwise
using StatsBase: mean, std, sem
import StatsBase: sem
import Base.size
import Base.names

struct ChainDataFrame
    df::DataFrame
end

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

names(c::ChainDataFrame) = names(c.df)
size(c::ChainDataFrame, args...) = size(c.df, args...)

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
