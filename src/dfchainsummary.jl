using DataFrames: colwise
using StatsBase: mean, std, sem
import StatsBase: sem
#import MCMCChains: mcse

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
    
    return sum_df
end
