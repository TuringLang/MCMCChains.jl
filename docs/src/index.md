# MCMCChains

Implementation of Julia types for summarizing MCMC simulations and utility functions for [diagnostics](@ref Diagnostics) and [visualizations](@ref StatsPlots.jl).

## Chains type

```@docs
Chains
```

## Indexing and parameter Names

Chains can be constructed with parameter names.
For example, to create a chains object with

- 500 samples,
- 2 parameters (named `:a` and `:b`)
- 3 chains

use

```@example index
using MCMCChains # hide
using Random # hide
Random.seed!(0) # hide
val = rand(500, 2, 3)
chn = Chains(val, [:a, :b])
```

By default, parameters will be given the name `:param_i`, where `i` is the parameter
number:

```@example index
chn = Chains(rand(100, 2, 2))
```

We can set and get indexes for parameter 2:

```@example index
chn_param2 = chn[1:5,2,:];
```

```@example index
chn[:,2,:] = repeat([4], inner = (100, 1, 2))
chn
```

## Rename Parameters

Parameter names can be changed with the function `replacenames`:

```@docs
replacenames
```

## Sections

Chains parameters are sorted into sections that represent groups of parameters, see 
[`MCMCChains.group`](@ref).
By default, every chain contains a `:parameters` section, to which all unassigned parameters are
assigned to. Chains can be assigned a named map during construction:

```@example index
chn = Chains(rand(100, 4, 2), [:a, :b, :c, :d])
```

The [`MCMCChains.set_section`](@ref) function returns a new `Chains` object:

```@example index
chn2 = set_section(chn, Dict(:internals => [:c, :d]))
```

Note that only `:a` and `:b` are being shown. You can explicity retrieve
an array of the summary statistics and the quantiles of the `:internals` section by
calling `describe(chn; sections = :internals)`, or of all variables with
`describe(chn; sections = nothing)`. Many functions such as [`MCMCChains.summarize`](@ref) or
[`MCMCChains.gelmandiag`](@ref) support the `sections` keyword argument.

## Groups of parameters

You can access the names of all parameters in a `chain` that belong to the group `:name` by using

```@docs
namesingroup
```

## The `get` Function

MCMCChains also provides a [`get`](@ref) function designed to make it easier to access
parameters:

```@example get
using MCMCChains # hide
val = rand(6, 3, 1)
chn = Chains(val, [:a, :b, :c]);

x = get(chn, :a)
```

You can also access the variables via `getproperty`:

```@example get
x.a
```

`get` also accepts vectors of things to retrieve, so you can call 

```@example get
x = get(chn, [:a, :b])
```

## Saving and Loading Chains

Like any Julia object, a `Chains` object can be saved using `Serialization.serialize`
and loaded back by `Serialization.deserialize` as identical as possible.
Note, however, that in general
[this process will not work if the reading and writing are done by different versions of Julia, or an instance of Julia with a different system image](https://docs.julialang.org/en/v1/stdlib/Serialization/#Serialization-1).
You might want to consider [JLSO](https://github.com/invenia/JLSO.jl) for saving metadata
such as the Julia version and the versions of all packages installed as well.

```julia
# Save a chain.
using Serialization
serialize("chain-file.jls", chn)

# Read a chain.
chn2 = deserialize("chain-file.jls")
```

## Exporting Chains

A few utility export functions have been provided to convers `Chains` objects to either an Array or a DataFrame:

```julia
# Several examples of creating an Array object:
Array(chns)
Array(chns[:s])
Array(chns, [:parameters])
Array(chns, [:parameters, :internals])

# By default chains are appended. This can be disabled
# using the append_chains keyword argument:
Array(chns, append_chains=false)

# This will return an `Array{Array, 1}` object containing
# an Array for each chain.

# A final option is:
Array(chns, remove_missing_union=false)

# This will not convert the Array columns from a
# `Union{Missing, Real}` to a `Vector{Real}`.
```

Similarly, for DataFrames:

```julia
DataFrame(chns)
DataFrame(chns[:s])
DataFrame(chns, [:parameters])
DataFrame(chns, [:parameters, :internals])
DataFrame(chns, append_chains=false)
DataFrame(chns, remove_missing_union=false)
```

See also `?DataFrame` and `?Array` for more help.

## Sampling Chains

MCMCChains overloads several `sample` methods as defined in StatsBase:

```julia
# Sampling `n` samples from the chain `a`. Optionally
# weighting the samples using `wv`.
sample([rng], a, [wv::AbstractWeights], n::Integer)

# As above, but supports replacing and ordering.
sample([rng], a, [wv::AbstractWeights], n::Integer; replace=true,
  ordered=false)
```

See also `?sample` for additional help. Alternatively, you can construct
and sample from a kernel density estimator using the KernelDensity package:

```julia
using KernelDensity

# Construct a kernel density estimator
c = kde(Array(chn[:s]))

# Generate 10000 weighted samples from the grid points
chn_weighted_sample = sample(c.x, Weights(c.density), 100000)
```
