# Getting started

## Chains type

```@docs
Chains
```

## Indexing and parameter Names

Chains can be constructed with parameter names.
For example, to create a chains object with

- 500 samples,
- 3 parameters (named `a`, `b`, and `c`)
- 3 chains

use

```@example index
using MCMCChains # hide
using Random # hide
Random.seed!(0) # hide
val = rand(500, 3, 3)
chn = Chains(val, [:a, :b, :c])
```

Select parameters by name with a vector to keep a `Chains` object for plotting or summaries:

```@example index
selected = chn[[:a, :c]]
size(selected)
```

A single name returns a matrix with iterations in rows and chains in columns.
A one-element vector keeps all three dimensions:

```@example index
size(chn[:a]), size(chn[[:a]])
```

Selected parameters follow their order in the original chain.

By default, parameters will be given the name `param_i`, where `i` is the parameter
number:

```@example index
chn = Chains(rand(100, 2, 2))
```

We can set and get indexes for parameter 2:

```@example index
chn_param2 = chn[1:5,2,:];
```

```@example index
chn[:,2,:] = fill(4, 100, 1, 2)
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
By default, every chain contains a `parameters` section, to which all unassigned parameters are
assigned to. Chains can be assigned a named map during construction:

```@example index
chn = Chains(rand(100, 4, 2), [:a, :b, :c, :d])
```

The [`MCMCChains.set_section`](@ref) function returns a new `Chains` object:

```@example index
chn2 = set_section(chn, Dict(:internals => [:c, :d]))
```

Note that only `a` and `b` are being shown. You can explicity retrieve
an array of the summary statistics and the quantiles of the `:internals` section by
calling `describe(chn; sections = :internals)`, or of all variables with
`describe(chn; sections = nothing)`. Many functions such as [`MCMCChains.summarize`](@ref)
support the `sections` keyword argument.

## Groups of parameters

You can access the names of all parameters in a `chain` that belong to the group `name` by using

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

Use `flatten=true` to keep each parameter's full draws under its own name when reading a section:
Only names in the selected section appear in the result.
For a single draw from one chain, the value is a scalar.

```@example get
get(chn; section=:parameters, flatten=true).a == chn[:a]
```

## Add derived parameters

Use `hcat` to store quantities computed from each draw:

```@example derived
using MCMCChains
chn = Chains(reshape(1:24, 4, 3, 2), [:a, :b, :c]; start=101, thin=2)
total = chn[:a] .+ chn[:b]
extended = hcat(chn, :total => total)
```

For a vector quantity, pass an iteration-by-chain matrix of vectors or an iteration-by-parameter-by-chain array.
The components are named `quantity[1]`, `quantity[2]`, and so on.
An iteration-by-chain matrix of named tuples, such as model return values, adds each named field as a parameter.
The fields may appear in different orders, but every draw must have the same field names.
Names must be distinct from those already in the chain.
Unlabeled arrays use the chain's iteration and chain order.
An `AxisArray` must have `:iter` and `:chain` axes whose labels match the chain's iteration labels and chain IDs.

```@example derived
generated = map(chn[:a], chn[:b]) do a, b
    (sum = a + b, prediction = [a, a + b])
end
with_generated = hcat(chn, generated)
namesingroup(with_generated, :prediction)
```

## Saving and Loading Chains

Like any Julia object, a `Chains` object can be saved using `Serialization.serialize`
and loaded back by `Serialization.deserialize` as identical as possible.
Note, however, that in general
[this process will not work if the reading and writing are done by different versions of Julia, or an instance of Julia with a different system image](https://docs.julialang.org/en/v1/stdlib/Serialization/#Serialization-1).
You might want to consider [JLSO](https://github.com/invenia/JLSO.jl) for saving metadata
such as the Julia version and the versions of all packages installed as well.

```julia
using Serialization

serialize("chain-file.jls", chn)
chn2 = deserialize("chain-file.jls")
```

The [MCMCChainsStorage.jl](https://github.com/farr/MCMCChainsStorage.jl) package also provides the ability to serialize/deserialize a chain to an HDF5 file across different versions of Julia and/or different system images.

## Exporting Chains

A few utility export functions have been provided to convert `Chains` objects to either an Array or a DataFrame:

```@example exporting
using MCMCChains # hide

chn = Chains(rand(3, 2, 2), [:a, :b])

Array(chn)
```

```@example exporting
Array(chn, [:parameters])
```

By default chains are appended. This can be disabled by using the `append_chains` keyword 
argument:

```@example exporting
A = Array(chn, append_chains=false)
```

which will return a matrix for each chain. For example, for the second chain:

```@example exporting
A[2]
```

Similarly, for DataFrames:

```@example exporting
using DataFrames

DataFrame(chn)
```

See also `?DataFrame` and `?Array` for more help.

## Sampling Chains

MCMCChains overloads several `sample` methods as defined in StatsBase:

```@docs
MCMCChains.sample(::Chains, ::Integer)
```

See `?sample` for additional help on sampling.
Alternatively, you can construct and sample from a kernel density estimator using
[KernelDensity.jl](https://github.com/JuliaStats/KernelDensity.jl),
see [test/sampling_tests.jl](https://github.com/TuringLang/MCMCChains.jl/blob/main/test/sampling_tests.jl).
