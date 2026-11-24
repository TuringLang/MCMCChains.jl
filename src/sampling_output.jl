using AbstractMCMC: SamplingOutput

"""
    SamplingChains(output::AbstractMCMC.SamplingOutput)

Wrap structured draws from a sampling run. The `output` field holds the supplied
`SamplingOutput` and its sampler metadata, and can be passed to methods that accept
`SamplingOutput`. Available with AbstractMCMC 5.18 or later.

Index by draw and chain position. Two integer indices return a draw. A selection of draws
or chains wraps a view of the samples, and `range` returns their original iteration numbers.
Saved sampler states still refer to the end of each original run.

```jldoctest
julia> using AbstractMCMC: SamplingOutput

julia> output = SamplingOutput(reshape([(x = i,) for i in 1:3], 3, 1); iterations = 2:2:6);

julia> chain = SamplingChains(output);

julia> chain[2, 1]
(x = 2,)

julia> range(chain[1:2:3, :])
2:4:6
```
"""
struct SamplingChains{O<:SamplingOutput} <: AbstractMCMC.AbstractChains
    output::O
end

Base.size(chain::SamplingChains, args...) = size(chain.output, args...)
Base.firstindex(chain::SamplingChains, dim::Integer) = firstindex(chain.output, dim)
Base.lastindex(chain::SamplingChains, dim::Integer) = lastindex(chain.output, dim)
Base.range(chain::SamplingChains) = chain.output.iterations

Base.getindex(chain::SamplingChains, draw::Integer, run::Integer) = chain.output[draw, run]

function Base.getindex(
    chain::SamplingChains,
    draws::Union{Colon,Integer,AbstractVector{<:Integer}},
    runs::Union{Colon,Integer,AbstractVector{<:Integer}},
)
    return SamplingChains(chain.output[draws, runs])
end
