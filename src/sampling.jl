"""
    sample(chn::Chains, [wv::AbstractWeights,] n; replace=true, kwargs...)
    sample(rng::Random.AbstractRNG, chn::Chains, [wv::AbstractWeights,] n; kwargs...)

Sample `n` samples from `chn`; see also [`subset`](@ref).
Optionally, the samples can be weighted using `wv`.
Here, `kwargs` defaults to `replace=true` and `ordered=false`.
"""
sample(chn::Chains, n; kwargs...) = sample(Random.GLOBAL_RNG, chn, n; kwargs...)
function sample(chn::Chains, wv::StatsBase.AbstractWeights, n)
    return sample(Random.GLOBAL_RNG, chn, wv, n)
end

function sample(
    rng::Random.AbstractRNG,
    chn::Chains,
    n;
    replace = true,
    ordered = false
)
    return subset(chn, sample(rng, range(chn), n; replace = replace, ordered = ordered))
end

function sample(
    rng::Random.AbstractRNG,
    chn::Chains,
    wv::StatsBase.AbstractWeights,
    n
) 
    return subset(chn, sample(rng, range(chn), wv, n))
end

"""
    subset(chn::Chains, samples)

Return a subset of samples.
"""
function subset(chn::Chains, samples)
    data = AxisArray(chn.value[samples, :, :].data;
                     iter = 1:length(samples), var = names(chn), chain = chains(chn))

    return Chains(data, missing, chn.name_map, chn.info)
end
