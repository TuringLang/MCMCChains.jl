sample(chn::Chains, n; kwargs...) = sample(Random.GLOBAL_RNG, chn, n; kwargs...)
function sample(chn::Chains, wv::AbstractWeights, n)
    return sample(Random.GLOBAL_RNG, chn, wv, n)
end

function sample(
    rng::Random.AbstractRNG,
    chn::Chains,
    n;
    replace = true,
    ordered = false
) 
    indxs = sample(rng, range(chn), n; replace = replace, ordered = ordered)

    return Chains(chn.value[indxs, :, :], names(chn), chn.name_map; info = chn.info)
end

function sample(
    rng::Random.AbstractRNG,
    chn::Chains,
    wv::AbstractWeights,
    n
) 
    indxs = sample(rng, range(chn), wv, n)

    return Chains(chn.value[indxs, :, :], names(chn), chn.name_map, info = chn.info)
end