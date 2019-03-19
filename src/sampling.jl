import StatsBase: sample
using Random, KernelDensity

function sample(rng::AbstractRNG, chn::MCMCChains.AbstractChains, n;
    replace=true, ordered=false) 
  indxs = sample(rng, 
    range(chn),
    n,
    replace=replace, ordered=ordered)
  return Chains(chn.value[indxs, :, :],
    names(chn), 
    chn.name_map,
    info = chn.info)
end

function sample(chn::MCMCChains.AbstractChains, n;
    replace=true, ordered=false)
  indxs = sample(Random.GLOBAL_RNG,
    range(chn),
    n,
    replace=replace, ordered=ordered)
  return Chains(chn.value[indxs, :, :],
    names(chn),
    chn.name_map,
    info = chn.info)
end

function sample(rng::AbstractRNG, chn::MCMCChains.AbstractChains,
    wv::AbstractWeights, n) 
  indxs = sample(rng, 
    range(chn), 
    wv, 
    n)
  return Chains(chn.value[indxs, :, :],
    names(chn),
    chn.name_map,
    info = chn.info)
end

function sample(chn::MCMCChains.AbstractChains,
    wv::AbstractWeights, n) 
  indxs = sample(Random.GLOBAL_RNG, 
    range(chn), 
    wv, 
    n)
    return Chains(chn.value[indxs, :, :],
      names(chn),
      chn.name_map,
      info = chn.info)
end
