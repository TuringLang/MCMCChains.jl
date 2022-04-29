"""
    sample([rng,] chn::Chains, [wv::AbstractWeights,] n; replace=true, ordered=false)

Sample `n` samples from the pooled (!) chain `chn`.

The keyword arguments `replace` and `ordered` determine whether sampling is performed with replacement and whether the sample is ordered, respectively.
If specified, sampling probabilities are proportional to weights `wv`.

!!! note
    If `chn` contains multiple chains, they are pooled (i.e., appended) before sampling.
    This ensures that even in this case exactly `n` samples are returned:
    ```jldoctest
    julia> chn = Chains(randn(11, 4, 3));

    julia> size(sample(chn, 7)) == (7, 4)
    true
    ```
"""
function sample(chn::Chains, n::Integer; replace::Bool=true, ordered::Bool=false)
    return sample(Random.GLOBAL_RNG, chn, n; replace=replace, ordered=ordered)
end
function sample(chn::Chains, wv::AbstractWeights, n::Integer; replace::Bool=true, ordered::Bool=false)
    return sample(Random.GLOBAL_RNG, chn, wv, n; replace=replace, ordered=ordered)
end

function sample(
    rng::Random.AbstractRNG,
    chn::Chains,
    n::Integer;
    replace::Bool = true,
    ordered::Bool = false
)
    return _sample(rng, chn, n; replace=replace, ordered=ordered)
end


function sample(
    rng::Random.AbstractRNG,
    chn::Chains,
    wv::AbstractWeights,
    n::Integer;
    replace::Bool = true,
    ordered::Bool = false
) 
    return _sample(rng, chn, wv, n; replace=replace, ordered=ordered)
end

# Internal implementation with generic arguments (possibly including weights) and keyword arguments
# This is not exposed to avoid method ambiguities and to accept incorrect (keyword) arguments 
function _sample(rng, chn, args...; kwargs...)
    data = chn.value.data
    pool_data = reshape(PermutedDimsArray(data, (1, 3, 2)), :, size(data, 2), 1)
    idxs = sample(rng, axes(pool_data, 1), args...; kwargs...)

    samples = AxisArray(pool_data[idxs, :, :]; iter=1:length(idxs), var=names(chn), chain=1:1)

    return Chains(samples, missing, chn.name_map, chn.info)
end

