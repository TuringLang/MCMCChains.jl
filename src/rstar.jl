"""
    rstar([rng ,]classif::Supervised, chains::Chains; kwargs...)

Compute the ``R^*`` convergence diagnostic of MCMC for the `chains`.

The keyword arguments supported here are the same as those in `rstar` for arrays of samples
and chain indices.
"""
function InferenceDiagnostics.rstar(
    classif::MLJModelInterface.Supervised, chn::Chains; kwargs...
)
    return InferenceDiagnostics.rstar(Random.GLOBAL_RNG, classif, chn; kwargs...)
end

function InferenceDiagnostics.rstar(
    rng::Random.AbstractRNG, classif::MLJModelInterface.Supervised, chn::Chains; kwargs...
)
    nchains = size(chn, 3)
    nchains <= 1 && throw(DimensionMismatch())

    # collect data
    x = Array(chn)
    y = repeat(chains(chn); inner = size(chn,1))

    return InferenceDiagnostics.rstar(rng, classif, x, y; kwargs...)
end
