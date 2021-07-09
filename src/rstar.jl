"""
    rstar([rng ,]classif::Supervised, chains::Chains; kwargs...)

Compute the ``R^*`` convergence diagnostic of MCMC for the `chains`.

The keyword arguments supported here are the same as those in `rstar` for arrays of samples
and chain indices.
"""
function MCMCDiagnosticTools.rstar(
    classif::MLJModelInterface.Supervised, chn::Chains; kwargs...
)
    return MCMCDiagnosticTools.rstar(Random.GLOBAL_RNG, classif, chn; kwargs...)
end

function MCMCDiagnosticTools.rstar(
    rng::Random.AbstractRNG, classif::MLJModelInterface.Supervised, chn::Chains; kwargs...
)
    nchains = size(chn, 3)
    nchains <= 1 && throw(DimensionMismatch())

    # collect data
    x = Array(chn)
    y = repeat(chains(chn); inner = size(chn,1))

    return MCMCDiagnosticTools.rstar(rng, classif, x, y; kwargs...)
end
