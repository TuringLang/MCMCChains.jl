function InferenceDiagnostics.gelmandiag(
    chains::Chains{<:Real};
    sections = _default_sections(chains),
    transform = false,
    kwargs...,
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Compute the potential scale reduction factor.
    psi = transform ? link(_chains) : _chains.value.data
    results = InferenceDiagnostics.gelmandiag(psi; kwargs...)

    # Create a data frame with the results.
    df = ChainDataFrame(
        "Gelman, Rubin, and Brooks diagnostic",
        merge((parameters = names(_chains),), results),
    )

    return df
end

function InferenceDiagnostics.gelmandiag_multivariate(
    chains::Chains{<:Real};
    sections = _default_sections(chains),
    transform = true,
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Compute the potential scale reduction factor.
    psi = transform ? link(_chains) : _chains.value.data
    results = InferenceDiagnostics.gelmandiag_multivariate(psi; kwargs...)

    # Create a data frame with the results.
    df = ChainDataFrame(
        "Gelman, Rubin, and Brooks diagnostic",
        (parameters = names(_chains), psrf = results.psrf, psrfci = results.psrfci),
    )

    return df, results.psrfmultivariate
end
