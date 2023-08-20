function MCMCDiagnosticTools.gelmandiag(
    chains::Chains{<:Real};
    sections = _default_sections(chains),
    transform = false,
    kwargs...,
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Compute the potential scale reduction factor.
    psi = transform ? link(_chains) : _chains.value.data
    results = MCMCDiagnosticTools.gelmandiag(_permutedims_diagnostics(psi); kwargs...)

    # Create a data frame with the results.
    stats = SummaryStats(
        "Gelman, Rubin, and Brooks diagnostic",
        merge((parameter = names(_chains),), results),
    )

    return stats
end

function MCMCDiagnosticTools.gelmandiag_multivariate(
    chains::Chains{<:Real};
    sections = _default_sections(chains),
    transform = true,
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Compute the potential scale reduction factor.
    psi = transform ? link(_chains) : _chains.value.data
    results = MCMCDiagnosticTools.gelmandiag_multivariate(
        _permutedims_diagnostics(psi);
        kwargs...,
    )

    # Create SummaryStats with the results.
    stats = SummaryStats(
        "Gelman, Rubin, and Brooks diagnostic",
        (parameter = names(_chains), psrf = results.psrf, psrfci = results.psrfci),
    )

    return stats, results.psrfmultivariate
end
