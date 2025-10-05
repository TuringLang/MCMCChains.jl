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
        results;
        name = "Gelman, Rubin, and Brooks diagnostic",
        labels = names(_chains),
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
        (psrf = results.psrf, psrfci = results.psrfci);
        name = "Gelman, Rubin, and Brooks diagnostic",
        labels = names(_chains),
    )

    return stats, results.psrfmultivariate
end
