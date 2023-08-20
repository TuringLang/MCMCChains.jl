"""
    mcse(chains::Chains; duration=compute_duration, kwargs...)

Estimate the Monte Carlo standard error.
"""
function MCMCDiagnosticTools.mcse(
    chains::Chains;
    sections = _default_sections(chains), kwargs...
)
    # Subset the chain
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Estimate the effective sample size
    mcse = MCMCDiagnosticTools.mcse(
        _permutedims_diagnostics(_chains.value.data);
        kwargs...,
    )

    nt = merge((parameter = names(_chains),), (; mcse))

    return SummaryStats("MCSE", nt)
end
