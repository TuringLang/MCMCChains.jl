"""
    discretediag(chains::Chains{<:Real}; sections, kwargs...)

Discrete diagnostic where `method` can be
`[:weiss, :hangartner, :DARBOOT, MCBOOT, :billinsgley, :billingsleyBOOT]`.
"""
function MCMCDiagnosticTools.discretediag(
    chains::Chains{<:Real};
    sections = _default_sections(chains),
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Compute statistics.
    between_chain_vals, within_chain_vals = MCMCDiagnosticTools.discretediag(
        _permutedims_diagnostics(_chains.value.data); kwargs...
    )

    # Create SummaryStats
    param_names = names(_chains)
    between_chain_stats = SummaryStats(
        between_chain_vals;
        name = "Chisq diagnostic - Between chains",
        labels = param_names,
    )
    within_chain_stats = map(1:size(_chains, 3)) do i
        vals = map(val -> val[:, i], within_chain_vals)
        return SummaryStats(
            vals;
            name = "Chisq diagnostic - Chain $i",
            labels = param_names,
        )
    end
    stats = vcat([between_chain_stats], within_chain_stats)

    return stats
end
