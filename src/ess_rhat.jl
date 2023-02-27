"""
    ess(chains::Chains; duration=compute_duration, kwargs...)

Estimate the effective sample size.

ESS per second options include `duration=MCMCChains.compute_duration` (the default)
and `duration=MCMCChains.wall_duration`.
"""
function MCMCDiagnosticTools.ess(
    chains::Chains;
    sections = _default_sections(chains), duration = compute_duration, kwargs...
)
    # Subset the chain
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Estimate the effective sample size
    ess = MCMCDiagnosticTools.ess(
        _permutedims_diagnostics(_chains.value.data);
        kwargs...,
    )

    # Calculate ESS/minute if available
    dur = duration(chains)

    # Convert to NamedTuple
    nt = if dur === missing
        merge((parameters = names(_chains),), (; ess))
    else
        ess_per_sec = ess/dur
        merge((parameters = names(_chains),), (; ess, ess_per_sec))
    end

    return ChainDataFrame("ESS", nt)
end

"""
    rhat(chains::Chains; kwargs...)

Estimate the ``\\widehat{R}`` diagnostic.
"""
function MCMCDiagnosticTools.rhat(
    chains::Chains;
    sections = _default_sections(chains), kwargs...
)
    # Subset the chain
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Estimate the rhat
    rhat = MCMCDiagnosticTools.rhat(
        _permutedims_diagnostics(_chains.value.data);
        kwargs...,
    )

    # Convert to NamedTuple
    nt = merge((parameters = names(_chains),), (; rhat))

    return ChainDataFrame("R-hat", nt)
end

"""
    ess_rhat(chains::Chains; duration=compute_duration, kwargs...)

Estimate the effective sample size and the ``\\widehat{R}`` diagnostic

ESS per second options include `duration=MCMCChains.compute_duration` (the default)
and `duration=MCMCChains.wall_duration`.
"""
function MCMCDiagnosticTools.ess_rhat(
    chains::Chains;
    sections = _default_sections(chains), duration = compute_duration, kwargs...
)
    # Subset the chain
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Estimate the effective sample size and rhat
    ess_rhat = MCMCDiagnosticTools.ess_rhat(
        _permutedims_diagnostics(_chains.value.data);
        kwargs...,
    )

    # Calculate ESS/minute if available
    dur = duration(chains)

    # Convert to NamedTuple
    nt = if dur === missing
        merge((parameters = names(_chains),), ess_rhat)
    else
        ess_per_sec = ess_rhat.ess/dur
        merge((parameters = names(_chains),), ess_rhat, (; ess_per_sec))
    end

    return ChainDataFrame("ESS/R-hat", nt)
end
