function MCMCDiagnosticTools.rafterydiag(
    chains::Chains; sections = _default_sections(chains), kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Compute statistics for each chain.
    chains_array = _chains.value.data
    results = map(axes(chains_array, 3)) do chain
        vec_of_namedtuples = map(axes(chains_array, 2)) do param
            return MCMCDiagnosticTools.rafterydiag(
                cskip(@view(_chains.value.data[:, param, chain])); kwargs...
            )
        end
        namedtuple_of_vecs = Tables.columntable(vec_of_namedtuples)
        return namedtuple_of_vecs
    end

    # Create SummaryStats.
    stats = [
        SummaryStats(
            result;
            name = "Raftery and Lewis diagnostic - Chain $i",
            labels = names(_chains),
        ) for (i, result) in enumerate(results)
    ]

    return stats
end
