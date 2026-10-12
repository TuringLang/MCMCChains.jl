@recipe function f(p::EnergyPlot; kind = :density)
    chains = p.args[1]
    _apply_chrome!(plotattributes)

    if kind ∉ (:density, :histogram)
        error("`kind` must be one of `:density` or `:histogram`")
    end

    internal_names = names(chains, :internals)
    required_params = [:hamiltonian_energy, :hamiltonian_energy_error]
    for param in required_params
        if param ∉ internal_names
            error(
                "`$param` not found in chain's internal parameters. Energy plots are only available for HMC/NUTS samplers.",
            )
        end
    end

    pooled = pool_chain(chains)
    energy = vec(pooled[:, :hamiltonian_energy, :])
    energy_error = vec(pooled[:, :hamiltonian_energy_error, :])

    mean_energy = mean(energy)
    std_energy = std(energy)
    centered_energy = (energy .- mean_energy) ./ std_energy
    scaled_energy_error = energy_error ./ std_energy

    title := "Energy Plot"
    xaxis := "Standardized Energy"
    yaxis := "Density"
    legend := :topright

    @series begin
        seriestype := kind
        label := "Marginal Energy"
        fillrange --> 0
        fillalpha --> FILL_ALPHA
        normalize --> true
        bins --> 50
        centered_energy
    end

    @series begin
        seriestype := kind
        label := "Energy Transition"
        fillrange --> 0
        fillalpha --> FILL_ALPHA
        normalize --> true
        bins --> 50
        scaled_energy_error
    end
end
