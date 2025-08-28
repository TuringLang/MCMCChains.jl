# Chains

The methods listed below are defined in `src/chains.jl`.

```@autodocs
Modules = [MCMCChains]
Pages = ["chains.jl"]
```

## Internals Fields

When using MCMC samplers, especially Hamiltonian Monte Carlo (HMC) and No-U-Turn Sampler (NUTS), the `Chains` object contains additional diagnostic information in the `:internals` section. These fields provide valuable insights into the sampler's behavior and can help diagnose sampling issues.

### Common Internal Fields

The following internal fields are commonly available when using HMC/NUTS samplers:

- **`acceptance_rate`**: The acceptance rate of the sampler. For NUTS, this is the average acceptance probability across the trajectory.

- **`hamiltonian_energy`**: The Hamiltonian energy at each sample point. This represents the total energy (kinetic + potential) of the system and is used in energy plots for diagnosing sampler performance.

- **`hamiltonian_energy_error`**: The error in Hamiltonian energy due to numerical integration. Large values may indicate issues with the step size or numerical precision.

- **`is_accept`**: Whether each sample was accepted. For NUTS, this is typically always `true` since NUTS uses a different acceptance mechanism.

- **`log_density`**: The log probability density of the target distribution at each sample point.

- **`lp`**: An alias for log probability density (may be identical to `log_density`).

- **`max_hamiltonian_energy_error`**: The maximum energy error encountered during the trajectory.

- **`n_steps`**: The number of leapfrog steps taken for each sample. In NUTS, this varies as the algorithm adapts the trajectory length.

- **`nom_step_size`**: The nominal step size used by the integrator.

- **`numerical_error`**: Indicates whether numerical errors were encountered during sampling.

- **`step_size`**: The actual step size used by the integrator for each sample.

- **`tree_depth`**: The depth of the binary tree built by NUTS. Higher values indicate longer trajectories.

### Using Internal Fields

You can access internal fields using:

```julia
# Get all internal field names
names(chain, :internals)

# Extract specific internal field
chain[:, :step_size, :]

# Get summary statistics for internals
describe(chain, sections=:internals)
```

### Diagnostic Applications

These internal fields are particularly useful for:

- **Energy plots**: Use `hamiltonian_energy` and `hamiltonian_energy_error` with `energyplot()` to diagnose sampler efficiency
- **Step size analysis**: Monitor `step_size` and `n_steps` to understand sampler adaptation
- **Convergence diagnostics**: Check `acceptance_rate` and `tree_depth` for sampling quality
- **Troubleshooting**: Use `numerical_error` and energy errors to identify sampling problems

See the plotting documentation for examples of using these fields in diagnostic plots.
