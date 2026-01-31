# Serialization

MCMCChains supports saving and loading `Chains` objects in multiple formats:

- **JSON** - Full round-trip with metadata preservation
- **CSV** - Simple tabular data export
- **StanCSV** - Stan-compatible format for interoperability with Stan ecosystem tools

## JSON

The JSON format preserves all chain data including metadata, parameter sections, and iteration information.

```julia
using MCMCChains, JSON

chn = Chains(rand(100, 2, 2), [:a, :b])

# Serialize to string
json_str = JSON.json(chn)

# Save to file
JSON.json("chains.json", chn)

# Load from file
chn_loaded = JSON.parsefile("chains.json", Chains)

# Parse from string
chn_from_str = JSON.parse(json_str, Chains)
```

### What Gets Saved

- Parameter values across all iterations and chains
- Parameter names and section groupings (`:parameters`, `:internals`)
- Iteration range and thinning information
- Log evidence (if available)
- Chain metadata (`info` field)

## CSV

MCMCChains implements the Tables.jl interface, so you can use CSV.jl directly for simple data export.

```julia
using MCMCChains, CSV

chn = Chains(rand(100, 2, 2), [:a, :b])

# Save (includes iteration and chain columns)
CSV.write("chains.csv", chn)

# Load
chn_loaded = Chains(CSV.File("chains.csv"))
```

!!! note
    Simple CSV export flattens all chains into rows and adds `iteration` and `chain` columns.
    This is useful for quick data export but doesn't preserve all metadata.

## StanCSV Format

For interoperability with Stan ecosystem tools (CmdStan, ArviZ, etc.), MCMCChains provides dedicated StanCSV functions that handle the Stan-specific format including comment headers, parameter name conventions, and column ordering.

### Stan CSV Format Features

- **Adaptation comments**: Step size and mass matrix info
- **Timing comments**: Elapsed time for warmup/sampling  
- **Parameter naming**: `theta[1,2]` → `theta.1.2` (dot notation)
- **Column ordering**: Sampler parameters (`lp__`, `accept_stat__`, etc.) appear first
- **Internals detection**: Parameters ending in `__` auto-classified as internals

### Writing StanCSV

```julia
using MCMCChains, CSV

chn = Chains(rand(100, 3, 2), [:mu, :sigma, :lp__], 
             Dict(:internals => [:lp__]))

# Write single chain
write_stancsv("chain_1.csv", chn; chain_id=1)

# Write all chains to separate files (chain_1.csv, chain_2.csv, ...)
write_stancsv("chain.csv", chn, Val(:all))

# Write without comments (plain CSV with Stan column naming)
write_stancsv("chain.csv", chn; include_adaptation=false, include_timing=false)
```

### Reading StanCSV

```julia
using MCMCChains, CSV

# Read single file (handles CmdStan output directly)
chn = read_stancsv("output.csv")

# Read multiple chain files
chn = read_stancsv(["chain_1.csv", "chain_2.csv", "chain_3.csv"])
```

The reader:
- Skips all comment lines (starting with `#`)
- Converts Stan parameter names back to Julia format (`theta.1.2` → `theta[1,2]`)
- Extracts metadata from comments (model name, seed, warmup, step size)
- Classifies parameters ending in `__` as internals

### Example Output

```
# Adaptation terminated
# Step size = 0.8
# Diagonal elements of inverse mass matrix:
# 1.0
lp__,mu,sigma
-6.74827,0.247195,1.5
-6.74827,0.280619,1.3
...
#
#  Elapsed Time: 0.01 seconds (Warm-up)
#                0.02 seconds (Sampling)
#                0.03 seconds (Total)
#
```

## Format Comparison

| Feature | JSON | CSV | StanCSV |
|---------|------|-----|---------|
| Full metadata | ✓ | ✗ | Partial |
| Section info | ✓ | ✗ | Auto-detected |
| Iteration range | ✓ | ✗ | ✗ |
| Stan compatible | ✗ | ✗ | ✓ |
| Human readable | ✓ | ✓ | ✓ |
| Multi-chain in one file | ✓ | ✓ | ✗ |
| Requires package | JSON.jl | CSV.jl | CSV.jl |

