# Serialization

MCMCChains provides functionality to save and load `Chains` objects to and from JSON format.
This is useful for:

- Sharing MCMC results with collaborators
- Archiving simulation results for reproducibility
- Exporting data for analysis in other tools
- Creating portable analysis pipelines

!!! note
    JSON serialization requires the `JSON` package. Make sure to run `using JSON` before using these functions.

## Quick Start

```julia
using MCMCChains
using JSON

# Create or sample a chain
chn = Chains(rand(1000, 2, 3), [:param1, :param2])

# Save to a JSON file
JSON.json("results.json", chn)

# Load back from JSON
chn_loaded = JSON.parsefile("results.json", Chains)
```

## Saving Chains to JSON

### Examples

Save to a specific file:

```julia
using MCMCChains, JSON

chn = Chains(rand(500, 2, 2), [:a, :b])
JSON.json("my_chain.json", chn)
```

Get JSON as a string instead of writing to file:

```julia
json_string = JSON.json(chn)
```

## Loading Chains from JSON

### Example

```julia
using MCMCChains, JSON

# Load a previously saved chain
chn = JSON.parsefile("my_chain.json", Chains)

# Use the chain normally
describe(chn)
```

## What Gets Serialized

The JSON serialization captures:

- **Numerical data**: All parameter values across iterations and chains
- **Dimensions**: Number of iterations, parameters, and chains
- **Parameter names**: Names of all variables
- **Sections**: Parameter groupings (`:parameters`, `:internals`, etc.)
- **Metadata**: Chain info including model name, sampler details, timing information

## Data Format

The JSON structure is designed to be both human-readable and efficient:

```json
{
  "size": [1000, 5, 2],
  "value_flat": [...],
  "iterations": [1, 2, 3, ...],
  "parameters": ["a", "b", "c", "d", "e"],
  "chains": [1, 2],
  "logevidence": null,
  "name_map": {
    "parameters": ["a", "b"],
    "internals": ["c", "d", "e"]
  },
  "info": {
    "model_name": "mymodel",
    "start_time": 1234567890.0,
    "sampler": "NUTS"
  }
}
```

!!! info "Type Safety"
    Complex Julia types (e.g., `Symbol`, `Missing`, `NamedTuple`) are automatically converted to JSON-compatible types during serialization and restored during deserialization.

## Integration with Turing.jl

```julia
using Turing, MCMCChains, JSON

@model function gdemo(x)
    s ~ InverseGamma(2, 3)
    m ~ Normal(0, sqrt(s))
    for i in eachindex(x)
        x[i] ~ Normal(m, sqrt(s))
    end
end

# Sample from model
chn = sample(gdemo([1.5, 2.0]), NUTS(), 1000)

# Add model name for better organization
chn = setinfo(chn, merge(chn.info, (model_name = "gdemo",)))

# Save - will create "gdemo.json"
JSON.json("gdemo.json", chn)

# Load and analyze
chn_loaded = JSON.parsefile("gdemo.json", Chains)
describe(chn_loaded)
```
