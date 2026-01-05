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
write_json(chn, "results.json")

# Load back from JSON
chn_loaded = read_json("results.json")
```

## Saving Chains to JSON

```@docs
write_json
```

### Examples

Save to a specific file:

```julia
using MCMCChains, JSON

chn = Chains(rand(500, 2, 2), [:a, :b])
write_json(chn, "my_chain.json")
```

Save using the model name (if available in chain info):

```julia
# If chain has model_name in info, it will be used as filename
info = (model_name = "gdemo",)
chn = Chains(rand(500, 2, 2), [:s, :m]; info = info)
write_json(chn)  # Creates "gdemo.json"
```

Get JSON as a string instead of writing to file:

```julia
json_string = write_json(chn; as_string = true)
```

## Loading Chains from JSON

```@docs
read_json
```

### Example

```julia
using MCMCChains, JSON

# Load a previously saved chain
chn = read_json("my_chain.json")

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
  "iter": [1, 2, 3, ...],
  "param_names": ["a", "b", "c", "d", "e"],
  "chain_names": ["chain_1", "chain_2"],
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

## Round-Trip Fidelity

The serialization is designed to preserve all information needed to reconstruct the chain:

```julia
using MCMCChains, JSON

# Original chain
chn_original = Chains(
    rand(1000, 3, 2), 
    [:a, :b, :c];
    info = (model_name = "demo", sampler = "NUTS")
)

# Save and load
write_json(chn_original, "temp.json")
chn_loaded = read_json("temp.json")

# Verify data matches
using Test
@test chn_original.value.data ≈ chn_loaded.value.data
@test names(chn_original) == names(chn_loaded)
@test chn_loaded.info.model_name == "demo"
```

## Integration with Turing.jl

The JSON serialization works seamlessly with chains generated from Turing models:

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
write_json(chn)

# Load and analyze
chn_loaded = read_json("gdemo.json")
describe(chn_loaded)
```


