module MCMCChainsJSONExt

using MCMCChains
using JSON

import MCMCChains: write_json, read_json

# Helpers for JSON serialization
_clean_json(x) = string(x)
_clean_json(x::Union{Number,String,Bool,Nothing}) = x
_clean_json(x::Symbol) = string(x)
_clean_json(x::Missing) = nothing
_clean_json(x::AbstractArray) = map(_clean_json, x)
_clean_json(x::AbstractDict) =
    Dict{String,Any}(string(k) => _clean_json(v) for (k, v) in pairs(x))
_clean_json(x::NamedTuple) =
    Dict{String,Any}(string(k) => _clean_json(v) for (k, v) in pairs(x))

function write_json(
    c::Chains,
    filepath::Union{AbstractString,Nothing} = nothing;
    as_string::Bool = false,
)
    data = Dict{String,Any}()

    # Flatten data for consistent JSON structure
    data["size"] = collect(size(c))
    data["value_flat"] = vec(c.value.data)

    data["iterations"] = collect(range(c))
    data["parameters"] = map(string, names(c))
    data["chains"] = collect(chains(c))
    data["logevidence"] = ismissing(c.logevidence) ? nothing : c.logevidence

    data["name_map"] =
        Dict{String,Any}(string(k) => map(string, v) for (k, v) in pairs(c.name_map))

    data["info"] = _clean_json(c.info)

    if as_string
        return JSON.json(data)
    end

    if isnothing(filepath)
        raw_name = get(c.info, :model_name, "chain")
        name = string(raw_name)
        filepath = "$(name).json"
    end

    open(filepath, "w") do f
        JSON.print(f, data)
    end

    return filepath
end

function read_json(filepath::AbstractString)
    data = JSON.parsefile(filepath)

    dims = Tuple(data["size"])
    raw_vec = data["value_flat"]

    # Handle JSON nulls as missing
    convert_val(x) = x === nothing ? missing : x
    val_vec = map(convert_val, raw_vec)

    # Reconstruct 3D array
    val_typed = convert(Vector{Union{Missing,Float64}}, val_vec)
    val = reshape(val_typed, dims)

    nms = Symbol.(data["parameters"])

    # Reconstruct Name Map
    nm_dict = data["name_map"]
    nm_pairs = [Symbol(k) => Symbol.(v) for (k, v) in nm_dict]
    name_map = (; nm_pairs...)

    # Reconstruct Info
    # Note: Complex objects in info (like Models) are loaded as Strings due to sanitization
    info_dict = data["info"]
    info_pairs = [Symbol(k) => v for (k, v) in info_dict]
    info = (; info_pairs...)

    iters = data["iterations"]
    start_val = isempty(iters) ? 1 : iters[1]
    step_val = length(iters) > 1 ? iters[2] - iters[1] : 1

    logev = data["logevidence"] === nothing ? missing : data["logevidence"]

    return Chains(
        val,
        nms,
        name_map;
        start = start_val,
        thin = step_val,
        evidence = logev,
        info = info,
    )
end

end # module
