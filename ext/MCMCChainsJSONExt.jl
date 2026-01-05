module MCMCChainsJSONExt

using MCMCChains
using JSON

import MCMCChains: write_json, read_json

function JSON.lower(c::Chains)
    return (
        size = size(c),
        value_flat = vec(c.value.data),
        iterations = range(c),
        parameters = names(c),
        chains = chains(c),
        logevidence = c.logevidence,
        name_map = c.name_map,
        info = _lower_info(c.info),
    )
end

_lower_info(x) = Dict{String,Any}()
function _lower_info(info::NamedTuple)
    result = Dict{String,Any}()
    for (k, v) in pairs(info)
        result[string(k)] = _lower_value(v)
    end
    return result
end

_lower_value(x) = string(x)
_lower_value(x::Union{Number,String,Bool,Nothing}) = x
_lower_value(x::Symbol) = string(x)
_lower_value(x::Missing) = nothing
_lower_value(x::AbstractArray) = map(_lower_value, x)
_lower_value(x::AbstractDict) =
    Dict{String,Any}(string(k) => _lower_value(v) for (k, v) in pairs(x))
_lower_value(x::NamedTuple) =
    Dict{String,Any}(string(k) => _lower_value(v) for (k, v) in pairs(x))

function write_json(
    c::Chains,
    filepath::Union{AbstractString,Nothing} = nothing;
    as_string::Bool = false,
)
    if as_string
        return JSON.json(c)
    end

    if isnothing(filepath)
        raw_name = get(c.info, :model_name, "chain")
        name = string(raw_name)
        filepath = "$(name).json"
    end

    JSON.json(filepath, c)

    return filepath
end

function read_json(filepath::AbstractString)
    data = JSON.parsefile(filepath; null = missing)

    dims = Tuple(data["size"])
    raw_vec = data["value_flat"]

    # Reconstruct 3D array
    val_typed = convert(Vector{Union{Missing,Float64}}, raw_vec)
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

    logev = data["logevidence"]

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
