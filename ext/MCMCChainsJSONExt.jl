module MCMCChainsJSONExt

using MCMCChains
using JSON

const StructUtils = JSON.StructUtils

StructUtils.structlike(::JSON.JSONStyle, ::Type{<:Chains}) = false

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

function StructUtils.lift(::JSON.JSONStyle, ::Type{Chains}, d::AbstractDict)
    dims = Tuple(d["size"])
    raw_vec = d["value_flat"]
    val_typed = Vector{Union{Missing,Float64}}(undef, length(raw_vec))
    for (i, x) in enumerate(raw_vec)
        val_typed[i] = x === nothing ? missing : x
    end
    val = reshape(val_typed, dims)

    nms = Symbol.(d["parameters"])

    nm_dict = d["name_map"]
    nm_pairs = [Symbol(k) => Symbol.(v) for (k, v) in nm_dict]
    name_map = (; nm_pairs...)
    info_dict = d["info"]
    info_pairs = [Symbol(k) => v for (k, v) in info_dict]
    info = (; info_pairs...)

    iters = d["iterations"]
    start_val = isempty(iters) ? 1 : iters[1]
    step_val = length(iters) > 1 ? iters[2] - iters[1] : 1

    logev = d["logevidence"]

    return Chains(
        val,
        nms,
        name_map;
        start = start_val,
        thin = step_val,
        evidence = logev,
        info = info,
    ), nothing
end

end # module
