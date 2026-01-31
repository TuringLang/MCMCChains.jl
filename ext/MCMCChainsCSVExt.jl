module MCMCChainsCSVExt

using MCMCChains
using CSV
using Tables
using Dates

const STAN_SAMPLER_PARAMS = Set([
    :lp__,
    :accept_stat__,
    :stepsize__,
    :treedepth__,
    :n_leapfrog__,
    :divergent__,
    :energy__,
])

function _convert_param_name_to_stan(name::Symbol)
    s = string(name)
    s = replace(s, "[" => ".")
    s = replace(s, "]" => "")
    s = replace(s, "," => ".")
    return s
end

function _convert_param_name_from_stan(name::AbstractString)
    if !contains(name, ".")
        return Symbol(name)
    end
    parts = split(name, ".")
    base = parts[1]
    idx_parts = parts[2:end]
    if !isempty(idx_parts) && all(p -> !isempty(p) && all(isdigit, p), idx_parts)
        indices = join(idx_parts, ",")
        return Symbol("$base[$indices]")
    end
    return Symbol(name)
end

function _get_stan_column_order(chn::Chains)
    all_names = names(chn)
    internals = get(chn.name_map, :internals, Symbol[])
    parameters = get(chn.name_map, :parameters, Symbol[])

    sampler_cols = [n for n in all_names if n in STAN_SAMPLER_PARAMS && n in internals]
    other_internals = [n for n in internals if n ∉ sampler_cols]
    ordered_params = [n for n in all_names if n in parameters]

    return vcat(sampler_cols, ordered_params, other_internals)
end

function _write_adaptation(io::IO, chn::Chains, chain_id::Int)
    println(io, "# Adaptation terminated")
    if :stepsize__ in names(chn)
        stepsize = chn[end, :stepsize__, chain_id]
        println(io, "# Step size = ", stepsize)
    end
    println(io, "# Diagonal elements of inverse mass matrix:")
    println(io, "# 1.0")
end

function _write_timing(io::IO, chn::Chains, chain_id::Int)
    println(io, "# ")
    start_t = MCMCChains.start_times(chn)
    stop_t = MCMCChains.stop_times(chn)

    has_timing =
        start_t !== missing &&
        stop_t !== missing &&
        chain_id <= length(start_t) &&
        chain_id <= length(stop_t) &&
        start_t[chain_id] !== missing &&
        stop_t[chain_id] !== missing

    if has_timing
        duration = Dates.value(stop_t[chain_id] - start_t[chain_id]) / 1000
        warmup = get(chn.info, :warmup_time, 0.0)
        println(io, "#  Elapsed Time: ", warmup, " seconds (Warm-up)")
        println(io, "#                ", duration, " seconds (Sampling)")
        println(io, "#                ", warmup + duration, " seconds (Total)")
    else
        println(io, "#  Elapsed Time: unknown")
    end
    println(io, "# ")
end

struct StanCSVTable
    chn::Chains
    chain_id::Int
    col_order::Vector{Symbol}
end

Tables.istable(::Type{StanCSVTable}) = true
Tables.columnaccess(::Type{StanCSVTable}) = true
Tables.columns(t::StanCSVTable) = t

function Tables.columnnames(t::StanCSVTable)
    return Tuple(Symbol(_convert_param_name_to_stan(n)) for n in t.col_order)
end

function Tables.getcolumn(t::StanCSVTable, i::Int)
    return Tables.getcolumn(t, t.col_order[i])
end

function Tables.getcolumn(t::StanCSVTable, nm::Symbol)
    stan_name = string(nm)
    for orig_name in t.col_order
        if _convert_param_name_to_stan(orig_name) == stan_name
            return vec(t.chn.value[:, orig_name, t.chain_id])
        end
    end
    return vec(t.chn.value[:, _convert_param_name_from_stan(stan_name), t.chain_id])
end

function Tables.schema(t::StanCSVTable)
    nms = Tables.columnnames(t)
    T = eltype(t.chn.value)
    return Tables.Schema(nms, ntuple(_ -> T, length(nms)))
end

function MCMCChains.write_stancsv(
    file::AbstractString,
    chn::Chains;
    chain_id::Int = 1,
    include_adaptation::Bool = true,
    include_timing::Bool = true,
    kwargs...,
)
    col_order = _get_stan_column_order(chn)
    table = StanCSVTable(chn, chain_id, col_order)

    open(file, "w") do io
        include_adaptation && _write_adaptation(io, chn, chain_id)
    end

    CSV.write(file, table; append = true, header = true, kwargs...)

    if include_timing
        open(file, "a") do io
            _write_timing(io, chn, chain_id)
        end
    end

    return file
end

function MCMCChains.write_stancsv(
    basename::AbstractString,
    chn::Chains,
    ::Val{:all};
    kwargs...,
)
    n_chains = size(chn, 3)
    files = String[]

    for chain_id = 1:n_chains
        base, ext = splitext(basename)
        ext = isempty(ext) ? ".csv" : ext
        file = n_chains == 1 ? basename : "$(base)_$(chain_id)$(ext)"
        MCMCChains.write_stancsv(file, chn; chain_id, kwargs...)
        push!(files, file)
    end

    return files
end

function MCMCChains.read_stancsv(file::AbstractString; kwargs...)
    lines = readlines(file)
    header_idx = findfirst(l -> !startswith(l, "#"), lines)
    header_idx === nothing && error("No header found in StanCSV file")

    data_lines =
        filter(l -> !startswith(l, "#") && !isempty(strip(l)), lines[header_idx:end])
    csv_data = CSV.File(IOBuffer(join(data_lines, "\n")); kwargs...)

    csv_names = Tables.columnnames(csv_data)
    col_names = [_convert_param_name_from_stan(string(n)) for n in csv_names]

    n_iter, n_params = length(csv_data), length(col_names)
    val = Matrix{Union{Missing,Float64}}(undef, n_iter, n_params)

    for (j, csv_name) in enumerate(csv_names)
        col = Tables.getcolumn(csv_data, csv_name)
        for (i, v) in enumerate(col)
            val[i, j] = ismissing(v) ? missing : Float64(v)
        end
    end

    internals =
        [n for n in col_names if endswith(string(n), "__") || n in STAN_SAMPLER_PARAMS]
    parameters = [n for n in col_names if n ∉ internals]
    name_map = (parameters = parameters, internals = internals)

    info_dict = Dict{Symbol,Any}()
    for line in lines[1:(header_idx-1)]
        m = match(r"^#\s*model\s*=\s*(.+?)(?:\s+\(Default\))?$", line)
        m !== nothing && (info_dict[:model_name] = strip(m[1]); continue)
        m = match(r"^#\s*Step size\s*=\s*([\d.eE+-]+)", line)
        m !== nothing && (info_dict[:stepsize] = parse(Float64, m[1]); continue)
        m = match(r"^#\s+seed\s*=\s*(\d+)", line)
        m !== nothing && (info_dict[:seed] = parse(Int, m[1]); continue)
        m = match(r"^#\s+num_warmup\s*=\s*(\d+)", line)
        m !== nothing && (info_dict[:num_warmup] = parse(Int, m[1]); continue)
    end

    return Chains(
        reshape(val, n_iter, n_params, 1),
        col_names,
        name_map;
        info = NamedTuple(info_dict),
    )
end

function MCMCChains.read_stancsv(files::AbstractVector{<:AbstractString}; kwargs...)
    return MCMCChains.chainscat([MCMCChains.read_stancsv(f; kwargs...) for f in files]...)
end

function MCMCChains.Chains(file::CSV.File)
    col_names = collect(Tables.columnnames(file))
    n_iter, n_params = length(file), length(col_names)

    val = Matrix{Union{Missing,Float64}}(undef, n_iter, n_params)
    for (j, nm) in enumerate(col_names)
        col = Tables.getcolumn(file, nm)
        for (i, v) in enumerate(col)
            val[i, j] = ismissing(v) ? missing : Float64(v)
        end
    end

    return Chains(val, Symbol.(col_names))
end

end
