struct ChainDataFrame{NT<:NamedTuple}
    name::String
    nt::NT
    nrows::Int
    ncols::Int

    function ChainDataFrame(name::String, nt::NamedTuple)
        lengths = length(first(nt))
        all(x -> length(x) == lengths, nt) || error("Lengths must be equal.")

        return new{typeof(nt)}(name, nt, lengths, length(nt))
    end
end

ChainDataFrame(nt::NamedTuple) = ChainDataFrame("", nt)

Base.size(c::ChainDataFrame) = (c.nrows, c.ncols)
Base.names(c::ChainDataFrame) = collect(keys(c.nt))

# Display

function Base.show(io::IO, df::ChainDataFrame)
    print(io, df.name, " (", df.nrows, " x ", df.ncols, ")")
end

function Base.show(io::IO, ::MIME"text/plain", df::ChainDataFrame)
    digits = get(io, :digits, 4)
    formatter = PrettyTables.ft_printf("%.$(digits)f")

    println(io, df.name)
    # Support for PrettyTables 0.9 (`borderless`) and 0.10 (`tf_borderless`)
    PrettyTables.pretty_table(
        io, df.nt;
        formatters = formatter,
        tf = isdefined(PrettyTables, :borderless) ? PrettyTables.borderless : PrettyTables.tf_borderless,
    )
end

Base.isequal(c1::ChainDataFrame, c2::ChainDataFrame) = isequal(c1, c2)

# Index functions
function Base.getindex(c::ChainDataFrame, s::Union{Colon, Integer, UnitRange}, g::Union{Colon, Integer, UnitRange})
    convert(Array, getindex(c, c.nt[:parameters][s], collect(keys(c.nt))[g]))
end

Base.getindex(c::ChainDataFrame, s::Vector{Symbol}, ::Colon) = getindex(c, s)
function Base.getindex(c::ChainDataFrame, s::Union{Symbol, Vector{Symbol}})
    getindex(c, s, collect(keys(c.nt)))
end

function Base.getindex(c::ChainDataFrame, s::Union{Colon, Integer, UnitRange}, ks)
    getindex(c, c.nt[:parameters][s], ks)
end

# dispatches involing `String` and `AbstractVector{String}`
Base.getindex(c::ChainDataFrame, s::String, ks) = getindex(c, Symbol(s), ks)
function Base.getindex(c::ChainDataFrame, s::AbstractVector{String}, ks)
    return getindex(c, Symbol.(s), ks)
end

# dispatch for `Symbol`
Base.getindex(c::ChainDataFrame, s::Symbol, ks) = getindex(c, [s], ks)

function Base.getindex(c::ChainDataFrame, s::AbstractVector{Symbol}, ks::Symbol)
    return getindex(c, s, [ks])
end

function Base.getindex(
    c::ChainDataFrame,
    s::AbstractVector{Symbol},
    ks::AbstractVector{Symbol}
)
    ind = indexin(s, c.nt[:parameters])

    not_found = map(x -> x === nothing, ind)

    any(not_found) && error("Cannot find parameters $(s[not_found]) in chain")

    # If there are multiple columns, return a new CDF.
    if length(ks) > 1
        if !(:parameters in ks)
            ks = vcat(:parameters, ks)
        end
        nt = NamedTuple{tuple(ks...)}(tuple([c.nt[k][ind] for k in ks]...))
        return ChainDataFrame(c.name, nt)
    else
        # Otherwise, return a vector if there's multiple parameters
        # or just a scalar if there's one parameter.
        if length(s) == 1
            return c.nt[ks[1]][ind][1]
        else
            return c.nt[ks[1]][ind]
        end
    end
end

function Base.lastindex(c::ChainDataFrame, i::Integer)
    if i == 1
        return c.nrows
    elseif i ==2
        return c.ncols
    else
        error("No such dimension")
    end
end

function Base.convert(::Type{Array}, c::C) where C<:ChainDataFrame
    T = promote_eltype_namedtuple_tail(c.nt)
    arr = Array{T, 2}(undef, c.nrows, c.ncols - 1)

    for (i, k) in enumerate(Iterators.drop(keys(c.nt), 1))
        arr[:, i] = c.nt[k]
    end

    return arr
end

function Base.convert(::Type{Array}, cs::Array{ChainDataFrame{T},1}) where T<:NamedTuple
    return mapreduce((x, y) -> cat(x, y; dims = Val(3)), cs) do c
        reshape(convert(Array, c), Val(3))
    end
end

"""
    summarize(
        chains[, stats_funs...];
        append_chains=true,
        name="SummaryStats",
        [sections, var_names],
    )

Summarize `chains` in a `PosteriorStats.SummaryStats`.

`stats_funs` is a collection of functions that reduces a matrix with shape `(draws, chains)`
to a scalar or a collection of scalars. Alternatively, an item in `stats_funs` may be a
`Pair` of the form `name => fun` specifying the name to be used for the statistic or of the
form `(name1, ...) => fun` when the function returns a collection. When the function returns
a collection, the names in this latter format must be provided.

If no stats functions are provided, then those specified in [`default_summary_stats`](@ref)
are computed.

`var_names` specifies the names of the parameters in data. If not provided, the names are
inferred from data.

# Examples

* `summarize(chns)` : Complete chain summary
* `summarize(chns[[:parm1, :parm2]])` : Chain summary of selected parameters
* `summarize(chns; sections=[:parameters])`  : Chain summary of :parameters section
* `summarize(chns; sections=[:parameters, :internals])` : Chain summary for multiple sections
"""
function PosteriorStats.summarize(
    chains::Chains, funs...;
    sections = _default_sections(chains),
    append_chains::Bool = true,
    var_names=nothing,
    kwargs...
)
    # Generate a chain to work on.
    chn = Chains(chains, _clean_sections(chains, sections))

    # Obtain names of parameters.
    names_of_params = var_names === nothing ? names(chn) : var_names

    if append_chains
        # Evaluate the functions.
        data = _permutedims_diagnostics(chains.value.data)
        summarize(data, funs...; var_names=names_of_params, kwargs...)
    else
        # Evaluate the functions.
        data = to_vector_of_matrices(chn)
        return map(data) do x
            z = reshape(x, size(x, 1), 1, size(x, 2))
            summarize(z, funs...; var_names=names_of_params, kwargs...)
        end
    end
end
