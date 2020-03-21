_as_vec(::Type{T}, x::T) where {T} = T[x]
_as_vec(::Type{T}, x::NTuple{N,T}) where {N,T} = T[x...]
_as_vec(::Type{T}, x::AbstractArray{T}) where {T} = vec(x)

struct TableChains{T<:Chains,S<:Vector{Symbol},C<:Vector{<:Integer},M}
    data::T
    sections::S
    chains::C
    coltypes::M
end

function TableChains(
    data::Chains;
    sections = sections(data),
    chains = chains(data),
    coltypes = NamedTuple(),
)
    sections = _clean_sections(data, _as_vec(Symbol, sections))
    chains = _as_vec(Int, chains) ∩ collect(MCMCChains.chains(data))
    return TableChains{typeof.((data, sections, chains, coltypes))...}(
        data,
        sections,
        chains,
        coltypes,
    )
end
@inline TableChains(t::TableChains) = t

const TableOrChains = Union{Chains,TableChains}

Tables.istable(::Type{<:TableOrChains}) = true

Tables.columnaccess(::Type{<:TableOrChains}) = true

Tables.columns(t::TableOrChains) = TableChains(t)

function Tables.columnnames(t::TableOrChains)
    t = TableChains(t)
    return Symbol[:Iteration; :Chain; Symbol.(names(t.data, t.sections))]
end

function Tables.getcolumn(t::TableOrChains, nm::Symbol)
    t = TableChains(t)
    if nm == :Iteration
        niter = size(t.data, 1)
        nchains = length(t.chains)
        return repeat(1:niter, nchains)
    elseif nm == :Chain
        niter = size(t.data, 1)
        return vcat((fill(c, niter) for c in t.chains)...)
    else
        val = getindex(t.data, :, nm, t.chains).value
        T = get(t.coltypes, nm, nothing)
        T === nothing && return vec(val)
        return vec(convert(Array{T}, val))
    end
end

Tables.rowaccess(::Type{<:TableOrChains}) = true

Tables.rows(t::TableOrChains) = TableChains(t)

Tables.rowtable(t::TableOrChains) = Tables.rowtable(Tables.columntable(t))

function Tables.namedtupleiterator(t::TableOrChains)
    return Tables.namedtupleiterator(Tables.columntable(t))
end

function Tables.schema(t::TableOrChains)
    t = TableChains(t)
    nms = Tables.columnnames(t)
    types = DataType[Int; Int; get.(Ref(t.coltypes), nms[3:end], Ref(eltype(t.data.value)))]
    return Tables.Schema(nms, types)
end

IteratorInterfaceExtensions.isiterable(::TableOrChains) = true
function IteratorInterfaceExtensions.getiterator(t::TableOrChains)
    return Tables.datavaluerows(Tables.columntable(t))
end

TableTraits.isiterabletable(::TableOrChains) = true
