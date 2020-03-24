# Tables and TableTraits interface


Tables.istable(::Type{<:Chains}) = true

Tables.columnaccess(::Type{<:Chains}) = true

Tables.columns(chn::Chains) = chn

function Tables.columnnames(chn::Chains)
    return Symbol[:Iteration; :Chain; Symbol.(names(chn))]
end

function Tables.getcolumn(chn::Chains, i::Integer)
    return Tables.getcolumn(chn, Tables.columnnames(chn)[i])
end
function Tables.getcolumn(chn::Chains, nm::Symbol)
    chainids = chains(chn)
    if nm == :Iteration
        niter = size(chn, 1)
        return repeat(1:niter, length(chainids))
    elseif nm == :Chain
        niter = size(chn, 1)
        return vcat((fill(c, niter) for c in chainids)...)
    else
        nchains = length(chainids)
        val = getindex(chn, :, nm, 1:nchains).value
        return vec(val)
    end
end

Tables.rowaccess(::Type{<:Chains}) = true

Tables.rows(chn::Chains) = chn

Tables.rowtable(chn::Chains) = Tables.rowtable(Tables.columntable(chn))

function Tables.namedtupleiterator(chn::Chains)
    return Tables.namedtupleiterator(Tables.columntable(chn))
end

function Tables.schema(chn::Chains)
    nms = Tables.columnnames(chn)
    T = eltype(chn.value)
    types = (Int, Int, ntuple(_ -> T, size(chn, 2))...)
    return Tables.Schema(nms, types)
end

IteratorInterfaceExtensions.isiterable(::Chains) = true
function IteratorInterfaceExtensions.getiterator(chn::Chains)
    return Tables.datavaluerows(Tables.columntable(chn))
end

TableTraits.isiterabletable(::Chains) = true
