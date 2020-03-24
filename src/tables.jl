# Tables and TableTraits interface

## Chains

Tables.istable(::Type{<:Chains}) = true

Tables.columnaccess(::Type{<:Chains}) = true

Tables.columns(chn::Chains) = chn

Tables.columnnames(chn::Chains) = (:Iteration, :Chain, Symbol.(names(chn))...)

function Tables.getcolumn(chn::Chains, i::Int)
    return Tables.getcolumn(chn, Tables.columnnames(chn)[i])
end
function Tables.getcolumn(chn::Chains, nm::Symbol)
    if nm == :Iteration
        niter, _, nchains = size(chn)
        return repeat(Base.OneTo(niter), nchains)
    elseif nm == :Chain
        chainids = chains(chn)
        niter = size(chn, 1)
        return vcat(fill.(chainids, niter)...)
    else
        return vec(getindex(chn, nm).value)
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

## ChainDataFrame

Tables.istable(::Type{<:ChainDataFrame}) = true

Tables.columnaccess(::Type{<:ChainDataFrame}) = true

Tables.columns(cdf::ChainDataFrame) = cdf

Tables.columnnames(cdf::ChainDataFrame) = keys(cdf.nt)

function Tables.getcolumn(cdf::ChainDataFrame, i::Int)
    return Tables.getcolumn(cdf, keys(cdf.nt)[i])
end
Tables.getcolumn(cdf::ChainDataFrame, nm::Symbol) = getproperty(cdf.nt, nm)

Tables.rowaccess(::Type{<:ChainDataFrame}) = true

Tables.rows(cdf::ChainDataFrame) = cdf

Tables.rowtable(cdf::ChainDataFrame) = Tables.rowtable(Tables.columntable(cdf))

function Tables.namedtupleiterator(cdf::ChainDataFrame)
    return Tables.namedtupleiterator(Tables.columntable(cdf))
end

function Tables.schema(cdf::ChainDataFrame)
    return Tables.Schema(keys(cdf.nt), eltype.(values(cdf.nt)))
end

IteratorInterfaceExtensions.isiterable(::ChainDataFrame) = true
function IteratorInterfaceExtensions.getiterator(cdf::ChainDataFrame)
    return Tables.datavaluerows(Tables.columntable(cdf))
end

TableTraits.isiterabletable(::ChainDataFrame) = true
