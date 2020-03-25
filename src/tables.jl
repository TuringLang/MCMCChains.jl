# Tables and TableTraits interface

## Chains

function _check_columnnames(chn::Chains)
    for name in names(chn)
        symname = Symbol(name)
        if symname === :iteration || symname === :chain
            error("'$(name)' is a reserved column name. Please rename the parameter.")
        end
    end
end

Tables.istable(::Type{<:Chains}) = true

Tables.columnaccess(::Type{<:Chains}) = true

function Tables.columns(chn::Chains)
    _check_columnnames(chn)
    return chn
end

Tables.columnnames(chn::Chains) = (:iteration, :chain, Symbol.(names(chn))...)

function Tables.getcolumn(chn::Chains, i::Int)
    return Tables.getcolumn(chn, Tables.columnnames(chn)[i])
end
function Tables.getcolumn(chn::Chains, nm::Symbol)
    if nm == :iteration
        iterations = range(chn)
        nchains = size(chn, 3)
        return repeat(iterations, nchains)
    elseif nm == :chain
        chainids = chains(chn)
        niter = size(chn, 1)
        return repeat(chainids; inner = niter)
    else
        return vec(getindex(chn, nm).value)
    end
end

Tables.rowaccess(::Type{<:Chains}) = true

function Tables.rows(chn::Chains)
    _check_columnnames(chn)
    return chn
end

Tables.rowtable(chn::Chains) = Tables.rowtable(Tables.columntable(chn))

function Tables.namedtupleiterator(chn::Chains)
    return Tables.namedtupleiterator(Tables.columntable(chn))
end

function Tables.schema(chn::Chains)
    _check_columnnames(chn)
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

Tables.columnnames(::ChainDataFrame{<:NamedTuple{names}}) where {names} = names

Tables.getcolumn(cdf::ChainDataFrame, i::Int) = cdf.nt[i]
Tables.getcolumn(cdf::ChainDataFrame, nm::Symbol) = cdf.nt[nm]

Tables.rowaccess(::Type{<:ChainDataFrame}) = true

Tables.rows(cdf::ChainDataFrame) = cdf

Tables.rowtable(cdf::ChainDataFrame) = Tables.rowtable(Tables.columntable(cdf))

function Tables.namedtupleiterator(cdf::ChainDataFrame)
    return Tables.namedtupleiterator(Tables.columntable(cdf))
end

function Tables.schema(::ChainDataFrame{NamedTuple{names,T}}) where {names,T}
    return Tables.Schema(names, eltype.(T.parameters))
end

IteratorInterfaceExtensions.isiterable(::ChainDataFrame) = true
function IteratorInterfaceExtensions.getiterator(cdf::ChainDataFrame)
    return Tables.datavaluerows(Tables.columntable(cdf))
end

TableTraits.isiterabletable(::ChainDataFrame) = true
