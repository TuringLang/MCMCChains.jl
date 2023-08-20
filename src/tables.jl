# Tables and TableTraits interface

####
#### Chains
####

function _check_columnnames(chn::Chains)
    for name in names(chn)
        symname = Symbol(name)
        if symname === :iteration || symname === :chain
            error("'$(name)' is a reserved column name. Please rename the parameter.")
        end
    end
end

#### Tables interface

Tables.istable(::Type{<:Chains}) = true

# AbstractColumns interface

Tables.columnaccess(::Type{<:Chains}) = true

function Tables.columns(chn::Chains)
    _check_columnnames(chn)
    return chn
end

Tables.columnnames(chn::Chains) = (:iteration, :chain, names(chn)...)

function Tables.getcolumn(chn::Chains, i::Int)
    return Tables.getcolumn(chn, Tables.columnnames(chn)[i])
end
function Tables.getcolumn(chn::Chains, nm::Symbol)
    if nm === :iteration
        iterations = range(chn)
        nchains = size(chn, 3)
        return repeat(iterations, nchains)
    elseif nm === :chain
        chainids = chains(chn)
        niter = size(chn, 1)
        return repeat(chainids; inner = niter)
    else
        return vec(chn[nm])
    end
end

# row access

Tables.rowaccess(::Type{<:Chains}) = true

Tables.rows(chn::Chains) = Tables.rows(Tables.columntable(chn))

# optional Tables overloads

function Tables.schema(chn::Chains)
    _check_columnnames(chn)
    nms = Tables.columnnames(chn)
    T = eltype(chn.value)
    types = (Int, Int, ntuple(_ -> T, size(chn, 2))...)
    return Tables.Schema(nms, types)
end

#### TableTraits interface

IteratorInterfaceExtensions.isiterable(::Chains) = true
function IteratorInterfaceExtensions.getiterator(chn::Chains)
    return Tables.datavaluerows(Tables.columntable(chn))
end

TableTraits.isiterabletable(::Chains) = true
