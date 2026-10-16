function _value_counts(val::AbstractMatrix)
    T = Base.nonmissingtype(eltype(val))
    frequencies = [Dict{T,Int}() for _ in axes(val, 2)]
    levels = Set{T}()
    for (counts, draws) in zip(frequencies, eachcol(val))
        for v in skipmissing(draws)
            isfinite(v) || throw(ArgumentError("countplot requires finite sample values"))
            v = iszero(v) ? zero(v) : v
            counts[v] = get(counts, v, 0) + 1
            push!(levels, v)
        end
    end
    isempty(levels) && throw(ArgumentError("countplot requires at least one sample value"))
    ordered = sort!(collect(levels))
    counts = [get(frequency, v, 0) for v in ordered, frequency in frequencies]
    return ordered, counts
end
