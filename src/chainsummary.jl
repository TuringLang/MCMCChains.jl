#################### ChainSummary ####################

#################### Types and Constructors ####################
struct ChainSummary
    value::Array{Float64, 3}
    rownames::Vector{AbstractString}
    colnames::Vector{AbstractString}
    header::AbstractString
    sorted::Bool

    function ChainSummary(value::Array{Float64, 3},
            rownames::Vector{AbstractString},
            colnames::Vector{AbstractString},
            header::AbstractString,
            sorted = false)
        dim = size(value)
        length(rownames) == dim[1] ||
            throw(DimensionMismatch("numbers of rownames and rows differ"))
        length(colnames) == dim[2] ||
            throw(DimensionMismatch("numbers of colnames and columns differ"))
        return new(value, rownames, colnames, header, sorted)
    end
end

function ChainSummary(value::Array{Float64, 3},
                      rownames::Vector{T},
                      colnames::Vector{U},
                      header::AbstractString,
                      sorted = false) where {T<:AbstractString, U<:AbstractString}
  ChainSummary(copy(value), AbstractString[rownames...],
               AbstractString[colnames...], header, sorted)
end

function ChainSummary(value::Matrix{Float64},
                      rownames::Vector{T},
                      colnames::Vector{U},
                      header::AbstractString,
                      sorted = false) where {T<:AbstractString, U<:AbstractString}
  dim = size(value)
  ChainSummary(reshape(value, dim[1], dim[2], 1), AbstractString[rownames...],
               AbstractString[colnames...], header, sorted)
end

struct ChainSummaries
    summaries::Vector{ChainSummary}
    ChainSummaries(v::Vector) = new(v)
end

#################### Base Methods ####################

# Return a new, sorted ChainSummary.
function sort(c::ChainSummary)
    sorted_names = sort(c.rownames, lt=natural)
    indices = indexin(sorted_names, c.rownames)
    new_rows = sorted_names
    new_value = copy(c.value)
    for i in eachindex(indices)
        new_value[i, :] = c.value[indices[i], :]
    end
    return ChainSummary(new_value, new_rows, c.colnames, c.header, c.sorted)
end

# Returns the maximum width of a chain summary when printed.
function max_width(s::ChainSummary)
    rnwid = map(length, s.rownames)
    mxrnwid = maximum(rnwid)
    charv = mapslices(showoff, s.value, dims = [1])
    cnwid = map(length, s.colnames)
    colwid = 1 .+ max.(cnwid, vec(maximum(map(length, charv), dims = [1, 3] )))
    return sum(colwid) + mxrnwid
end

## write n ' ' characters to io
wrtsp(io::IO, n) = while (n -= 1) >= 0 write(io, ' ') end

function Base.show(io::IO, cs::ChainSummaries)
    if length(cs.summaries) == 1
        show(io, cs.summaries[1])
    else
        for i in cs.summaries
            println(i.header)
            show(io, i)
            println()
        end
    end
end

function Base.show(io::IO, s::ChainSummary)
    # Sort the summary if needed.
    if s.sorted
        s = sort(s)
    end

    ## rowname widths
    rnwid = map(length, s.rownames)
    mxrnwid = maximum(rnwid)
    ## column name widths
    cnwid = map(length, s.colnames)
    ## s.value as right alignable strings
    charv = mapslices(showoff, s.value, dims = [1])
    colwid = 1 .+ max.(cnwid, vec(maximum(map(length, charv), dims = [1, 3] )))
    m, n, f = size(charv)
    for k in 1:f
        ## write the column headers centered on the column widths
        wrtsp(io, mxrnwid)
        for j in 1:n
            ## do not count the leading space
            nspace = colwid[j] - cnwid[j] - 1
            ## divide by 2 rounding down
            nright = nspace >> 1
            wrtsp(io, 1 + nspace - nright)
            print(io, s.colnames[j])
            wrtsp(io, nright)
        end
        println(io)
        for i in 1:m
            wrtsp(io, mxrnwid - rnwid[i])
            print(io, s.rownames[i])
            for j in 1:n
                wrtsp(io, colwid[j] - length(charv[i, j, k]))
                print(io, charv[i, j, k])
            end
            println(io)
        end
        println(io)
    end
end
