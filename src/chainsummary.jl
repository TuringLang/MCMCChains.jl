#################### ChainSummary ####################

#################### Types and Constructors ####################

immutable ChainSummary
  value::Array{Float64, 3}
  rownames::Vector{AbstractString}
  colnames::Vector{AbstractString}
  header::AbstractString

  function ChainSummary(value::Array{Float64, 3},
                        rownames::Vector{AbstractString},
                        colnames::Vector{AbstractString},
                        header::AbstractString)
    dim = size(value)
    length(rownames) == dim[1] ||
      throw(DimensionMismatch("numbers of rownames and rows differ"))
    length(colnames) == dim[2] ||
      throw(DimensionMismatch("numbers of colnames and columns differ"))
    new(value, rownames, colnames, header)
  end
end

function ChainSummary(value::Array{Float64, 3}, rownames::Vector{T},
                     colnames::Vector{U}, header::AbstractString) where {T<:AbstractString, U<:AbstractString}
  ChainSummary(copy(value), AbstractString[rownames...],
               AbstractString[colnames...], header)
end

function ChainSummary(value::Matrix{Float64}, rownames::Vector{T},
                     colnames::Vector{U}, header::AbstractString) where {T<:AbstractString, U<:AbstractString}
  dim = size(value)
  ChainSummary(reshape(value, dim[1], dim[2], 1), AbstractString[rownames...],
               AbstractString[colnames...], header)
end


#################### Base Methods ####################

function Base.showall(io::IO, s::ChainSummary)
  println(io, s.header)
  show(io, s)
end

## write n ' ' characters to io
wrtsp(io::IO, n) = while (n -= 1) >= 0 write(io, ' ') end

function Base.show(io::IO, s::ChainSummary)
  ## rowname widths
  rnwid = map(length, s.rownames)
  mxrnwid = maximum(rnwid)
  ## column name widths
  cnwid = map(length, s.colnames)
  ## s.value as right alignable strings
  charv = mapslices(showoff, s.value, 1)
  colwid = 1 + max.(cnwid, vec(maximum(map(length, charv), [1, 3])))
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
