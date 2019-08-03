#################### Geweke Diagnostic ####################

function gewekediag(x::Vector{T}; first::Real=0.1, last::Real=0.5,
                             etype=:imse, args...) where {T<:Real}
  if !(0.0 < first < 1.0)
    throw(ArgumentError("first is not in (0, 1)"))
  elseif !(0.0 < last < 1.0)
    throw(ArgumentError("last is not in (0, 1)"))
  elseif first + last > 1.0
    throw(ArgumentError("first and last proportions overlap"))
  end
  n = length(x)
  x1 = x[1:round(Int, first * n)]
  x2 = x[round(Int, n - last * n + 1):n]
  z = (mean(x1) - mean(x2)) /
      sqrt(mcse(x1, etype; args...)^2 + mcse(x2, etype; args...)^2)
  [round(z, digits = 3), round(1.0 - erf(abs(z) / sqrt(2.0)), digits = 4)]
end

function gewekediag(chn::AbstractChains; first::Real=0.1, last::Real=0.5,
                    etype=:imse,
                    sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
                    showall=false,
                    sorted=true, 
                    args...)
    c = showall ?
        sorted ? sort(chn) : chn :
        Chains(chn, _clean_sections(chn, sections); sorted=sorted)

    _, p, m = size(c.value)
    diags = [vals = Array{Float64}(undef, p, 2) for _ in 1:m]
    for j in 1:p, k in 1:m
        diags[k][j,:] = gewekediag(
                            collect(skipmissing(c.value[:, j, k])),
                            first=first,
                            last=last,
                            etype=etype;
                            args...
                        )
    end

    name = "Geweke Diagnostic\nFirst Window Fraction = $first\n" *
        "Second Window Fraction = $last\n"

    pnames = Symbol.(names(c))
    colnames = [:parameters, :z_score, :p_value]
    columns = [[pnames, d[:,1], d[:,2]] for d in diags]
    dfs = [ChainDataFrame("Geweke Diagnostic - Chain $i",
                          DataFrame(columns[i], colnames)) for i in 1:m]

    return dfs
end
