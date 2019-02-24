#################### Mathematical Operators ####################
function cummean(x::AbstractArray)
    return mapslices(cummean, x, dims = 1)
end

function cummean(x::AbstractVector)
    y = similar(x, Float64)
    xs = 0.0
    fill!(y, xs)

    c = 0
    for i in 1:length(x)
        if !ismissing(x[i])
            c += 1
            xs += x[i]
        end
        y[i] = xs / c
    end
    return y
end

## Csorgo S and Faraway JJ. The exact and asymptotic distributions of the
## Cramer-von Mises statistic. Journal of the Royal Statistical Society,
## Series B, 58: 221-234, 1996.
function pcramer(q::Real)
    p = 0.0
    for k in 0:3
        c1 = 4.0 * k + 1.0
        c2 = c1^2 / (16.0 * q)
        p += gamma(k + 0.5) / factorial(k) * sqrt(c1) * exp(-c2) * besselk(0.25, c2)
    end
    return p / (pi^1.5 * sqrt(q))
end

# This sorting stack was sourced from https://rosettacode.org/wiki/Natural_sorting#Julia
splitbynum(x) = split(x, r"(?<=\D)(?=\d)|(?<=\d)(?=\D)")
numstringtonum(arr) = [(n = tryparse(Float32, e)) != nothing ? n : e for e in arr]
function natural(x, y)
    xarr = numstringtonum(splitbynum(x))
    yarr = numstringtonum(splitbynum(y))
    for i in 1:min(length(xarr), length(yarr))
        if typeof(xarr[i]) != typeof(yarr[i])
            a = string(xarr[i]); b = string(yarr[i])
        else
             a = xarr[i]; b = yarr[i]
        end
        if a == b
            continue
        else
            return a < b
        end
    end
    return length(xarr) < length(yarr)
end

function _dict2namedtuple(d::Dict)
    t_keys = ntuple(x -> Symbol(collect(keys(d))[x]), length(keys(d)))
    t_vals = ntuple(x -> collect(values(d))[x], length(values(d)))
    return NamedTuple{t_keys}(t_vals)
end

function _namedtuple2dict(n::NamedTuple)
    return Dict(pairs(n))
end
