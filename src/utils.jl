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
  p / (pi^1.5 * sqrt(q))
end
