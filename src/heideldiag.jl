#################### Heidelberger and Welch Diagnostic ####################

function heideldiag(x::Vector{T}; alpha::Real=0.05, eps::Real=0.1,
                             etype=:imse, start::Integer=1, args...) where {T<:Real}
  n = length(x)
  delta = trunc(Int, 0.10 * n)
  y = x[trunc(Int, n / 2):end]
  S0 = length(y) * mcse(y, etype; args...)^2
  i, pvalue, converged, ybar = 1, 1.0, false, NaN
  while i < n / 2
    y = x[i:end]
    m = length(y)
    ybar = mean(y)
    B = cumsum(y) - ybar * collect(1:m)
    Bsq = (B .* B) ./ (m * S0)
    I = sum(Bsq) / m
    pvalue = 1.0 - pcramer(I)
    converged = pvalue > alpha
    if converged
      break
    end
    i += delta
  end
  halfwidth = sqrt(2.0) * erfinv(1.0 - alpha) * mcse(y, etype; args...)
  passed = halfwidth / abs(ybar) <= eps
  [i + start - 2, converged, round(pvalue, digits = 4), ybar, halfwidth, passed]
end

function heideldiag(c::AbstractChains;
                    alpha = 0.05,
                    eps = 0.1,
                    etype = :imse,
                    args...
                   )
    # Preallocate.
    _, p, m = size(c.value)
    vals = Array{Float64}(undef, p, 6, m)

    # Perform tests.
    for j in 1:p, k in 1:m
        vals[j, :, k] = heideldiag(
                            collect(skipmissing(c.value[:, j, k])),
                            alpha=alpha,
                            eps=eps,
                            etype=etype,
                            start=first(c.range);
                            args...
                           )
    end

    # Set header.
    hdr = header(c) * "\nHeidelberger and Welch Diagnostic:\n" *
        "Target Halfwidth Ratio = $eps\nAlpha = $alpha\n"

    # Round values.
    vals = map(x -> round(x, digits=4), vals)

    # Make a ChainSummary object.
    ChainSummary(vals, c.names, ["Burn-in", "Stationarity", "p-value", "Mean",
                               "Halfwidth", "Test"], hdr)
end
