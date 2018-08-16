#################### Raftery and Lewis Diagnostic ####################

function rafterydiag{T<:Real}(x::Vector{T}; q::Real=0.025, r::Real=0.005,
                              s::Real=0.95, eps::Real=0.001,
                              range::Range=1:length(x))
  nx = length(x)
  phi = sqrt(2.0) * erfinv(s)
  nmin = ceil(Int, q * (1.0 - q) * (phi / r)^2)
  if nmin > nx
    warn("At least $nmin samples are needed for specified q, r, and s")
    kthin = burnin = total = NaN
  else
    dichot = Int[(x .<= quantile(x, q))...]
    kthin = 0
    bic = 1.0
    local test, ntest
    while bic >= 0.0
      kthin += 1
      test = dichot[1:kthin:nx]
      ntest = length(test)
      temp = test[1:(ntest - 2)] + 2 * test[2:(ntest - 1)] + 4 * test[3:ntest]
      trantest = reshape(counts(temp, 0:7), 2, 2, 2)
      g2 = 0.0
      for i1 in 1:2, i2 in 1:2, i3 in 1:2
        tt = trantest[i1, i2, i3]
        if tt > 0
          fitted = sum(trantest[:, i2, i3]) * sum(trantest[i1, i2, :]) /
                   sum(trantest[:, i2, :])
          g2 += 2.0 * tt * log(tt / fitted)
        end
      end
      bic = g2 - 2.0 * log(ntest - 2.0)
    end
    tranfinal = counts(test[1:(ntest - 1)] + 2 * test[2:ntest], 0:3)
    alpha = tranfinal[3] / (tranfinal[1] + tranfinal[3])
    beta = tranfinal[2] / (tranfinal[2] + tranfinal[4])
    kthin *= step(range)
    m = log(eps * (alpha + beta) / max(alpha, beta)) /
        log(abs(1.0 - alpha - beta))
    burnin = kthin * ceil(m) + start(range) - 1
    n = ((2.0 - alpha - beta) * alpha * beta * phi^2) /
        (r^2 * (alpha + beta)^3)
    keep = kthin * ceil(n)
    total = burnin + keep
  end
  [kthin, burnin, total, nmin, total / nmin]
end

function rafterydiag(c::AbstractChains; q::Real=0.025, r::Real=0.005,
                     s::Real=0.95, eps::Real=0.001)
  _, p, m = size(c.value)
  vals = Array{Float64}(p, 5, m)
  for j in 1:p, k in 1:m
    vals[j, :, k] = rafterydiag(c.value[:, j, k], q=q, r=r, s=s, eps=eps,
                                range=c.range)
  end
  hdr = header(c) * "\nRaftery and Lewis Diagnostic:\n" *
        "Quantile (q) = $q\nAccuracy (r) = $r\nProbability (s) = $s\n"
  ChainSummary(vals, c.names, ["Thinning", "Burn-in", "Total", "Nmin",
                               "Dependence Factor"], hdr)
end
