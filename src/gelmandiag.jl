#################### Gelman, Rubin, and Brooks Diagnostics ####################

function gelmandiag(c::AbstractChains; alpha::Real=0.05, mpsrf::Bool=false,
                    transform::Bool=false)
  
    @assert !any(ismissing.(c.value)) "Gelman diagnostics doesn't support missing values"
    
    n, p, m = size(c.value)
  m >= 2 ||
    throw(ArgumentError("less than 2 chains supplied to gelman diagnostic"))

  psi = transform ? link(c) : c.value

  S2 = mapslices(cov, psi, dims = [1, 2])
  W = dropdims(mapslices(mean, S2, dims = [3]), dims = 3)

  psibar = reshape(mapslices(mean, psi, dims = [1]), p, m)'
  B = n * cov(psibar)

  w = diag(W)
  b = diag(B)
  s2 = reshape(mapslices(diag, S2, dims = [1, 2]), p, m)'
  psibar2 = vec(mapslices(mean, psibar, dims = [1]))

  var_w = vec(mapslices(var, s2, dims = [1])) / m
  var_b = (2.0 / (m - 1)) * b.^2
  var_wb = (n / m) * (diag(cov(s2, psibar.^2))
                      - 2.0 * psibar2 .* diag(cov(s2, psibar)))

  V = ((n - 1) / n) * w + ((m + 1) / (m * n)) * b
  var_V = ((n - 1)^2 * var_w + ((m + 1) / m)^2 * var_b +
           (2.0 * (n - 1) * (m + 1) / m) * var_wb) / n^2
  df = 2.0 * V.^2 ./ var_V
  B_df = m - 1
  W_df = 2.0 * w.^2 ./ var_w

  psrf = Array{Float64}(undef, p, 2)
  R_fixed = (n - 1) / n
  R_random_scale = (m + 1) / (m * n)
  q = 1.0 - alpha / 2.0
  for i in 1:p
    correction = (df[i] + 3.0) / (df[i] + 1.0)
    R_random = R_random_scale * b[i] / w[i]
    psrf[i, 1] = sqrt(correction * (R_fixed + R_random))
    if !isnan(R_random)
      R_random *= quantile(FDist(B_df, W_df[i]), q)
    end
    psrf[i, 2] = sqrt(correction * (R_fixed + R_random))
  end
  psrf_labels = ["PSRF", string(100 * q) * "%"]
  psrf_names = c.names

  if mpsrf
    x = isposdef(W) ?
      R_fixed + R_random_scale * eigmax(inv(cholfact(W)) * B) :
      NaN
    psrf = vcat(psrf, [x NaN])
    psrf_names = [psrf_names; "Multivariate"]
  end

  hdr = header(c) * "\nGelman, Rubin, and Brooks Diagnostic:"
  ChainSummary(round.(psrf, digits = 3), psrf_names, psrf_labels, hdr)
end
