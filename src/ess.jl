# function ess(draw::Vector;
#         n = min(size.(x, 3)...)
#         m = n*2
#         maxlags = min(250, n-1),
#         )
#     # Collect all the lags.
#     lags = collect(0:maxlags)
#
#     allchain = mean(vcat([d for d in draw]...))
#     eachchain = Vector(undef, m)
#     s = Vector(undef, m)
#     for j in 1:m
#         eachchain[j] = mean(draw[j])
#         s[j] = (1/(n-1)) * sum((draw[j] .- eachchain[j]).^2)
#     end
#     B = (n / (m - 1)) * sum((eachchain .- allchain).^2)
#     W = (1/m) * sum(s)
#     vhat = (n-1)/n * W[i] + (1/n) * B[i]
#     Rhat = sqrt(varhat[i] / W[i])
#
#     V = Vector(undef, length(lags))
#     ρ = Vector(undef, length(lags))
#
#     for t in eachindex(lags)
#         sum_autocor =
#         # Old way.
#         # lag = lags[t]
#         # range1 = lag+1:n
#         # range2 = 1:(n-lag)
#         #
#         # draw = parameter_vec[i]
# 		# p = param[i]
#         # z = [draw[j][range1] .- draw[j][range2] for j in 1:m]
# 		# z = sum([zi .^ 2 for zi in z])
#         # V[i][t] = 1 / (m * (n-lag)) * sum(z)
#         # autocors = [autocor(draw[j], [lag])[1] for j in 1:m]
#         # # ρ[i][t] = 1 - (W[i] - sum(autocors))/varhat[i]
# 		# ρ[i][t] = 1 - V[i][t] / (2 * varhat[i])
#     end
# end
#
# function variogram(x::Vector{Vector{T}}, t::Integer, n::Integer, m::Integer) where T<:Real
#     rng1 = 1:(n-t)
#     rng2 = (1+t):n
#     Vt = 0.0
#     for j in 1:length(x)
#         for i in (t+1):n
#             Vt += (x[j][i] - x[j][i-t])^2
#         end
#     end
#     return Vt / (m* (n-t))
# end
#
# function varhat(x::Vector{Vector{T}}, n::Integer, m::Integer) where T<:Real
#     chain_means = [mean(x[j]) for j in 1:m]
#     total_mean = mean(chain_means)
#     s = [sum((x[j] .- chain_means[j]).^2) for j in 1:m] ./ (n-1)
#
#     B = (n / (m-1)) * sum((chain_means .- total_mean).^2)
#     W = mean(s)
#
#     return (n-1)/n * W + 1/n*B
# end
#
# function ρ_hat(x::Vector{Vector{T}}, vhat, t::Integer, n::Integer, m::Integer) where T<:Real
#     return 1 - variogram(x, t, n, m) / (2*vhat)
# end
