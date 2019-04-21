function ess(x::Vector; maxlags = 250)
    m = size(x, 1)
    n = min(size.(x, 3)...)
    maxlags = min(maxlags, n-1)
    lags = collect(0:maxlag)
    #
    # V = Vector(undef, length(lags))
    # for t in lags
    #     V[t] = variogram(x, t, n)
    # end

    vhat = varhat(x, n, m)
    ρ = [ρ_hat(x, varhat, t, n) for t in lags]
    P = ρ_hat(x, varhat, 0, n)
    for T in 1:2:n
        next_rho = ρ_hat(x, varhat, T, n) + ρ_hat(x, varhat, T+1, n)

        if next_rho < 0.0
            P += next_rho
        else
            P += next_rho
            break
        end
    end
    return m*n / (1 + 2*P)
end

function variogram(x::Vector{Vector}, t::Integer, n::Integer)
    rng1 = 1:(n-t)
    rng2 = (1+t):n
    Vt = 0.0
    for j in 1:length(x)
        for i in (t+1):n
            Vt += (x[j][i] - x[j][i-t])^2
        end
    end
    return Vt / (m* (n-t))
end

function varhat(x::Vector{Vector}, n::Integer, m::Integer)
    chain_means = [mean(x[j]) for j in 1:m]
    total_mean = mean(chain_means)
    s = [sum((x[j] .- chain_means[j]).^2)] ./ (n-1)

    B = (n / (m-1)) * sum((chain_means .- total_mean).^2)
    W = mean(s)

    return (n-1)/n * W + 1/n*B
end

function ρ_hat(x::Vector{Vector}, varhat, t::Integer, n::Integer)
    return 1 - variogram(x, t, n) / (2*varhat)
end
