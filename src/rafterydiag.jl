#################### Raftery and Lewis Diagnostic ####################

"""
    rafterydiag(x::Vector{<:Real}; q, r, s, eps, range)
    rafterydiag(chains::Chains; sections, q, r, s, eps)

Raftery and Lewis Diagnostic.
"""
function rafterydiag(
                     x::Vector{<:Real};
                     q = 0.025,
                     r = 0.005,
                     s = 0.95,
                     eps = 0.001,
                     range = 1:length(x)
                    )

    nx = length(x)
    phi = sqrt(2.0) * erfinv(s)
    nmin = ceil(Int, q * (1.0 - q) * (phi / r)^2)
    if nmin > nx
        @warn "At least $nmin samples are needed for specified q, r, and s"
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
        m = log(eps * (alpha + beta) / max(alpha, beta)) / log(abs(1.0 - alpha - beta))
        burnin = kthin * ceil(m) + first(range) - 1
        n = ((2.0 - alpha - beta) * alpha * beta * phi^2) / (r^2 * (alpha + beta)^3)
        keep = kthin * ceil(n)
        total = burnin + keep
    end
    return [kthin, burnin, total, nmin, total / nmin]
end

function rafterydiag(
    chains::Chains;
    sections = _default_sections(chains),
    q = 0.025,
    r = 0.005,
    s = 0.95,
    eps = 0.001
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    _, p, m = size(_chains.value)
    vals = [Array{Float64}(undef, p, 5) for i in 1:m]
    for j in 1:p, k in 1:m
        vals[k][j, :] = rafterydiag(
            collect(skipmissing(_chains.value[:, j, k])),
            q=q,
            r=r,
            s=s,
            eps=eps,
            range=range(_chains)
        )
    end

    # Retrieve columns.
    data = [[vals[k][:, i] for i in 1:5] for k in 1:m]

    # Obtain names of parameters.
    names_of_params = names(_chains)

    # Compute data frames.
    vector_of_df = [
        ChainDataFrame(
            "Raftery and Lewis Diagnostic - Chain $i",
            (parameters = names_of_params, thinning = columns[1], burnin = columns[2],
             total = columns[3], nmin = columns[4], dependencefactor = columns[5])
        )
        for (i, columns) in enumerate(data)
    ]

    return vector_of_df
end
