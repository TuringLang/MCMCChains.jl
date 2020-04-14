#################### Posterior Statistics ####################
"""
    autocor(chains[; lags = [1, 5, 10, 50], demean = true, append_chains = true, kwargs...])

Compute the autocorrelation of each parameter for the chain.

Setting `append_chains=false` will return a vector of dataframes containing the
autocorrelations for each chain.
"""
function autocor(
    chains::Chains;
    lags::AbstractVector{<:Integer} = [1, 5, 10, 50],
    demean::Bool = true,
    append_chains = false,
    kwargs...
)
    funs = Function[]
    func_names = @. Symbol("lag ", lags)
    for i in lags
        push!(funs, x -> autocor(x, [i], demean=demean)[1])
    end

    return summarize(
        chains, funs...;
        func_names = func_names,
        append_chains = append_chains,
        name = "Autocorrelation",
        kwargs...
    )
end

"""
    cor(chains[; sections = :parameters, append_chains = true, kwargs...])

Compute the Pearson correlation matrix for the chain.

Setting `append_chains=false` will return a vector of dataframes containing a correlation
matrix for each chain.
"""
function cor(
    chains::Chains;
    sections = :parameters,
    append_chains = true,
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Obstain names of parameters.
    names_of_params = names(_chains)

    if append_chains
        df = chaindataframe_cor("Correlation", names_of_params, to_matrix(_chains))
        return df
    else
        vector_of_df = [
            chaindataframe_cor(
                "Correlation - Chain $i", names_of_params, data
            )
            for (i, data) in enumerate(to_vector_of_matrices(_chains))
        ]
        return vector_of_df
    end
end

function chaindataframe_cor(name, names_of_params, chains::AbstractMatrix; kwargs...)
    # Compute the correlation matrix.
    cormat = cor(chains)

    # Summarize the results in a named tuple.
    nt = (; parameters = names_of_params,
          zip(names_of_params, (cormat[:, i] for i in axes(cormat, 2)))...)

    # Create a ChainDataFrame.
    return ChainDataFrame(name, nt; kwargs...)
end

"""
    changerate(chains[; sections = :parameters, append_chains = true, kwargs...])

Compute the change rate for the chain.

Setting `append_chains=false` will return a vector of dataframes containing the change
rates for each chain.
"""
function changerate(
    chains::Chains{<:Real};
    sections = :parameters,
    append_chains = true,
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Obstain names of parameters.
    names_of_params = names(_chains)

    if append_chains
        df = chaindataframe_changerate("Change Rate", names_of_params, _chains.value.data)
        return df
    else
        vector_of_df = [
            chaindataframe_changerate(
                "Change Rate - Chain $i", names_of_params, data
            )
            for (i, data) in enumerate(to_vector_of_matrices(_chains))
        ]
        return vector_of_df
    end
end

function chaindataframe_changerate(name, names_of_params, chains; kwargs...)
    # Compute the change rates.
    changerates, mvchangerate = chains

    # Summarize the results in a named tuple.
    nt = (; zip(names_of_params, changerates)..., multivariate = mvchangerate)

    # Create a ChainDataFrame.
    return ChainDataFrame(name, nt; kwargs...)
end

changerate(chains::AbstractMatrix{<:Real}) = changerate(reshape(chains, Val(3)))
function changerate(chains::AbstractArray{<:Real,3})
    niters, nparams, nchains = size(chains)

    changerates = zeros(nparams)
    mvchangerate = 0.0

    for chain in 1:nchains, iter in 2:niters
        isanychanged = false

        for param in 1:nparams
            # update if the sample is different from the one in the previous iteration
            if chains[iter-1, param, chain] != chains[iter, param, chain]
                changerates[param] += 1
                isanychanged = true
            end
        end

        mvchangerate += isanychanged
    end

    factor = nchains * (niters - 1)
    changerates ./= factor
    mvchangerate /= factor

    changerates, mvchangerate
end

describe(c::Chains; args...) = describe(stdout, c; args...)

"""
    describe(io, chains[;
             q = [0.025, 0.25, 0.5, 0.75, 0.975],
             etype = :bm,
             kwargs...])

Print the summary statistics and quantiles for the chain.
"""
function describe(
    io::IO,
    chains::Chains;
    q = [0.025, 0.25, 0.5, 0.75, 0.975],
    etype = :bm,
    kwargs...
)
    dfs = vcat(summarystats(chains; etype = etype, kwargs...),
               quantile(chains; q = q, kwargs...))
    return dfs
end

function _hpd(x::AbstractVector{<:Real}; alpha::Real=0.05)
    n = length(x)
    m = max(1, ceil(Int, alpha * n))

    y = sort(x)
    a = y[1:m]
    b = y[(n - m + 1):n]
    _, i = findmin(b - a)

    return [a[i], b[i]]
end

function hpd(chn::Chains; alpha::Real=0.05, kwargs...)
    labels = [:upper, :lower]
    u(x) = _hpd(x, alpha=alpha)[1]
    l(x) = _hpd(x, alpha=alpha)[2]
    return summarize(chn, u, l; name = "HPD", func_names = labels, kwargs...)
end

"""
    quantile(chains[; q = [0.025, 0.25, 0.5, 0.75, 0.975], append_chains = true, kwargs...])

Compute the quantiles for each parameter in the chain.

Setting `append_chains=false` will return a vector of dataframes containing the quantiles
for each chain.
"""
function quantile(
    chains::Chains;
    q::AbstractVector = [0.025, 0.25, 0.5, 0.75, 0.975],
    append_chains = true,
    kwargs...
)
    # compute quantiles
    funs = Function[]
    func_names = @. Symbol(100 * q, :%)
    for i in q
        push!(funs, x -> quantile(cskip(x), i))
    end

    return summarize(
        chains, funs...;
        func_names = func_names,
        append_chains = append_chains,
        name = "Quantiles",
        kwargs...
    )
end

"""
	  ess(chn::Chains[;	sections = :parameters, kwargs...])

Compute the effective sample size (ESS) of a chain.

More information can be found in the Gelman et al. (2014) book "Bayesian Data Analysis", or
in [this article](https://arxiv.org/abs/1903.08008).
"""
function ess(
    chains::Chains;
    sections = :parameters,
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Compute the ESS.
    nt = merge((parameters = names(_chains),), ess(_chains.value.data; kwargs...))

	  return ChainDataFrame("ESS", nt)
end

function ess(
    chn::AbstractArray{<:Union{Real,Missing},3};
    maxlag = 250
)
    n, nparams, n_chain_orig = size(chn)

	# Split the chains.
	parameter_vec = Vector(undef, nparams)
    midpoint = Int32(floor(size(chn, 1) / 2))
    for i in 1:nparams
		parameter_vec[i] = []
		for j in 1:n_chain_orig
            c1 = vec(cskip(chn[1:midpoint, i, j]))
            c2 = vec(cskip(chn[midpoint+1:end, i, j]))
            n = min(n, length(c1), length(c2))
            push!(parameter_vec[i], c1, c2)
		end
	end

    # Misc allocations.
    m = n_chain_orig * 2
    maxlag = min(maxlag, n-1)
    lags = 0:maxlag

    # Preallocate B, W, varhat, and Rhat vectors for each param.
    B = Vector(undef, nparams)
    W = Vector(undef, nparams)
    varhat = Vector(undef, nparams)
    Rhat = Vector{Float64}(undef, nparams)

    # calculate B, W, varhat, and Rhat for each param.
    for i in 1:nparams
		draw = parameter_vec[i]
        allchain = mean(vcat([d for d in draw]...))
        eachchain = [mean(draw[j]) for j in 1:m]
        s = [sum((draw[j] .- eachchain[j]).^2) / (n-1) for j in 1:m]
        B[i] = (n / (m - 1)) * sum((eachchain .- allchain).^2)
        W[i] = sum(s) / m
        varhat[i] = (n-1)/n * W[i] + (B[i] / n)
        Rhat[i] = sqrt(varhat[i] / W[i])
    end

	V = Vector(undef, nparams)
    ρ = Vector(undef, nparams)
    for p in eachindex(V)
        V[p] = Vector(undef, length(lags))
		ρ[p] = Vector(undef, length(lags))
    end

    # Calculate ρ
    c_autocor = Vector(undef, nparams)
    for i in 1:nparams
        c_autocor[i] = [0.0]
    end

    for t in eachindex(lags)
        lag = lags[t]
        range1 = lag+1:n
        range2 = 1:(n-lag)
        for i in 1:nparams
            draw = parameter_vec[i]
            z = [draw[j][range1] .- draw[j][range2] for j in 1:m]
			z = sum([zi .^ 2 for zi in z])
            V[i][t] = 1 / (m * (n-lag)) * sum(z)
            autocors = [_autocorrelation(draw[j], lag) for j in 1:m]
			ρ[i][t] = 1 - V[i][t] / (2 * varhat[i])
        end
    end

	# Find first odd positive integer where ρ[p][T+1] + ρ[p][T+2] is negative
    P = Vector(undef, nparams)
    ess = Vector{Float64}(undef, nparams)
	for i in 1:nparams
        big_P = 0.0
		ρ_val = Float64.(ρ[i])

        # Big P.
        P[i] = Float64[ρ_val[1]]
        k = tprime = 1
        for tprime in 1:Int(floor((length(lags)/2 - 1)))
            sumvals = ρ_val[2*tprime] + ρ_val[2*tprime+1]
            push!(P[i], sumvals)
            k = tprime
            
            if sumvals < 0
                break
            end
        end

        # Create monotone.
        P_monotone = [min(P[i][t], P[i][1:t]...) for t in 1:length(P[i])]

        ess[i] = (n*m) / (-1 + 2*sum(P_monotone))
    end

    return (ess = ess, rhat = Rhat)
end

# this function is sourced from https://github.com/tpapp/MCMCDiagnostics.jl/blob/master/src/MCMCDiagnostics.jl
function _autocorrelation(x::AbstractVector, k::Integer, v = var(x), kwargs...)
    x1 = @view(x[1:(end-k)])
    x2 = @view(x[(1+k):end])
    V = sum((x1 .- x2).^2) / length(x1)
    1 - V / (2*v)
end

"""
    summarystats(chains[;
                 sections = :parameters, append_chains = true, maxlag = 250, etype = :bm,
                 kwargs...])

Compute the mean, standard deviation, naive standard error, Monte Carlo standard error,
and effective sample size for each parameter in the chain.

Setting `append_chains=false` will return a vector of dataframes containing the summary
statistics for each chain.

The remaining keyword arguments are passed to the `msce` function.
"""
function summarystats(
    chains::Chains;
    sections = :parameters,
    append_chains::Bool = true,
    maxlag = 250,
    etype = :bm,
    kwargs...
)
    # Make some functions.
    df_mcse(x) = length(x) < 200 ?
        missing :
        mcse(cskip(x), etype; kwargs...)

    # Store everything.
    funs = [mean∘cskip, std∘cskip, sem∘cskip, df_mcse]
    func_names = [:mean, :std, :naive_se, :mcse]

    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Caluclate ESS separately.
    ess_df = ess(_chains; sections = nothing, maxlag = maxlag)

    # Summarize.
    summary_df = summarize(
        _chains, funs...;
        func_names = func_names,
        append_chains = append_chains,
        additional_df = ess_df,
        name = "Summary Statistics",
        sections = nothing
    )

    return summary_df
end

"""
    mean(chains[, params; kwargs...])

Calculate the mean of a chain.
"""
function mean(chains::Chains; kwargs...)
    # Store everything.
    funs = [mean∘cskip]
    func_names = [:mean]

    # Summarize.
    summary_df = summarize(
        chains, funs...;
        func_names = func_names,
        name = "Mean",
        kwargs...
    )

    return summary_df
end

mean(chn::Chains, syms) = mean(chn[:, syms, :])
# resolve method ambiguity with `mean(f, ::AbstractArray)`
mean(chn::Chains, syms::AbstractVector) = mean(chn[:, syms, :])
