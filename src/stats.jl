#################### Posterior Statistics ####################
"""
    autocor(chn;
        lags=[1, 5, 10, 50],
        demean=true,
        relative=true
        showall=false,
        append_chains=true,
        digits=missing,
        sections=[:parameters])

Compute the autocorrelation of each parameter for the chain. Setting `append_chains=false` will return a vector of dataframes containing the autocorrelations for each chain.

The `digits` keyword may be a(n)
- `Integer`, which sets rounds all numerical columns to the value)
- `NamedTuple`, which only rounds the named column to the specified digits, as with `(mean=2, sd=3)`. This would round the `mean` column to 2 digits and the `sd` column to 3 digits.
- `Dict`, with a similar structure as `NamedTuple`. `Dict(mean => 2, sd => 3)` would set `mean` to two digits and `sd` to three digits.
"""
function autocor(chn::AbstractChains;
        lags::Vector=[1, 5, 10, 50],
        demean::Bool=true,
        relative::Bool=true,
        showall=false,
        append_chains = false,
        sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
        digits=missing,
        kwargs...)
    funs = Function[]
    func_names = String[]
    for i in lags
        push!(funs, x -> autocor(x, [i], demean=demean)[1])
        push!(func_names, "lag $i")
    end
    return summarize(chn, funs...;
        func_names = func_names,
        showall = showall,
        append_chains = append_chains,
        name = "Autocorrelation",
        digits=digits)
end

"""
    cor(chn; showall=false, append_chains=true, sections=[:parameters], digits=missing)

Compute the Pearson correlation matrix for the chain. Setting `append_chains=false` will
return a vector of dataframes containing a correlation matrix for each chain.

The `digits` keyword may be a(n)
- `Integer`, which sets rounds all numerical columns to the value)
- `NamedTuple`, which only rounds the named column to the specified digits, as with `(mean=2, sd=3)`. This would round the `mean` column to 2 digits and the `sd` column to 3 digits.
- `Dict`, with a similar structure as `NamedTuple`. `Dict(mean => 2, sd => 3)` would set `mean` to two digits and `sd` to three digits.
"""
function cor(chn::AbstractChains;
        showall=false,
        append_chains=true,
        sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
        digits=missing,
        kwargs...)
    df = DataFrame(chn, sections, showall=showall, append_chains=append_chains)

    # Generate
    cormat = if append_chains
        cor(convert(Matrix, df[:, 1:end]))
    else
        [cor(convert(Matrix, i[:, 1:end])) for i in df]
    end
    nms = append_chains ? names(df) : names(df[1])
    columns = if append_chains
        [nms, [cormat[:, i] for i in 1:size(cormat, 2)]...]
    else
        [[nms, [cm[:, i] for i in 1:size(cm, 2)]...] for cm in cormat]
    end
    colnames = vcat([:parameters], nms...)
    df_summary = if append_chains
        ChainDataFrame("Correlation", DataFrame(columns, colnames), digits=digits)
    else
        [ChainDataFrame("Correlation", DataFrame(c, colnames), digits=digits)
            for c in columns]
    end
    return df_summary
end

"""
    changerate(chn;
        append_chains=true,
        showall=false,
        sections=[:parameters],
        digits=missing)

Computes the change rate for the chain. Setting `append_chains=false` will
return a vector of dataframes containing the change rates for each chain.

The `digits` keyword may be a(n)
- `Integer`, which sets rounds all numerical columns to the value)
- `NamedTuple`, which only rounds the named column to the specified digits, as with `(mean=2, sd=3)`. This would round the `mean` column to 2 digits and the `sd` column to 3 digits.
- `Dict`, with a similar structure as `NamedTuple`. `Dict(mean => 2, sd => 3)` would set `mean` to two digits and `sd` to three digits.
"""
function changerate(chn::AbstractChains;
    append_chains=true,
    showall=false,
    sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
    digits=missing,
    kwargs...)
    # Check for missing values.
    @assert !any(ismissing.(chn.value)) "Change rate comp. doesn't support missing values."

    df = DataFrame(chn, append_chains=true, showall=showall)
    n, p = size(df[1])
    m = length(chains(chn))

    r = zeros(Float64, p, 1, 1)
    r_mv = 0.0
    delta = Array{Bool}(undef, p)
    for k in 1:m
        dfk = df[k]
        prev = convert(Vector, dfk[1,:])
        for i in 2:n
            for j in 1:p
                x = dfk[i, j]
                dx = x != prev[j]
                r[j] += dx
                delta[j] = dx
                prev[j] = x
            end
            r_mv += any(delta)
        end
    end
    vals = round.([r..., r_mv] / (m * (n - 1)), digits = 3)
    rownames = push!(names(df[1]), :multivariate)
    return DataFrame(parameters = rownames,
        change_rate = vals,
        name="Change Rate",
        digits=digits)
end

describe(c::AbstractChains; args...) = describe(stdout, c; args...)

"""
    describe(io::IO,
        c::AbstractChains;
        q = [0.025, 0.25, 0.5, 0.75, 0.975],
        etype=:bm,
        showall=false,
        sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
        digits=missing,
        args...)

Prints the summary statistics and quantiles for the chain.

The `digits` keyword may be a(n)
- `Integer`, which sets rounds all numerical columns to the value)
- `NamedTuple`, which only rounds the named column to the specified digits, as with `(mean=2, sd=3)`. This would round the `mean` column to 2 digits and the `sd` column to 3 digits.
- `Dict`, with a similar structure as `NamedTuple`. `Dict(mean => 2, sd => 3)` would set `mean` to two digits and `sd` to three digits.
"""
function describe(io::IO,
                  c::AbstractChains;
                  q = [0.025, 0.25, 0.5, 0.75, 0.975],
                  etype=:bm,
                  showall::Bool=false,
                  sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
                  digits=missing,
                  args...
                 )
    dfs = [summarystats(c,
                showall=showall,
                sections=sections,
                etype=etype,
                digits=digits,
                args...),
           quantile(c,
                showall=showall,
                sections=sections,
                q=q,
                digits=digits)]
    return dfs
end

function _hpd(x::Vector{T}; alpha::Real=0.05) where {T<:Real}
    n = length(x)
    m = max(1, ceil(Int, alpha * n))

    y = sort(x)
    a = y[1:m]
    b = y[(n - m + 1):n]
    _, i = findmin(b - a)

    return [a[i], b[i]]
end

function hpd(chn::AbstractChains; alpha::Real=0.05,
        append_chains=true,
        showall=false,
        sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
        digits=nothing)
    labels = [:upper, :lower]
    u(x) = _hpd(x, alpha=alpha)[1]
    l(x) = _hpd(x, alpha=alpha)[2]
    return summarize(chn, u, l;
        func_names = labels,
        showall=showall,
        name="HPD",
        digits=digits)
end

"""
    quantile(chn;
        q=[0.025, 0.25, 0.5, 0.75, 0.975],
        append_chains=true,
        showall=false,
        sections=[:parameters],
        digits=missing)

Computes the quantiles for each parameter in the chain. Setting `append_chains=false` will
return a vector of dataframes containing the quantiles for each chain.

The `digits` keyword may be a(n)
- `Integer`, which sets rounds all numerical columns to the value)
- `NamedTuple`, which only rounds the named column to the specified digits, as with `(mean=2, sd=3)`. This would round the `mean` column to 2 digits and the `sd` column to 3 digits.
- `Dict`, with a similar structure as `NamedTuple`. `Dict(mean => 2, sd => 3)` would set `mean` to two digits and `sd` to three digits.
"""
function quantile(chn::AbstractChains;
        q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975],
        append_chains=true,
        showall=false,
        sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
        digits=missing)
    # compute quantiles
    funs = Function[]
    func_names = String[]
    for i in q
        push!(funs, x -> quantile(cskip(x), i))
        push!(func_names, "$(string(100*i))%")
    end
    return summarize(chn, funs...;
        func_names=func_names,
        showall=showall,
        sections=sections,
        name = "Quantiles",
        digits=digits)
end

"""
	ess(chn::AbstractChains;
		showall=false,
		sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
		maxlag = 250,
        digits=missing)

Compute a chain's number of effective samples. More information can be found in the Gelman et al. (2014) book "Bayesian Data Analysis", or in [this article](https://arxiv.org/abs/1903.08008).

The `digits` keyword may be a(n)
- `Integer`, which sets rounds all numerical columns to the value)
- `NamedTuple`, which only rounds the named column to the specified digits, as with `(mean=2, sd=3)`. This would round the `mean` column to 2 digits and the `sd` column to 3 digits.
- `Dict`, with a similar structure as `NamedTuple`. `Dict(mean => 2, sd => 3)` would set `mean` to two digits and `sd` to three digits.
"""
function ess(chn::AbstractChains;
    showall=false,
    sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
    maxlag = 250,
    digits=missing
)
	param = showall ? names(chn) : names(chn, sections)
	n_chain_orig = size(chn, 3)

	# Split the chains.
	parameter_vec = Vector(undef, length(param))
	midpoint = Int32(floor(size(chn, 1) / 2))
	for i in 1:length(param)
		parameter_vec[i] = []
		for j in 1:n_chain_orig
            c1 = vec(cskip(chn[1:midpoint, param[i], j].value.data))
            c2 = vec(cskip(chn[midpoint+1:end, param[i], j].value.data))
            push!(parameter_vec[i], c1, c2)
		end
	end

    # Misc allocations.
    m = n_chain_orig * 2
    n = min(length.(parameter_vec[1])...)
    maxlag = min(maxlag, n-1)
    lags = collect(0:maxlag)

    # Preallocate B, W, varhat, and Rhat vectors for each param.
    B = Vector(undef, length(param))
    W = Vector(undef, length(param))
    varhat = Vector(undef, length(param))
    Rhat = Vector(undef, length(param))

    # calculate B, W, varhat, and Rhat for each param.
    for i in 1:length(param)
		draw = parameter_vec[i]
		p = param[i]
        allchain = mean(vcat([d for d in draw]...))
        eachchain = [mean(draw[j]) for j in 1:m]
        s = [sum((draw[j] .- eachchain[j]).^2) / (n-1) for j in 1:m]
        B[i] = (n / (m - 1)) * sum((eachchain .- allchain).^2)
        W[i] = sum(s) / m
        varhat[i] = (n-1)/n * W[i] + (B[i] / n)
        Rhat[i] = sqrt(varhat[i] / W[i])
    end

	V = Vector(undef, length(param))
    ρ = Vector(undef, length(param))
    for p in eachindex(V)
        V[p] = Vector(undef, length(lags))
		ρ[p] = Vector(undef, length(lags))
    end

    # Calculate ρ
    c_autocor = Vector(undef, length(param))
    for i in 1:length(param)
        c_autocor[i] = [0.0]
    end

    for t in eachindex(lags)
        lag = lags[t]
        range1 = lag+1:n
        range2 = 1:(n-lag)
        for i in 1:length(param)
            draw = parameter_vec[i]
			p = param[i]
            z = [draw[j][range1] .- draw[j][range2] for j in 1:m]
			z = sum([zi .^ 2 for zi in z])
            V[i][t] = 1 / (m * (n-lag)) * sum(z)
            autocors = [_autocorrelation(draw[j], lag) for j in 1:m]
			ρ[i][t] = 1 - V[i][t] / (2 * varhat[i])
        end
    end

	# Find first odd positive integer where ρ[p][T+1] + ρ[p][T+2] is negative
    P = Vector(undef, length(param))
    ess = Vector(undef, length(param))
	for i in 1:length(param)
        big_P = 0.0
		ρ_val = Float64.(ρ[i])

        # Big P.
        P[i] = Float64[ρ_val[1]]
        k = tprime = 1
        for tprime in 1:Int(floor((length(lags)/2 - 1)))
            sumvals = ρ_val[2*tprime] + ρ_val[2*tprime+1]
            if sumvals < 0
                break
            else
                push!(P[i], sumvals)
                k = tprime
            end
        end

        # Create monotone.
        P_monotone = [min(P[i][t], P[i][1:t]...) for t in 1:length(P[i])]

        ess[i] = (n*m) / (-1 + 2*sum(P_monotone))
	end

    df = DataFrame(parameters = Symbol.(param), ess = ess, r_hat = Rhat)
	return ChainDataFrame("ESS", df, digits=digits)
end

# this function is sourced from https://github.com/tpapp/MCMCDiagnostics.jl/blob/master/src/MCMCDiagnostics.jl
function _autocorrelation(x::AbstractVector, k::Integer, v = var(x), kwargs...)
    x1 = @view(x[1:(end-k)])
    x2 = @view(x[(1+k):end])
    V = sum((x1 .- x2).^2) / length(x1)
    1 - V / (2*v)
end

"""
    summarystats(chn;
        append_chains=true,
        showall=false,
        sections=[:parameters],
        digits=missing,
        args...)

Computes the mean, standard deviation, naive standard error, Monte Carlo standard error, and effective sample size for each parameter in the chain. Setting `append_chains=false` will return a vector of dataframes containing the summary statistics for each chain. `args...` is passed to the `msce` function.

The `digits` keyword may be a(n)
- `Integer`, which sets rounds all numerical columns to the value)
- `NamedTuple`, which only rounds the named column to the specified digits, as with `(mean=2, sd=3)`. This would round the `mean` column to 2 digits and the `sd` column to 3 digits.
- `Dict`, with a similar structure as `NamedTuple`. `Dict(mean => 2, sd => 3)` would set `mean` to two digits and `sd` to three digits.
"""
function summarystats(chn::MCMCChains.AbstractChains;
        append_chains::Bool=true,
        showall::Bool=false,
        sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters],
        etype=:bm,
        digits=missing,
        args...
    )

    # Make some functions.
    sem(x) = sqrt(var(x) / length(x))
    df_mcse(x) = length(x) < 200 ?
        missing :
        mcse(cskip(x), etype, args...)

    # Store everything.
    funs = [mean∘cskip, std∘cskip, sem∘cskip, df_mcse]
    func_names = [:mean, :std, :naive_se, :mcse]

    # Caluclate ESS separately.
    ess_df = ess(chn, sections=sections, showall=showall).df

    # Summarize.
    summary_df = summarize(chn, funs...;
        sections=sections,
        func_names=func_names,
        showall=showall,
        name="Summary Statistics",
        additional_df = ess_df,
        digits=digits)

    return summary_df
end
