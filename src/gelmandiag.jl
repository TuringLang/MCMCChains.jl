#################### Gelman, Rubin, and Brooks Diagnostics ####################

function _gelmandiag(
    psi::AbstractArray{<:Real,3};
    alpha::Real = 0.05
)
    niters, nparams, nchains = size(psi)
    nchains > 1 || error("Gelman diagnostic requires at least 2 chains")

    rfixed = (niters - 1) / niters
    rrandomscale = (nchains + 1) / (nchains * niters)

    S2 = mapslices(cov, psi, dims = [1, 2])
    W = dropdims(mapslices(mean, S2, dims = [3]), dims = 3)

    psibar = reshape(mapslices(mean, convert(Array, psi), dims = [1]), nparams, nchains)'
    B = niters .* cov(psibar)

    w = diag(W)
    b = diag(B)
    s2 = reshape(mapslices(diag, S2, dims = [1, 2]), nparams, nchains)'
    psibar2 = vec(mapslices(mean, psibar, dims = [1]))

    var_w = vec(mapslices(var, s2, dims = [1])) ./ nchains
    var_b = (2 / (nchains - 1)) .* b.^2
    var_wb = (niters / nchains) .*
        (diag(cov(s2, psibar.^2)) .- 2 .* psibar2 .* diag(cov(s2, psibar)))

    V = @. rfixed * w + rrandomscale * b
    var_V = rfixed^2 * var_w + rrandomscale^2 * var_b + 2 * rfixed * rrandomscale * var_wb
    
    df = @. 2 * V^2 / var_V

    B_df = nchains - 1
    W_df = @. 2 * w^2 / var_w

    estimates = Array{Float64}(undef, nparams)
    upperlimits = Array{Float64}(undef, nparams)

    q = 1 - alpha / 2
    for i in 1:nparams
        correction = (df[i] + 3) / (df[i] + 1)
        rrandom = rrandomscale * b[i] / w[i]

        estimates[i] = sqrt(correction * (rfixed + rrandom))

        if !isnan(rrandom)
            rrandom *= quantile(FDist(B_df, W_df[i]), q)
        end
        upperlimits[i] = sqrt(correction * (rfixed + rrandom))
    end

    return estimates, upperlimits, W, B
end

"""
    gelmandiag(chains::AbstractArray{<:Real,3}; kwargs...)

Gelman, Rubin and Brooks diagnostics.
"""
function gelmandiag(chains::AbstractArray{<:Real,3}; kwargs...)
    estimates, upperlimits = _gelmandiag(chains; kwargs...)

    return (psrf = estimates, psrfci = upperlimits)
end

"""
    gelmandiag_multivariate(chains::AbstractArray{<:Real,3}; kwargs...)

Multivariate Gelman, Rubin and Brooks diagnostics.
"""
function gelmandiag_multivariate(
    chains::AbstractArray{<:Real,3};
    kwargs...
)
    niters, nparams, nchains = size(chains)
    if nparams < 2
        error("computation of the multivariate potential scale reduction factor requires ",
              "at least two variables")
    end

    estimates, upperlimits, W, B = _gelmandiag(chains; kwargs...)

    rfixed = (niters - 1) / niters
    rrandomscale = (nchains + 1) / (nchains * niters)
    multivariate = rfixed + rrandomscale * LinearAlgebra.eigmax(W \ B)

    return (psrf = estimates, psrfci = upperlimits, psrfmultivariate = multivariate)
end

function gelmandiag(
    chains::Chains{<:Real};
    sections = _default_sections(chains),
    transform = false,
    alpha = 0.05,
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Compute the potential scale reduction factor.
    psi = transform ? link(_chains) : _chains.value.data
    results = gelmandiag(psi; alpha = alpha, kwargs...)

    # Create a named tuple with the results.
    colnames = (:psrf, Symbol(100 * (1 - alpha / 2), :%))
    nt = (; parameters = names(_chains), zip(colnames, results)...)

    return ChainDataFrame("Gelman, Rubin, and Brooks Diagnostic", nt)
end

function gelmandiag_multivariate(
    chains::Chains{<:Real};
    sections = _default_sections(chains),
    transform = true,
    alpha = 0.05,
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Compute the potential scale reduction factor.
    psi = transform ? link(_chains) : _chains.value.data
    results = gelmandiag_multivariate(psi; alpha = alpha, kwargs...)

    # Create a named tuple with the results.
    colnames = (:psrf, Symbol(100 * (1 - alpha / 2), :%))
    nt = (; parameters = names(_chains), zip(colnames, (results.psrf, results.psrfci))...)

    return ChainDataFrame("Gelman, Rubin, and Brooks Diagnostic", nt),
        results.multivariate
end
