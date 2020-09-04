using .XGBoost
"""
    rstar(chains::Chains; subset = 0.8, niter = 1_000, eta = 0.5, XGBoostParams)
    rstar(chains::Chains, iterations::Int; subset = 0.8, niter = 1_000, eta = 0.5, XGBoostParams)
    rstar(x::AbstractMatrix, y::AbstractVector, nchains::Int, iterations::Int; subset = 0.8, niter = 1_000, eta = 0.5, XGBoostParams)

Compute the R* statistic for convergence diagnostic of MCMC. This implementation is an adaption of Algorithm 1 & 2, described in [Lambert & Vehtari]. Note that the correctness of the statistic depends on the convergence of the classifier used internally in the statistic. You can track if the training of the classifier converged by inspection of the printed RMSE values from the XGBoost backend. To adjust the number of iterations used to train the classifier set `niter` accordingly.

# Usage

```julia-repl
using XGBoost
# You need to load XGBoost before using MCMCChains.Rstar

...

chn = ...

# Compute R* using defaults settings for the  gradient boosting classifier used to compute the statistic.
R = rstar(chn)

# Compute 100 samples of the R* statistic using sampling from according to the prediction probabilities.
# This approach can be slow and results in a less accurate estimation of the R* statistic.
# See discussion in Section 3.1.3 in the paper.
Rs = rstar(chn, 100)

# estimate Rstar
R = mean(Rs)

# visualize distribution
histogram(Rs)
```

## References:
[Lambert & Vehtari] Ben Lambert and Aki Vehtari. "R∗: A robust MCMC convergence diagnostic with uncertainty using gradient-boostined machines." Arxiv 2020.
"""
function rstar(rng::Random.AbstractRNG, x::AbstractMatrix, y::AbstractVector, nchains::Int; iterations = 1, subset = 0.8, opt_iter = 1_000, xgboostparams...)

    size(x,1) != length(y) && throw(DimensionMismatch())

    N = length(y)

    # randomly sub-select training and testing set
    Ntrain = round(Int, N*subset)
    Ntest = N - Ntrain

    ids = Random.randperm(rng, N)
    train_ids = view(ids, 1:Ntrain)
    test_ids = view(ids, (Ntrain+1):N)

    # use predicted probabilities?
    mode = iterations > 1 ? "multi:softprob" : "multi:softmax"

    # train classifier using XGBoost
    classif = XGBoost.xgboost(x[train_ids,:], opt_iter; label = y[train_ids], 
                      objective = mode, num_class = nchains,
                      xgboostparams...)

    Rstats = ones(iterations) * Inf
    ytest = view(y, test_ids)
    for i in 1:iterations
        # predict labels for "test" data
        p = XGBoost.predict(classif, x[test_ids,:])

        pred = if length(p) == Ntest*nchains
            # Note: XGBoost outputs a nchains × Ntest array.
            # For some reasons the docu says something else.
            catrand(rng, reshape(p, nchains, Ntest)') .- 1
        else
            p
        end

        # compute statistic
        a = mean(((p,y),) -> p == y, zip(pred, ytest))
        Rstats[i] = nchains*a
    end

    return Rstats
end

function rstar(chn::Chains; kwargs...)
    nchains = size(chn, 3)
    nchains <= 1 && throw(DimensionMismatch())

    # collect data
    x = Array(chn)
    y = repeat(chains(chn); inner = size(chn,1)) .- 1

    return rstar(Random.default_rng(), x, y, nchains; kwargs...)
end

"""
    catrand(rng, probs)

A helper function to obtain random samples a bit more efficiently.
The `probs` matrix is expected to be N × K and doesn't have to be normalised.
"""
function catrand(rng::Random.AbstractRNG, probs::AbstractMatrix)
    N,K = size(probs)

    u = Random.rand(rng, N)
    z = zeros(Int, N)
    @inbounds for n in 1:N
        k = 1
        p = view(probs,n,:)
        c = p[1]
        while c < u[n] && k < K
            c += p[k += 1]
        end
        z[n] = k
    end
    return z
end
