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

# Compute R⋆ using defaults settings for the  gradient boosting classifier used to compute the statistic.
# This is the recomended use.
R = rstar(chn)

# Compute 100 samples of the R⋆ statistic using sampling from according to the prediction probabilities.
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
function rstar(x::AbstractMatrix, y::AbstractVector, nchains::Int, iterations::Int; subset = 0.8, niter = 1_000, eta = 0.5, xgboostparams...)

    N = length(y)

    # randomly sub-select training and testing set
    Ntrain = round(Int, N*subset)
    Ntest = N - Ntrain

    ids = Random.shuffle(collect(1:N))
    train_ids = ids[1:Ntrain]
    test_ids = ids[(Ntrain+1):end]

    @assert length(test_ids) == Ntest

    # use predicted probabilities?
    mode = iterations > 1 ? "multi:softprob" : "multi:softmax"

    @info "Training classifier"
    # train classifier using XGBoost
    classif = XGBoost.xgboost(x[train_ids,:], niter; label = y[train_ids], 
                      objective = mode, num_class = nchains,
                      xgboostparams...)

    @info "Computing R* statistics"
    Rstats = ones(iterations) * Inf
    for i in 1:iterations
        # predict labels for "test" data
        p = XGBoost.predict(classif, x[test_ids,:])

        pred = if length(p) == Ntest*nchains
            probs = reshape(p, Ntest, nchains)
            map(s -> rand(Categorical(s / sum(s))), eachrow(probs))
        else
            p
        end

        # compute statistic
        a = mean(pred .== y[test_ids])
        Rstats[i] = nchains*a
    end

    return Rstats
end

function rstar(chn::Chains, iterations; kwargs...)
    nchains = size(chn, 3)
    @assert nchains > 1

    # collect data
    x = mapreduce(c -> Array(chn[:,:,c]), vcat, chains(chn))
    y = mapreduce(c -> ones(Int, length(chn))*c, vcat, chains(chn)) .- 1

    return rstar(x, y, nchains, iterations, kwargs...)
end

rstar(chn::Chains; kwargs...) = first(rstar(chn, 1; kwargs...))
