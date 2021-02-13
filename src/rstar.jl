"""
    rstar([rng ,]classif::Supervised, chains::Chains; kwargs...)
    rstar([rng ,]classif::Supervised, x::AbstractMatrix, y::AbstractVector; kwargs...)

Compute the distribution of the R* convergence diagnostic of MCMC.

This implementation is an adaption of Algorithm 1 & 2, described in [^LambertVehtari]. In
contrast to Algorithm 2 it does return the analytical distribution of the statistic
instead of empirical samples.

Note that the correctness of the statistic depends on the convergence of the classifier used
internally in the statistic. You can inspect the training of the classifier by adjusting the
verbosity level.

# Keyword Arguments

- `subset = 0.8`: Subset used for training, i.e. 0.8 implies 80% of the samples are used.
- `verbosity = 0`: Verbosity level used during fitting of the classifier.

# Usage

```julia
# Load an MLJ classifier to compute the statistic, e.g., the XGBoost classifier.
using MLJModels
XGBoost = @load XGBoostClassifier

# Compute the distribution of the R* statistic.
rstar(XGBoost(), chn)
```

[^LambertVehtari]: Ben Lambert and Aki Vehtari. "R∗: A robust MCMC convergence diagnostic with uncertainty using gradient-boostined machines." Arxiv 2020.
"""
function rstar(
    rng::Random.AbstractRNG, classif::MMI.Supervised, x::AbstractMatrix, y::AbstractVector{Int};
    subset=0.8, verbosity::Int=0,
)
    N = length(y)
    size(x, 1) != N && throw(DimensionMismatch())

    # randomly sub-select training and testing set
    Ntrain = round(Int, N * subset)
    Ntest = N - Ntrain
    ids = Random.randperm(rng, N)
    train_ids = view(ids, 1:Ntrain)
    test_ids = view(ids, (Ntrain+1):N)

    # train classifier
    y_categorical = MMI.categorical(y)
    fitresult, _ = MMI.fit(
        classif,
        verbosity,
        Tables.table(x[train_ids, :]),
        y_categorical[train_ids],
    )

    # compute predictions on test data
    xtest = Tables.table(x[test_ids, :])
    pred = MMI.predict(classif, fitresult, xtest)

    # compute distribution of statistic
    ytest = y_categorical[test_ids]
    dist = rstar_distribution(pred, ytest)

    return dist
end

function rstar(
    classif::MMI.Supervised, x::AbstractMatrix, y::AbstractVector{Int}; kwargs...
)
    return rstar(Random.GLOBAL_RNG, classif, x, y; kwargs...)
end

function rstar(classif::MMI.Supervised, chn::Chains; kwargs...)
    return rstar(Random.GLOBAL_RNG, classif, chn; kwargs...)
end

function rstar(
    rng::Random.AbstractRNG, classif::MMI.Supervised, chn::Chains; kwargs...
)
    nchains = size(chn, 3)
    nchains <= 1 && throw(DimensionMismatch())

    # collect data
    x = Array(chn)
    y = repeat(chains(chn); inner = size(chn,1))

    return rstar(rng, classif, x, y; kwargs...)
end

function rstar_distribution(pred::T, ytest::T) where {T<:AbstractVector}
    mean_accuracy = mean(zip(pred, ytest)) do (p, y)
        p == y
    end
    nclasses = length(MMI.classes(ytest))
    return Dirac(nclasses * mean_accuracy)
end

function rstar_distribution(
    pred::AbstractVector{<:UnivariateDistribution},
    ytest::AbstractVector,
)
    # Compute probabilities of Poisson Binomial distribution with support `0:length(ytest)`.
    probs = Distributions.poissonbinomial_pdf(map(pdf, pred, ytest))
    nclasses = length(MMI.classes(ytest))
    return DiscreteNonParametric(
        range(0, nclasses; length=length(ytest) + 1), probs,
    )
end
