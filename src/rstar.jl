# Data structure required to store samples for MLJ.
struct RStarTable{T<:AbstractFloat}
    data::Matrix{T}
    names::Vector{Symbol}
end
RStarTable(chn::Chains) = RStarTable(Array(chn), names(chn, :parameters))

Base.size(table::RStarTable) = size(table.data)
Base.size(table::RStarTable, i::Integer) = size(table.data, i)
Base.getindex(table::RStarTable, i::Integer) = RStarTable(table.data[i,:], table.names)
Base.getindex(table::RStarTable, i::AbstractVector{<:Integer}) = RStarTable(table.data[i,:], table.names)

Tables.istable(::Type{<:RStarTable}) = true

Tables.names(table::RStarTable) = table.names
Tables.matrix(table::RStarTable) = table.data
Tables.schema(table::RStarTable{T}) where {T} = Tables.Schema(table.names, fill(eltype(T), size(table.data,2)))

Tables.columnaccess(::Type{<:RStarTable}) = true
Tables.columns(table::RStarTable) = table

Tables.columnnames(table::RStarTable) = Tables.names(table)
Tables.getcolumn(table::RStarTable, nm::Symbol) = table.data[:,findfirst(table.names.==nm)]
Tables.getcolumn(table::RStarTable, i::Int) = table.data[:,i]

"""
    rstar(chains::Chains, model::MLJModelInterface.Supervised; subset = 0.8, iterations = 10, verbosity = 0)
    rstar(rng::Random.AbstractRNG, chains::Chains, model::MLJModelInterface.Supervised; subset = 0.8, iterations = 10, verbosity = 0)
    rstar(rng::Random.AbstractRNG, model::MLJModelInterface.Supervised, x::MCMCChains.RStarTable, y::AbstractVector; subset = 0.8, iterations = 10, verbosity = 0)

Compute the R* statistic for convergence diagnostic of MCMC. This implementation is an adaption of Algorithm 1 & 2, described in [Lambert & Vehtari]. Note that the correctness of the statistic depends on the convergence of the classifier used internally in the statistic. You can track if the training of the classifier converged by inspection of the printed RMSE values from the XGBoost backend. To adjust the number of iterations used to train the classifier set `niter` accordingly.

# Keyword Arguments
* `subset = 0.8` ... Subset used to train the classifier, i.e. 0.8 implies 80% of the samples are used.
* `iterations = 10` ... Number of iterations used to estimate the statistic. If the classifier is not probabilistic, i.e. does not return class probabilities, it is advisable to use a value of one.
* `verbosity = 0` ... Verbosity level used during fitting of the classifier.

# Usage

```julia-repl
using MLJ, MLJModels
# You need to load MLJBase and the respective package your are using for classification first.

# Select a classification model to compute the Rstar statistic.
# For example the XGBoost classifier.
model = @load XGBoostClassifier()

# Compute 100 samples of the R* statistic using sampling from according to the prediction probabilities.
Rs = rstar(chn, model, iterations = 20)

# estimate Rstar
R = mean(Rs)

# visualize distribution
histogram(Rs)
```

## References:
[Lambert & Vehtari] Ben Lambert and Aki Vehtari. "R∗: A robust MCMC convergence diagnostic with uncertainty using gradient-boostined machines." Arxiv 2020.
"""
function rstar(rng::Random.AbstractRNG, classif::MLJModelInterface.Supervised, x::RStarTable, y::AbstractVector{Int}; iterations = 10, subset = 0.8, verbosity = 0)

    size(x,1) != length(y) && throw(DimensionMismatch())
    iterations >= 1 && ArgumentError("Number of iterations has to be positive!")

    if iterations > 1 && classif isa Deterministic
        @warn("Classifier is not a probabilistic classifier but number of iterations is > 1.")
    elseif iterations == 1 && classif isa Probabilistic
        @warn("Classifier is probabilistic but number of iterations is equal to one.")
    end

    N = length(y)
    K = length(unique(y))

    # randomly sub-select training and testing set
    Ntrain = round(Int, N*subset)
    Ntest = N - Ntrain

    ids = Random.randperm(rng, N)
    train_ids = view(ids, 1:Ntrain)
    test_ids = view(ids, (Ntrain+1):N)

    # train classifier using XGBoost
    fitresult, _ = MLJModelInterface.fit(classif, verbosity, x[train_ids], MLJModelInterface.categorical(y[train_ids]))

    Rstats = fill(Inf, iterations)
    xtest = x[test_ids]
    ytest = MLJModelInterface.categorical(view(y, test_ids))

    for i in 1:iterations

        # predict labels for "test" data
        p = MLJModelInterface.predict(classif, fitresult, xtest)
        pred = similar(ytest)

        if classif isa Deterministic
            pred[:] = p
        elseif classif isa Probabilistic
            pred[:] = get.(rand.(rng, p))
        else
            throw(ErrorException("Unknown type of classifier"))
        end

        # compute statistic
        a = mean(((p,y),) -> p == y, zip(pred, ytest))
        Rstats[i] = K*a
    end

    return Rstats
end

function rstar(chn::Chains, model::MLJModelInterface.Supervised; kwargs...)
    return rstar(Random.default_rng(), chn, model; kwargs...)
end

function rstar(rng::Random.AbstractRNG, chn::Chains, model::MLJModelInterface.Supervised; kwargs...)
    nchains = size(chn, 3)
    nchains <= 1 && throw(DimensionMismatch())

    # collect data
    x = RStarTable(chn)
    y = repeat(chains(chn); inner = size(chn,1))

    return rstar(rng, model, x, y; kwargs...)
end
