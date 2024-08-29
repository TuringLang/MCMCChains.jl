"""
    rstar(rng=Random.GLOBAL_RNG, classifier, chains::Chains; kwargs...)

Compute the ``R^*`` convergence diagnostic of the MCMC `chains` with the `classifier`.

The keyword arguments supported here are the same as those in `rstar` for arrays of samples
and chain indices.

# Examples

```jldoctest rstar; setup = :(using Random; Random.seed!(200))
julia> using MLJBase, MLJDecisionTreeInterface, Statistics

julia> chains = Chains(fill(4.0, 100, 2, 3));
```

One can compute the distribution of the ``R^*`` statistic with the probabilistic classifier.

```jldoctest rstar
julia> distribution = rstar(DecisionTreeClassifier(), chains);

julia> isapprox(mean(distribution), 1; atol=0.1)
true
```

For deterministic classifiers, a single ``R^*`` statistic is returned.

```jldoctest rstar
julia> decisiontree_deterministic = Pipeline(
           DecisionTreeClassifier();
           operation=predict_mode,
       );

julia> value = rstar(decisiontree_deterministic, chains);

julia> isapprox(value, 1; atol=0.2)
true
```
"""
function MCMCDiagnosticTools.rstar(
    classif, chn::Chains; kwargs...
)
    return MCMCDiagnosticTools.rstar(Random.GLOBAL_RNG, classif, chn; kwargs...)
end

function MCMCDiagnosticTools.rstar(
    rng::Random.AbstractRNG,
    classif,
    chn::Chains;
    sections = _default_sections(chn),
    kwargs...
)
    nchains = size(chn, 3)
    nchains <= 1 && throw(DimensionMismatch())

    _chn = Chains(chn, _clean_sections(chn, sections))
    x = _permutedims_diagnostics(_chn.value.data)

    return MCMCDiagnosticTools.rstar(rng, classif, x; kwargs...)
end
