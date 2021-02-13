using MCMCChains
using Distributions
using MLJModels
using Test

N = 1000
val = rand(N, 8, 4)
colnames = ["a", "b", "c", "d", "e", "f", "g", "h"]
internal_colnames = ["c", "d", "e", "f", "g", "h"]
chn = Chains(val, colnames, Dict(:internals => internal_colnames))
X = randn(N, 2)
y = rand(1:3, N)

XGBoost = @load XGBoostClassifier
ThresholdPredictor = @load BinaryThresholdPredictor

@testset "R star test" begin
    for classif in (XGBoost(), ThresholdPredictor(wrapped_model=XGBoost()))
        # Compute R* statistic for a mixed chain.
        R = rstar(classif, X, y)

        # Resulting R value should be close to one, i.e. the classifier does not perform better than random guessing.
        if classif isa MCMCChains.MMI.Probabilistic
            @test R isa DiscreteNonParametric
            @test extrema(R) == (0, 3)
        else
            @test R isa Dirac
        end
        @test mean(R) ≈ 1 atol=0.15

        # Compute R* statistic for a mixed chain.
        R = rstar(classif, chn)

        # Resulting R value should be close to one, i.e. the classifier does not perform better than random guessing.
        if classif isa MCMCChains.MMI.Probabilistic
            @test R isa DiscreteNonParametric
            @test extrema(R) == (0, size(chn, 3))
        else
            @test R isa Dirac
        end
        @test mean(R) ≈ 1 atol=0.15

        # Compute R* statistic for a non-mixed chain.
        niter = 1000
        val = hcat(sin.(1:niter), cos.(1:niter))
        val = cat(val, hcat(cos.(1:niter)*100, sin.(1:niter)*100), dims=3)
        chn_notmixed = Chains(val)

        # Restuling R value should be close to two, i.e. the classifier should be able to learn an almost perfect decision boundary between chains.
        R = rstar(classif, chn_notmixed)
        if classif isa MCMCChains.MMI.Probabilistic
            @test R isa DiscreteNonParametric
            @test extrema(R) == (0, size(chn_notmixed, 3))
        else
            @test R isa Dirac
        end
        @test mean(R) ≈ 2 atol=0.1
    end
end
