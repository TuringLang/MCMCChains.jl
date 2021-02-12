using MCMCChains
using MLJModels
using Test

N = 1000
val = rand(N, 8, 4)
colnames = ["a", "b", "c", "d", "e", "f", "g", "h"]
internal_colnames = ["c", "d", "e", "f", "g", "h"]
chn = Chains(val, colnames, Dict(:internals => internal_colnames))

XGBoost = @load XGBoostClassifier
ThresholdPredictor = @load BinaryThresholdPredictor

@testset "R star test" begin
    @testset "Probabilistic classifier" begin
        classif = XGBoost()

        # Compute R* statistic for a mixed chain.
        R = rstar(classif, randn(N,2), rand(1:3,N))

        # Resulting R value should be close to one, i.e. the classifier does not perform better than random guessing.
        @test mean(R) ≈ 1 atol=0.15

        # Compute R* statistic for a mixed chain.
        R = rstar(classif, chn)

        # Resulting R value should be close to one, i.e. the classifier does not perform better than random guessing.
        @test mean(R) ≈ 1 atol=0.15

        # Compute R* statistic for a non-mixed chain.
        niter = 1000
        val = hcat(sin.(1:niter), cos.(1:niter))
        val = cat(val, hcat(cos.(1:niter)*100, sin.(1:niter)*100), dims=3)
        chn_notmixed = Chains(val)

        # Restuling R value should be close to two, i.e. the classifier should be able to learn an almost perfect decision boundary between chains.
        R = rstar(classif, chn_notmixed)
        @test mean(R) ≈ 2 atol=0.1
    end

    @testset "Deterministic classifier" begin
        classif = BinaryThresholdPredictor(wrapped_model=XGBoost())

        # Compute R* statistic for a mixed chain.
        R = (
            @test_logs (
                :warn,
                "Classifier is not a probabilistic classifier but number of iterations is > 1.",
            ) match_mode=:any rstar(classif, randn(N,2), rand(1:3,N))
        )

        # Resulting R value should be close to one, i.e. the classifier does not perform better than random guessing.
        @test all(r == first(R) for r in R)
        @test first(R) ≈ 1 atol=0.15

        # Compute R* statistic for a mixed chain.
        R = (
            @test_logs (
                :warn,
                "Classifier is not a probabilistic classifier but number of iterations is > 1.",
            ) match_mode=:any rstar(classif, chn)
        )

        # Resulting R value should be close to one, i.e. the classifier does not perform better than random guessing.
        @test all(r == first(R) for r in R)
        @test first(R) ≈ 1 atol=0.15

        # Compute R* statistic for a non-mixed chain.
        niter = 1000
        val = hcat(sin.(1:niter), cos.(1:niter))
        val = cat(val, hcat(cos.(1:niter)*100, sin.(1:niter)*100), dims=3)
        chn_notmixed = Chains(val)

        # Restuling R value should be close to two, i.e. the classifier should be able to learn an almost perfect decision boundary between chains.
        R = (
            @test_logs (
                :warn,
                "Classifier is not a probabilistic classifier but number of iterations is > 1.",
            ) match_mode=:any rstar(classif, chn_notmixed)
        )
        @test all(r == first(R) for r in R)
        @test first(R) ≈ 2 atol=0.1
    end
end
