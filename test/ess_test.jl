using MCMCChains
using Distributions
using Random
using Test

Random.seed!(20)

@testset "ess tests" begin
    # Consider the model
    # ```math
    # θ ~ Beta(1, 1)
    # y ~ Bernoulli(θ)
    # ```
    # for variable ``y`` of successes and failures. Given 10 observations of 2 successes
    # and 8 failures, sample four chains with 1000 samples each from the posterior
    # ```math
    # θ ~ Beta(3, 9)
    # ```
    chain = Chains(rand(Beta(3, 9), 1000, 1, 4))

    # Compute the expected sample size
    ess_chain = ess(chain)

    @test all(x -> 1000 ≤ x ≤ 2000, ess_chain[:ess])
    @test all(x -> isapprox(x , 1.0, atol=0.1), ess_chain[:r_hat])
end
