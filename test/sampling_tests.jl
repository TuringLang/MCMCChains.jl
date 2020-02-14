using MCMCChains
using Distributions
using KernelDensity
using StatsBase

using Random
using Statistics
using Test

Random.seed!(20)

@testset "sampling api" begin
    # Sample some data.
    vals = hcat(rand(Normal(1, 0.01), 500), rand(Normal(5, 0.01), 500))
    chn = Chains(vals, ["m", "s"])
    
    chn_sample = sample(chn, 5)
    @test range(chn_sample) == 1:1:5

    c = kde(Array(chn[:s]))
    chn_weighted_sample = sample(c.x, Weights(c.density), 100000)
    @test mean(chn_weighted_sample) ≈ 5.0 atol=0.1
end
