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

    c = kde(vec(chn[:s]))
    chn_weighted_sample = sample(c.x, Weights(c.density), 100000)
    @test mean(chn_weighted_sample) ≈ 5.0 atol=0.1

    # issue #223
    @testset "multiple chains" begin
	chn = Chains(randn(11, 4, 3))
	wv = pweights(rand(33))

	for kwargs in ((), (replace=true,), (replace=false,), (ordered=false,), (ordered=true,))
            # without weights
	    @test sample(chn, 5; kwargs...) isa Chains
	    @test size(sample(chn, 5; kwargs...)) == (5, 4, 1)

            # with weights
	    @test sample(chn, wv, 5; kwargs...) isa Chains
	    @test size(sample(chn, wv, 5; kwargs...)) == (5, 4, 1)
	end
    end
end
