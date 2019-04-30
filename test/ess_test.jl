using MCMCChains
using Turing
using Random
using Test

@testset "ess tests" begin
    Random.seed!(20)

    dat = [0,1,0,0,0,0,0,0,0,1]

    @model bermodel(y) = begin
        theta ~ Beta(1,1)
        for n = 1:length(y)
            y[n] ~ Bernoulli(theta)
        end
    end

    # Make four chains.
    chn = mapreduce(x->sample(bermodel(dat),
        NUTS(2000, 1000, 0.65)),
        chainscat,
        1:4)

    e = ess(chn[1001:2000, :, :])
    display(e)

    @test all([z ≥ 1000 && z ≤ 2000 for z in e[:ess]])
    @test all(isapprox.(e[:r_hat], 1.0, atol=0.1))
end
