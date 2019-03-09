using MCMCChains
using Turing
using Test

@testset "get tests" begin
    @model model(x) = begin
        m = Vector{Real}(undef, 10)
        [m[i] ~ Normal(0, 0.5) for i in 1:10]
        s ~ InverseGamma(2,3)
        x ~ Normal(0, s)
    end

    model = model(2.0)
    sampler = HMC(500, 0.01, 5)
    chn = sample(model, sampler)

    get1 = get(chn, :m)
    get2 = get(chn, [:m, :s])
    
    @test length(get1.m) == 10
    @test length(get1.m[5]) == 500
    @test length(get2.s) == 500
    @test collect(keys(get2)) == [:m, :s]
end
