using MCMCChains
using Turing
using Test

@testset "get tests" begin
    @model model(x) = begin
        m = Vector{Real}(undef, 10)
        [m[i] ~ Normal(i, 0.1) for i in 1:10]
        s ~ InverseGamma(2,3)
        x ~ Normal(0, s)
    end

    n_samples = 1000
    model = model(2.0)
    sampler = MH(n_samples)
    chn = sample(model, sampler)

    get1 = get(chn, :m)
    get2 = get(chn, [:m, :s])

    @test length(get1.m) == 10
    @test length(get1.m[5]) == n_samples
    @test round(mean(get1.m[1])) ≈ 1.0
    @test round(mean(get1.m[10])) ≈ 10.0
    @test length(get2.s) == n_samples
    @test collect(keys(get2)) == [:m, :s]

    getsection = get(chn, section = :internals)
    @test length(keys(getsection)) == 2

    getall = get_params(chn)
    @test getall == get(chn, section=[:internals, :parameters])
    @test length(keys(getall)) == 4
end
