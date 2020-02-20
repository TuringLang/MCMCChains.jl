println("imports")
@time using MCMCChains
@time using Distributions
@time using Random
@time using Test

@time Random.seed!(20)

@testset "get tests" begin
    # Consider the model
    # ```math
    # m ~ Normal([1, 2, ..., 10], 0.1^2)
    # s ~ IG(2, 3)
    # x ~ Normal(0, s)
    # ```
    # for variable ``x``. Given observation ``x = 2``, sample one chain with 1000 samples
    # from the posterior
    # ```math
    # m ~ Normal([1, 2, ..., 10], 0.1^2)
    # s ∼ IG(2.5, 5)
    # ```
    n_samples = 1000
    println("sampling")
    @time vals = hcat(rand(MvNormal(1:10, 0.1), n_samples)',
                rand(InverseGamma(2.5, 5), n_samples))
    println("constructor")
    @time chn = Chains(vals, [("m[$i]" for i in 1:10)..., "s"])

    println("get")
    @time get1 = get(chn, :m)
    @time get2 = get(chn, [:m, :s])

    @test length(get1.m) == 10
    @test length(get1.m[5]) == n_samples
    @test round(mean(get1.m[1])) ≈ 1.0
    @test round(mean(get1.m[10])) ≈ 10.0
    @test length(get2.s) == n_samples
    @test collect(keys(get2)) == [:m, :s]

    println("getall")
    @time getall = get_params(chn)
    @time @test getall == get(chn, section = [:parameters])
    @test length(keys(getall)) == 2
end
