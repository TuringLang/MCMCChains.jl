using MCMCChains
using Distributions
using Random
using Test

Random.seed!(20)

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
    vals =
        hcat(rand(MvNormal(1:10, 0.1), n_samples)', rand(InverseGamma(2.5, 5), n_samples))
    chn = Chains(vals, [("m[$i]" for i = 1:10)..., "s"])

    get1 = get(chn, :m)
    get2 = get(chn, [:m, :s])

    @test length(get1.m) == 10
    @test length(get1.m[5]) == n_samples
    @test round(mean(get1.m[1])) ≈ 1.0
    @test round(mean(get1.m[10])) ≈ 10.0
    @test length(get2.s) == n_samples
    @test sort(collect(keys(get2))) == [:m, :s]

    getall = get_params(chn)
    @test getall == get(chn, section = [:parameters])
    @test length(keys(getall)) == 2
end

@testset "flattened sections" begin
    chn = Chains(
        reshape(1:24, 4, 3, 2),
        [:a, Symbol("x[1]"), Symbol("x[2]")],
        Dict(:parameters => [:a, Symbol("x[1]")], :internals => [Symbol("x[2]")]),
    )
    for (section, selected) in (
        (:parameters, [:a, Symbol("x[1]")]),
        ([:internals, :parameters], [Symbol("x[2]"), :a, Symbol("x[1]")]),
        ([:parameters, :parameters], [:a, Symbol("x[1]")]),
    )
        result = get(chn; section, flatten = true)
        @test keys(result) == Tuple(selected)
        @test result == get(chn, selected; flatten = true)
    end
    @test get_params(chn; flatten = true) == get(chn, names(chn); flatten = true)
    @test get(chn; section = Symbol[], flatten = true) == NamedTuple()
end
