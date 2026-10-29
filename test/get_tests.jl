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

@testset "ordered array parameter access" begin
    chn = Chains(reshape(1:39, 3, 13, 1), ["x[$i]" for i = 13:-1:1])
    components = get(chn, :x).x
    @test length(components) == 13
    @test components[1] == chn[Symbol("x[1]")]
    @test components[2] == chn[Symbol("x[2]")]
    @test components[13] == chn[Symbol("x[13]")]

    requested = [Symbol("x[1]"), Symbol("x[13]")]
    @test names(chn[requested]) == requested
    @test chn[requested].value.data == chn.value.data[:, [13, 1], :]
    @test names(chn[:, requested, :]) == requested
    @test names(chn[String.(requested)]) == requested
    @test names(chn[[requested[1], requested[1]]]) == requested[[1, 1]]
    @test_throws ArgumentError chn[[:missing]]
    scalar_chn = Chains(reshape(1:9, 3, 3, 1), [:a, :b, :c])
    selected = scalar_chn[[:c, :a]]
    @test selected.name_map.parameters == [:c, :a]
    @test collect(keys(get_params(selected))) == [:c, :a]

    matrix_chn = Chains(reshape(1:4, 1, 4, 1), ["θ[2,2]", "θ[1,2]", "θ[2,1]", "θ[1,1]"])
    @test map(only, get(matrix_chn, :θ).θ) == (4, 3, 2, 1)
    nonnumeric_chn = Chains(reshape(1:2, 1, 2, 1), ["x[b]", "x[a]"])
    @test map(only, get(nonnumeric_chn, :x).x) == (1, 2)
    mixed_rank_chn = Chains(reshape(1:2, 1, 2, 1), ["x[1,1]", "x[2]"])
    @test map(only, get(mixed_rank_chn, :x).x) == (1, 2)
    large_index_chn = Chains(reshape(1:2, 1, 2, 1), ["x[999999999999999999999]", "x[1]"])
    @test map(only, get(large_index_chn, :x).x) == (1, 2)

    summary = mean(chn)
    @test summary.nt.mean == [mean(chn[name]) for name in summary.nt.parameters]
end
