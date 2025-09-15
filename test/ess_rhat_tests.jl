using MCMCChains
using DataFrames
using FFTW
using PosteriorStats: SummaryStats

using Random
using Statistics
using Test

@testset "ESS per second" begin
    c1 = Chains(randn(100,5, 1), info = (start_time=time(), stop_time = time()+1))
    c2 = Chains(randn(100,5, 1), info = (start_time=time()+1, stop_time = time()+2))
    c = chainscat(c1, c2)

    wall = MCMCChains.wall_duration(c)
    compute = MCMCChains.compute_duration(c)

    @test round(wall, digits=1) ≈ round(c2.info.stop_time - c1.info.start_time, digits=1)
    @test compute ≈ (MCMCChains.compute_duration(c1) + MCMCChains.compute_duration(c2))

    for f in (ess, ess_rhat)
        s = f(c)
        @test s isa SummaryStats
        df = DataFrame(s)
        @test length(df[!, :ess_per_sec]) == 5
        @test all(map(!ismissing, df[!, :ess_per_sec]))
    end
end

@testset "ess/rhat/ess_rhat (chains)" begin
    x = rand(10_000, 40, 10)
    chain = Chains(x)

    for autocov_method in (AutocovMethod(), FFTAutocovMethod(), BDAAutocovMethod()), kind in (:bulk, :basic), f in (ess, ess_rhat, rhat)
        # analyze chain
        ess_stats = ess(chain; autocov_method = autocov_method, kind = kind)
        rhat_stats = rhat(chain; kind = kind)
        ess_rhat_stats = ess_rhat(chain; autocov_method = autocov_method, kind = kind)

        @test ess_stats isa SummaryStats
        @test ess_stats.name == "ESS"
        @test rhat_stats isa SummaryStats
        @test rhat_stats.name == "R-hat"
        @test ess_rhat_stats isa SummaryStats
        @test ess_rhat_stats.name == "ESS/R-hat"

        ess_df, rhat_df, ess_rhat_df = DataFrame.((ess_stats, rhat_stats, ess_rhat_stats))

        # analyze array
        ess_array, rhat_array = ess_rhat(
            permutedims(x, (1, 3, 2)); autocov_method = autocov_method, kind = kind,
        )
        @test ess_df[!, :ess] == ess_rhat_df[!, :ess] == ess_array
        @test rhat_df[!, :rhat] == ess_rhat_df[!, :rhat] == rhat_array
    end
end

@testset "ESS and R̂ (single sample)" begin # check that issue #137 is fixed
    val = rand(1, 5, 3)
    chain = Chains(val, ["a", "b", "c", "d", "e"])

    for autocov_method in (AutocovMethod(), FFTAutocovMethod(), BDAAutocovMethod())
        # analyze chain
        ess_stats = ess(chain; autocov_method = autocov_method)
        @test ess_stats isa SummaryStats
        @test ess_stats.name == "ESS"
        ess_df = DataFrame(ess_stats)
        @test isequal(ess_df[!, :ess], fill(NaN, 5))
        @test isequal(ess_df[!, :ess_per_sec], fill(missing, 5))
        
        ess_rhat_stats = ess_rhat(chain; autocov_method = autocov_method)
        @test ess_rhat_stats isa SummaryStats
        @test ess_rhat_stats.name == "ESS/R-hat"
        ess_rhat_df = DataFrame(ess_rhat_stats)
        @test isequal(ess_rhat_df[!, :ess], fill(NaN, 5))
        @test isequal(ess_rhat_df[!, :rhat], fill(NaN, 5))
        @test isequal(ess_rhat_df[!, :ess_per_sec], fill(missing, 5))
    end

    rhat_stats = rhat(chain)
    @test rhat_stats isa SummaryStats
    @test rhat_stats.name == "R-hat"
    rhat_df = DataFrame(rhat_stats)
    @test isequal(rhat_df[!, :rhat], fill(NaN, 5))
end
