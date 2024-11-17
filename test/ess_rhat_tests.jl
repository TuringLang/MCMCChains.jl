using MCMCChains
using FFTW

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
        @test length(s[:ess_per_sec]) == 5
        @test all(map(!ismissing, s[:ess_per_sec]))
    end
end

@testset "ess/rhat/ess_rhat (chains)" begin
    x = rand(10_000, 40, 10)
    chain = Chains(x)

    for autocov_method in (AutocovMethod(), FFTAutocovMethod(), BDAAutocovMethod()), kind in (:bulk, :basic), f in (ess, ess_rhat, rhat)
        # analyze chain
        ess_df = ess(chain; autocov_method = autocov_method, kind = kind)
        rhat_df = rhat(chain; kind = kind)
        ess_rhat_df = ess_rhat(chain; autocov_method = autocov_method, kind = kind)

        # analyze array
        ess_array, rhat_array = ess_rhat(
            permutedims(x, (1, 3, 2)); autocov_method = autocov_method, kind = kind,
        )
        @test ess_df[:ess] == ess_rhat_df[:ess] == ess_array
        @test rhat_df[:rhat] == ess_rhat_df[:rhat] == rhat_array
    end
end

@testset "ESS and R̂ (single sample)" begin # check that issue #137 is fixed
    val = rand(1, 5, 3)
    chain = Chains(val, ["a", "b", "c", "d", "e"])

    for autocov_method in (AutocovMethod(), FFTAutocovMethod(), BDAAutocovMethod())
        # analyze chain
        ess_df = ess(chain; autocov_method = autocov_method)
        @test isequal(ess_df[:ess], fill(NaN, 5))
        @test isequal(ess_df[:ess_per_sec], fill(missing, 5))
        
        ess_rhat_df = ess_rhat(chain; autocov_method = autocov_method)
        @test isequal(ess_rhat_df[:ess], fill(NaN, 5))
        @test isequal(ess_rhat_df[:rhat], fill(NaN, 5))
        @test isequal(ess_rhat_df[:ess_per_sec], fill(missing, 5))
    end

    rhat_df = rhat(chain)
    @test isequal(rhat_df[:rhat], fill(NaN, 5))
end
