using MCMCChains
using FFTW

using Random
using Statistics
using Test

@testset "ESS per second" begin
    c1 = Chains(randn(100,5, 1), info = (start_time=time(), stop_time = time()+1))
    c2 = Chains(randn(100,5, 1), info = (start_time=time(), stop_time = time()+1))
    c = chainscat(c1, c2)

    wall = MCMCChains.wall_duration(c)
    compute = MCMCChains.compute_duration(c)

    @test wall <= compute
    @test compute == (MCMCChains.compute_duration(c1) + MCMCChains.compute_duration(c2))
    
    s = ess(c)
    @test length(s[:,:ess_per_sec]) == 5
    @test all(map(!ismissing, s[:,:ess_per_sec]))
end

@testset "copy and split" begin
    # check a matrix with even number of rows
    x = rand(50, 20)

    # check incompatible sizes
    @test_throws DimensionMismatch MCMCChains.copyto_split!(similar(x, 25, 20), x)
    @test_throws DimensionMismatch MCMCChains.copyto_split!(similar(x, 50, 40), x)

    y = similar(x, 25, 40)
    MCMCChains.copyto_split!(y, x)
    @test reshape(y, size(x)) == x

    # check a matrix with odd number of rows
    x = rand(51, 20)

    # check incompatible sizes
    @test_throws DimensionMismatch MCMCChains.copyto_split!(similar(x, 25, 20), x)
    @test_throws DimensionMismatch MCMCChains.copyto_split!(similar(x, 51, 40), x)

    MCMCChains.copyto_split!(y, x)
    @test reshape(y, 50, 20) == x[vcat(1:25, 27:51), :]
end

@testset "ESS and R̂ (IID samples)" begin
    Random.seed!(20)

    rawx = randn(10_000, 40, 10)

    # Repeat tests with different scales
    for scale in (1, 50, 100)
        x = scale * rawx

        ess_standard, rhat_standard = MCMCChains.ess_rhat(x)
        ess_standard2, rhat_standard2 = MCMCChains.ess_rhat(x; method = ESSMethod())
        ess_fft, rhat_fft = MCMCChains.ess_rhat(x; method = FFTESSMethod())
        ess_bda, rhat_bda = MCMCChains.ess_rhat(x; method = BDAESSMethod())

        # check that we get (roughly) the same results
        @test ess_standard == ess_standard2
        @test ess_standard ≈ ess_fft
        @test rhat_standard == rhat_standard2 == rhat_fft == rhat_bda

        # check that the estimates are reasonable
        @test all(x -> isapprox(x, 100_000; atol = 5_000), ess_standard)
        @test all(x -> isapprox(x, 100_000; atol = 5_000), ess_bda)
        @test all(x -> isapprox(x, 1; atol = 0.1), rhat_standard)

        # BDA method fluctuates more
        @test var(ess_standard) < var(ess_bda)
    end
end

@testset "ESS and R̂ (identical samples)" begin
    x = ones(10_000, 40, 10)

    ess_standard, rhat_standard = MCMCChains.ess_rhat(x)
    ess_standard2, rhat_standard2 = MCMCChains.ess_rhat(x; method = ESSMethod())
    ess_fft, rhat_fft = MCMCChains.ess_rhat(x; method = FFTESSMethod())
    ess_bda, rhat_bda = MCMCChains.ess_rhat(x; method = BDAESSMethod())

    # check that the estimates are all NaN
    for ess in (ess_standard, ess_standard2, ess_fft, ess_bda)
        @test all(isnan, ess)
    end
    for rhat in (rhat_standard, rhat_standard2, rhat_fft, rhat_bda)
        @test all(isnan, rhat)
    end
end

@testset "ESS and R̂ (chains)" begin
    x = rand(10_000, 40, 10)
    chain = Chains(x)

    for method in (ESSMethod(), FFTESSMethod(), BDAESSMethod())
        # analyze chain
        ess_df = ess(chain; method = method)

        # analyze array
        ess_array, rhat_array = MCMCChains.ess_rhat(x; method = method)

        @test ess_df[:,2] == ess_array
        @test ess_df[:,3] == rhat_array
    end
end

@testset "ESS and R̂ (single sample)" begin # check that issue #137 is fixed
    val = rand(1, 5, 3)
    chain = Chains(val, ["a", "b", "c", "d", "e"])

    for method in (ESSMethod(), FFTESSMethod(), BDAESSMethod())
        # analyze chain
        ess_df = ess(chain; method = method)

        # analyze array
        ess_array, rhat_array = MCMCChains.ess_rhat(val; method = method)

        @test ismissing(ess_df[:,2][1]) # since min(maxlag, niter - 1) = 0
        @test ismissing(ess_df[:,3][1])
    end
end
