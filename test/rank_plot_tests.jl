using Test
using MCMCChains
using StatsPlots
using Random

@testset "Rank plots" begin
    Random.seed!(42)
    niter, nchains, nbins = 500, 3, 20

    @testset "_rank_bin_counts" begin
        val = randn(niter, nchains)
        edges, counts = MCMCChains._rank_bin_counts(val, nbins)

        @test length(edges) == nbins + 1
        @test size(counts) == (nbins, nchains)
        # Every draw is binned exactly once.
        @test sum(counts) == niter * nchains
        @test all(sum(counts; dims = 1) .== niter)

        # Chains from one distribution should sit near the uniform count. The bound is loose
        # because this is a random draw, but a chain that occupied one end of the ranks
        # would blow well past it.
        expected = niter / nbins
        @test maximum(abs.(counts .- expected)) < expected
    end

    @testset "a chain that has not mixed is visible" begin
        # Third chain shifted, so it should take the top ranks and vacate the bottom ones.
        val = hcat(randn(niter, 2), randn(niter) .+ 5)
        _, counts = MCMCChains._rank_bin_counts(val, nbins)

        @test counts[1, 3] == 0                       # none of the lowest ranks
        @test counts[end, 3] > niter / nbins          # more than its share of the highest
        @test counts[end, 3] > counts[end, 1]
    end

    @testset "ties do not lose draws" begin
        val = fill(1.0, 10, 2)
        _, counts = MCMCChains._rank_bin_counts(val, 5)
        @test sum(counts) == 20
    end

    @testset "recipe" begin
        chn = Chains(randn(niter, 2, nchains), [:a, :b])
        @test rankplot(chn) isa Plots.Plot
        @test plot(chn; seriestype = :rankplot) isa Plots.Plot
        @test :rankplot in MCMCChains.supportedplots
    end
end
