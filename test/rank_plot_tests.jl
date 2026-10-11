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

    @testset "ties do not order the chains by where they sit" begin
        # Every draw is one of two values, so almost every rank is tied. Chains drawn from
        # the same distribution have to land on the same bins, whatever order the array
        # happens to hold them in.
        Random.seed!(11)
        val = float.(rand(0:1, 2000, 4))
        _, counts = MCMCChains._rank_bin_counts(val, 20)
        @test sum(counts) == 8000
        # Breaking ties by position puts the first chain on the lowest ranks and the last
        # on the highest, so the two would occupy different bins entirely.
        @test findall(>(0), counts[:, 1]) == findall(>(0), counts[:, 4])
        @test vec(sum(counts; dims = 1)) == fill(2000, 4)
    end

    @testset "missing draws are skipped" begin
        val = Array{Union{Missing,Float64}}(randn(100, 2))
        val[5] = missing
        _, counts = MCMCChains._rank_bin_counts(val, 10)
        @test sum(counts) == 199
    end

    @testset "recipe" begin
        chn = Chains(randn(niter, 2, nchains), [:a, :b])
        @test rankplot(chn) isa Plots.Plot
        @test plot(chn; seriestype = :rankplot) isa Plots.Plot
        @test :rankplot in MCMCChains.supportedplots
    end

    @testset "each chain keeps its own colour and name" begin
        chn = Chains(randn(niter, 1, 3), [:a])
        p = rankplot(chn)
        steps = filter(s -> s[:seriestype] === :path, p.series_list)
        @test length(steps) == 3

        # The reference line must not take the first colour or the first label from under
        # the chains, which is what emitting it before them did.
        channels(c) =
            round.(
                Float64.((
                    Plots.Colors.red(c),
                    Plots.Colors.green(c),
                    Plots.Colors.blue(c),
                ));
                digits = 3,
            )
        palette = MCMCChains.CHAIN_PALETTE
        for (c, s) in enumerate(steps)
            @test s[:label] == "Chain $c"
            @test channels(s[:linecolor]) ==
                  channels(parse(Plots.Colors.Colorant, palette[c]))
        end
    end
end
