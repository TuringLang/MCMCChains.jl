using Test
using MCMCChains
using StatsPlots
using Random

@testset "Histograms" begin
    Random.seed!(42)

    @testset "how many bins" begin
        values = randn(2000)
        @test MCMCChains._bin_count(values, :sturges) == ceil(Int, log2(2000)) + 1

        # Freedman and Diaconis ask for more bins than Sturges on a sample this size, and
        # `:auto` takes the larger of the two.
        fd = MCMCChains._bin_count(values, :fd)
        @test fd > MCMCChains._bin_count(values, :sturges)
        @test MCMCChains._bin_count(values, :auto) == fd

        # More draws of the same thing, more bins.
        @test MCMCChains._bin_count(randn(20000), :fd) >
              MCMCChains._bin_count(randn(200), :fd)

        # Nothing to bin on, so one bin rather than an error or a division by zero.
        @test MCMCChains._bin_count(fill(2.0, 100), :auto) == 1
        @test MCMCChains._bin_count([1.0], :auto) == 1
        # A sample with no spread between its quartiles still gets Sturges.
        spiked = vcat(fill(0.0, 999), 10.0)
        @test MCMCChains._bin_count(spiked, :auto) ==
              MCMCChains._bin_count(spiked, :sturges)

        @test_throws ArgumentError MCMCChains._bin_count(values, :nonsense)
    end

    @testset "every chain is binned the same way" begin
        # Two chains that do not overlap, so binning them apart would be obvious.
        series = [randn(500), randn(500) .+ 8]
        edges = MCMCChains._histogram_edges(series, 20)

        # The edges span both chains rather than either one of them.
        @test first(edges) ≈ minimum(Iterators.flatten(series))
        @test last(edges) > maximum(Iterators.flatten(series))
        # One bin past the top, so the largest draw has somewhere to go.
        @test length(edges) == 22
        @test step(edges) ≈ (maximum(Iterators.flatten(series)) - first(edges)) / 20
    end

    @testset "a chain that never moved still has bins" begin
        edges = MCMCChains._histogram_edges([[2.0, 2.0, 2.0]], 5)
        @test first(edges) < 2.0 < last(edges)
    end

    @testset "what the caller asks for wins" begin
        edges = range(-1, 1; length = 5)
        @test MCMCChains._histogram_edges([randn(10)], edges) === edges

        series = [randn(1000), randn(1000)]
        pooled = vcat(series...)
        auto = MCMCChains._histogram_edges(series, :auto)
        @test length(auto) == MCMCChains._bin_count(pooled, :auto) + 2
        @test length(MCMCChains._histogram_edges(series, :sturges)) <= length(auto)
    end

    @testset "recipe" begin
        chn = Chains(randn(400, 2, 3), [:a, :b])
        p = histogram(chn)
        @test p isa Plots.Plot

        # Every chain of a parameter is drawn over one set of edges.
        edges = [s.plotattributes[:bins] for s in p.series_list]
        @test length(unique(edges)) == length(names(chn))
        @test all(e -> e isa AbstractRange, edges)

        @test histogram(chn, bins = 25) isa Plots.Plot
        @test histogram(Chains(randn(400, 1, 1), [:a])) isa Plots.Plot
        # A discrete parameter reaches the same path from mixeddensity.
        @test mixeddensity(Chains(float.(rand(1:2, 400, 1, 3)), [:d])) isa Plots.Plot
    end
end
