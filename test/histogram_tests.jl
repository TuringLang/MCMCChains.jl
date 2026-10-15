using Test
using MCMCChains
using StatsPlots
using Random

@testset "Histograms" begin
    Random.seed!(42)

    @testset "how many bins" begin
        # Every rule has to give the count Plots gives, so that passing `bins` keeps
        # meaning what it meant before these edges were worked out here.
        for values in (randn(100), randn(5000), 3 .* randn(731) .+ 2, rand(50))
            for rule in MCMCChains.BIN_RULES
                @test MCMCChains._bin_count(values, rule) ==
                      Plots._auto_binning_nbins((values,), 1; mode = rule)
            end
        end

        values = randn(2000)
        @test MCMCChains._bin_count(values, :auto) == MCMCChains._bin_count(values, :fd)
        # More draws of the same thing, more bins.
        @test MCMCChains._bin_count(randn(20000), :fd) >
              MCMCChains._bin_count(randn(200), :fd)

        # Nothing to bin on, so one bin rather than an error or a division by zero.
        @test MCMCChains._bin_count(fill(2.0, 100), :auto) == 1
        @test MCMCChains._bin_count([1.0], :auto) == 1
        # Quartiles that coincide, so the width would be a division by zero.
        spiked = vcat(fill(0.0, 999), 10.0)
        for rule in MCMCChains.BIN_RULES
            @test MCMCChains._bin_count(spiked, rule) >= 1
        end

        @test_throws ArgumentError MCMCChains._bin_count(values, :nonsense)
    end

    @testset "a parameter on a few values is binned on them" begin
        # One bin for each value, rather than a width read off the spread.
        @test MCMCChains._discrete_edges(float.(rand(1:2, 400))) == [0.5, 1.5, 2.5]
        @test MCMCChains._discrete_edges(float.(rand(0:3, 400))) ==
              [-0.5, 0.5, 1.5, 2.5, 3.5]
        # Values that are not whole numbers, and not evenly spaced.
        @test MCMCChains._discrete_edges(repeat([0.0, 0.25, 1.0], 100)) ==
              [-0.125, 0.125, 0.625, 1.375]

        # Not this kind of parameter.
        @test MCMCChains._discrete_edges(randn(400)) === nothing
        @test MCMCChains._discrete_edges(randn(15)) === nothing
        @test MCMCChains._discrete_edges(fill(2.0, 100)) === nothing
        @test MCMCChains._discrete_edges(Float64[]) === nothing
        # Too many levels to give each one a bar.
        @test MCMCChains._discrete_edges(float.(rand(1:50, 5000))) === nothing

        # Edges a histogram can actually be drawn over: finite, and each one past the last.
        drawable(edges) = all(isfinite, edges) && issorted(edges, lt = <=)

        # `-0.0` and `0.0` are the same point to a reader, so they share a bar rather than
        # leaving two levels with nothing between them. Every float type, including the
        # ones whose numbers are objects and so are not the same object as their own zero.
        for zero_pair in (
            repeat([-0.0, 0.0], 100),
            repeat([-0.0f0, 0.0f0], 100),
            repeat([BigFloat(-0.0), BigFloat(0.0)], 100),
        )
            @test MCMCChains._discrete_edges(zero_pair) === nothing
            @test drawable(MCMCChains._histogram_edges([zero_pair], :auto))
        end
        @test MCMCChains._discrete_edges(repeat([-0.0, 0.0, 1.0], 100)) == [-0.5, 0.5, 1.5]
        # And no empty bar left where the two were held apart.
        for signed in (
            repeat([-0.0, 0.0, 1.0], 100),
            repeat(BigFloat[BigFloat(-0.0), BigFloat(0.0), BigFloat(1.0)], 100),
        )
            edges = MCMCChains._histogram_edges([signed], :auto)
            @test StatsPlots.Plots._make_hist((signed,), edges).weights == [200, 100]
        end

        # The draws a bar for each value is offered on: at most twenty of them, over at
        # least five times as many draws as there are values. Anything wider is binned by a
        # rule, as the docs say.
        @test MCMCChains._discrete_edges(repeat([0.0, 1.0, 2.0], 5)) !== nothing
        @test MCMCChains._discrete_edges(repeat([0.0, 1.0, 2.0], 4)) === nothing
        @test MCMCChains._discrete_edges(repeat(collect(1.0:20.0), 5)) !== nothing
        @test MCMCChains._discrete_edges(repeat(collect(1.0:21.0), 5)) === nothing
        # Counted across the chain, not on each value, so a parameter that sits on one
        # value and visits the others rarely still gets a bar for each.
        lopsided = vcat(fill(0.0, 13), [1.0, 2.0])
        @test MCMCChains._discrete_edges(lopsided) !== nothing
        @test StatsPlots.Plots._make_hist(
            (lopsided,),
            MCMCChains._histogram_edges([lopsided], :auto),
        ).weights == [13, 1, 1]
        @test MCMCChains._discrete_edges(vcat([0.0, 1.0], fill(2.0, 12))) === nothing

        # Levels far enough out that the midpoint of the pair overflows if taken as a sum.
        @test drawable(MCMCChains._discrete_edges(repeat([1e308, 1.1e308], 100)))

        chn = Chains(float.(rand(1:2, 200, 1, 2)), [:d])
        edges = [s.plotattributes[:bins] for s in histogram(chn).series_list]
        @test all(e -> e == [0.5, 1.5, 2.5], edges)
        # A count the caller asked for still wins.
        @test all(
            s -> length(s.plotattributes[:bins]) == 26,
            histogram(chn, bins = 25).series_list,
        )
        # So does a rule. Only `:auto` leaves the choice to us.
        pooled = float.(rand(1:2, 400))
        for rule in filter(!=(:auto), collect(MCMCChains.BIN_RULES))
            @test MCMCChains._histogram_edges([pooled], rule) isa AbstractRange
        end
    end

    @testset "a rule we do not compute is left to Plots" begin
        # `:wand` is Plots', so it goes back untouched rather than raising.
        @test MCMCChains._histogram_edges([randn(50)], :wand) === :wand

        chn = Chains(randn(100, 1, 2), [:x])
        for rule in (MCMCChains.BIN_RULES..., :wand)
            @test histogram(chn; bins = rule) isa Plots.Plot
        end
    end

    @testset "every chain is binned the same way" begin
        # Two chains that do not overlap, so binning them apart would be obvious.
        series = [randn(500), randn(500) .+ 8]
        edges = MCMCChains._histogram_edges(series, 20)

        # The edges span both chains rather than either one of them.
        @test first(edges) ≈ minimum(Iterators.flatten(series))
        # Past the largest draw, which an open last bin would otherwise leave out.
        @test last(edges) > maximum(Iterators.flatten(series))
        # The count asked for is the count returned.
        @test length(edges) == 21
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
        @test length(auto) == MCMCChains._bin_count(pooled, :auto) + 1
        @test length(MCMCChains._histogram_edges(series, :sturges)) <= length(auto)
        # A count asked for is the count of bins, not one more.
        @test length(MCMCChains._histogram_edges(series, 60)) == 61
    end

    @testset "draws that are not finite are left out of the range" begin
        series = [[1.0, 2.0, Inf], [0.0, 1.5]]
        edges = MCMCChains._histogram_edges(series, 10)
        @test all(isfinite, edges)
        @test first(edges) ≈ 0.0
        @test last(edges) > 2.0
        # Nothing finite to bin on at all, so still a usable set of edges.
        @test all(isfinite, MCMCChains._histogram_edges([[NaN, Inf]], 5))
    end

    @testset "recipe" begin
        chn = Chains(randn(400, 2, 3), [:a, :b])
        p = histogram(chn)
        @test p isa Plots.Plot

        # Every chain of a parameter is drawn over one set of edges.
        edges = [s.plotattributes[:bins] for s in p.series_list]
        @test length(unique(edges)) == length(names(chn))
        @test all(e -> e isa AbstractRange, edges)

        # Parameters may be on scales that share no useful edges, so only chains pool.
        mixed = Chains(hcat(randn(300, 1, 2), 1000 .* randn(300, 1, 2)), [:small, :big])
        by_parameter = histogram(mixed, colordim = :parameter)
        @test all(s -> s.plotattributes[:bins] === :auto, by_parameter.series_list)

        @test histogram(chn, bins = 25) isa Plots.Plot
        @test histogram(Chains(randn(400, 1, 1), [:a])) isa Plots.Plot
        # A discrete parameter reaches the same path from mixeddensity.
        @test mixeddensity(Chains(float.(rand(1:2, 400, 1, 3)), [:d])) isa Plots.Plot
    end
end
