using Test
using MCMCChains
using StatsPlots
using Random

@testset "Histograms" begin
    Random.seed!(42)

    @testset "how many bins" begin
        for values in (randn(100), randn(5000), 3 .* randn(731) .+ 2, rand(50))
            for rule in MCMCChains.BIN_RULES
                @test MCMCChains._bin_count(values, rule) ==
                      Plots._auto_binning_nbins((values,), 1; mode = rule)
            end
        end

        values = randn(2000)
        @test MCMCChains._bin_count(values, :auto) == MCMCChains._bin_count(values, :fd)
        @test MCMCChains._bin_count(randn(20000), :fd) >
              MCMCChains._bin_count(randn(200), :fd)

        @test MCMCChains._bin_count(fill(2.0, 100), :auto) == 1
        @test MCMCChains._bin_count([1.0], :auto) == 1
        spiked = vcat(fill(0.0, 999), 10.0)
        for rule in MCMCChains.BIN_RULES
            @test MCMCChains._bin_count(spiked, rule) >= 1
        end

        @test_throws ArgumentError MCMCChains._bin_count(values, :nonsense)
    end

    @testset "a parameter on a few values is binned on them" begin
        @test MCMCChains._discrete_edges(float.(rand(1:2, 400))) == [0.5, 1.5, 2.5]
        @test MCMCChains._discrete_edges(float.(rand(0:3, 400))) ==
              [-0.5, 0.5, 1.5, 2.5, 3.5]
        @test MCMCChains._discrete_edges(repeat([0.0, 0.25, 1.0], 100)) ==
              [-0.125, 0.125, 0.625, 1.375]

        @test MCMCChains._discrete_edges(randn(400)) === nothing
        @test MCMCChains._discrete_edges(randn(15)) === nothing
        @test MCMCChains._discrete_edges(fill(2.0, 100)) === nothing
        @test MCMCChains._discrete_edges(Float64[]) === nothing
        @test MCMCChains._discrete_edges(float.(rand(1:50, 5000))) === nothing

        drawable(edges) = all(isfinite, edges) && issorted(edges, lt = <=)

        # Include BigFloat, whose signed zeros are distinct objects.
        for zero_pair in (
            repeat([-0.0, 0.0], 100),
            repeat([-0.0f0, 0.0f0], 100),
            repeat([BigFloat(-0.0), BigFloat(0.0)], 100),
        )
            @test MCMCChains._discrete_edges(zero_pair) === nothing
            @test drawable(MCMCChains._histogram_edges([zero_pair], :auto))
        end
        @test MCMCChains._discrete_edges(repeat([-0.0, 0.0, 1.0], 100)) == [-0.5, 0.5, 1.5]
        for signed in (
            repeat([-0.0, 0.0, 1.0], 100),
            repeat(BigFloat[BigFloat(-0.0), BigFloat(0.0), BigFloat(1.0)], 100),
        )
            edges = MCMCChains._histogram_edges([signed], :auto)
            @test StatsPlots.Plots._make_hist((signed,), edges).weights == [200, 100]
        end

        @test MCMCChains._discrete_edges(repeat([0.0, 1.0, 2.0], 5)) !== nothing
        @test MCMCChains._discrete_edges(repeat([0.0, 1.0, 2.0], 4)) === nothing
        @test MCMCChains._discrete_edges(repeat(collect(1.0:20.0), 5)) !== nothing
        @test MCMCChains._discrete_edges(repeat(collect(1.0:21.0), 5)) === nothing
        # The draw threshold applies across levels, including rarely visited ones.
        lopsided = vcat(fill(0.0, 13), [1.0, 2.0])
        @test MCMCChains._discrete_edges(lopsided) !== nothing
        @test StatsPlots.Plots._make_hist(
            (lopsided,),
            MCMCChains._histogram_edges([lopsided], :auto),
        ).weights == [13, 1, 1]
        @test MCMCChains._discrete_edges(vcat([0.0, 1.0], fill(2.0, 12))) === nothing

        @test drawable(MCMCChains._discrete_edges(repeat([1e308, 1.1e308], 100)))

        for crowded in (
            repeat([prevfloat(1.0), 1.0], 100),
            repeat([1.0, nextfloat(1.0)], 100),
            repeat([-1.0e308, 1.0e308], 100),
        )
            @test MCMCChains._discrete_edges(crowded) === nothing
        end
        kept = repeat([prevfloat(1.0), 1.0], 100)
        @test sum(
            StatsPlots.Plots._make_hist(
                (kept,),
                MCMCChains._histogram_edges([kept], :auto),
            ).weights,
        ) == length(kept)

        chn = Chains(float.(rand(1:2, 200, 1, 2)), [:d])
        edges = [s.plotattributes[:bins] for s in histogram(chn).series_list]
        @test all(e -> e == [0.5, 1.5, 2.5], edges)
        @test all(
            s -> length(s.plotattributes[:bins]) == 26,
            histogram(chn, bins = 25).series_list,
        )
        pooled = float.(rand(1:2, 400))
        for rule in filter(!=(:auto), collect(MCMCChains.BIN_RULES))
            @test MCMCChains._histogram_edges([pooled], rule) isa AbstractRange
        end
    end

    @testset "a rule we do not compute is left to Plots" begin
        # `:wand` is handled by Plots.
        @test MCMCChains._histogram_edges([randn(50)], :wand) === :wand

        chn = Chains(randn(100, 1, 2), [:x])
        for rule in (MCMCChains.BIN_RULES..., :wand)
            @test histogram(chn; bins = rule) isa Plots.Plot
        end
    end

    @testset "every chain is binned the same way" begin
        series = [randn(500), randn(500) .+ 8]
        edges = MCMCChains._histogram_edges(series, 20)

        @test first(edges) ≈ minimum(Iterators.flatten(series))
        @test last(edges) > maximum(Iterators.flatten(series))
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
        @test length(MCMCChains._histogram_edges(series, 60)) == 61
    end

    @testset "draws that are not finite are left out of the range" begin
        series = [[1.0, 2.0, Inf], [0.0, 1.5]]
        edges = MCMCChains._histogram_edges(series, 10)
        @test all(isfinite, edges)
        @test first(edges) ≈ 0.0
        @test last(edges) > 2.0
        @test all(isfinite, MCMCChains._histogram_edges([[NaN, Inf]], 5))
    end

    @testset "recipe" begin
        chn = Chains(randn(400, 2, 3), [:a, :b])
        p = histogram(chn)
        @test p isa Plots.Plot

        edges = [s.plotattributes[:bins] for s in p.series_list]
        @test length(unique(edges)) == length(names(chn))
        @test all(e -> e isa AbstractRange, edges)

        mixed = Chains(hcat(randn(300, 1, 2), 1000 .* randn(300, 1, 2)), [:small, :big])
        by_parameter = histogram(mixed, colordim = :parameter)
        @test all(s -> s.plotattributes[:bins] === :auto, by_parameter.series_list)

        @test histogram(chn, bins = 25) isa Plots.Plot
        @test histogram(Chains(randn(400, 1, 1), [:a])) isa Plots.Plot
        @test mixeddensity(Chains(float.(rand(1:2, 400, 1, 3)), [:d])) isa Plots.Plot
    end
end
