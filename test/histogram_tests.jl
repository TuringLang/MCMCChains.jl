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

        # `-0.0` and `0.0` are the same point to a reader, so they share a bar rather than
        # leaving two levels with nothing between them. Every float type, not only the one
        # where `-0.0` is its own object.
        for zero_pair in (
            repeat([-0.0, 0.0], 100),
            repeat([-0.0f0, 0.0f0], 100),
            repeat([BigFloat(-0.0), BigFloat(0.0)], 100),
        )
            @test MCMCChains._discrete_edges(zero_pair) === nothing
        end
        @test MCMCChains._discrete_edges(repeat([-0.0, 0.0, 1.0], 100)) == [-0.5, 0.5, 1.5]

        # Neighbours a single step apart have no value between them to put an edge on, so
        # they share a bar rather than becoming two bars one step wide, which nobody can
        # see. Everything else keeps a bar to itself.
        @test MCMCChains._separable([0.0, 1.0, nextfloat(1.0)]) == [0.0, 1.0]
        @test MCMCChains._separable([1.0, 2.0, 3.0]) == [1.0, 2.0, 3.0]
        # Too far apart to measure between is not the same as too close to fit between, so
        # a pair whose distance overflows keeps a bar each.
        @test MCMCChains._separable([-1.0e308, 1.0e308]) == [-1.0e308, 1.0e308]
        @test MCMCChains._separable(Float16[-60000, 60000]) == Float16[-60000, 60000]
        # A whole run of consecutive values is one bar, not every other member of it.
        function consecutive(start, n)
            values = [start]
            for _ = 2:n
                push!(values, nextfloat(last(values)))
            end
            return values
        end
        @test consecutive(1.0, 4) == [
            1.0,
            nextfloat(1.0),
            nextfloat(nextfloat(1.0)),
            nextfloat(nextfloat(nextfloat(1.0))),
        ]
        @test MCMCChains._separable(consecutive(1.0, 4)) == [1.0]
        @test MCMCChains._separable(vcat(consecutive(1.0, 4), 5.0)) == [1.0, 5.0]

        # Every bar holds draws, every draw is in a bar, and the bar count is the one
        # written down here rather than one read back off the code under test.
        for (values, bars) in (
            (float.(rand(1:2, 400)), 2),
            (float.(rand(0:3, 400)), 4),
            (repeat([0.0, 0.25, 1.0], 100), 3),
            (repeat([-0.0, 0.0, 1.0], 100), 2),
            (repeat(Float32[1, 2, 3], 100), 3),
            # Neighbours with nothing between them share a bar.
            (vcat(fill(0.0, 100), fill(1.0, 100), fill(nextfloat(1.0), 200)), 2),
            (repeat([1.0, nextfloat(1.0)], 100), 1),
            (repeat([0.0, nextfloat(0.0)], 100), 1),
            # A run of them, and a run with something clear of it.
            (repeat(consecutive(1.0, 4), 100), 1),
            (repeat(consecutive(0.0, 9), 50), 1),
            (repeat(vcat(consecutive(1.0, 4), 5.0), 80), 2),
            # Levels at the ends of the range, with nothing beyond them.
            (repeat([1.0e308, 1.1e308], 100), 2),
            (repeat([-1.0e308, 1.0e308], 100), 2),
            (repeat([0.0, floatmax(Float64)], 100), 2),
            (repeat(Float16[-60000, 60000], 100), 2),
            (repeat(Float16[-60000, 0, 60000], 100), 3),
        )
            edges = MCMCChains._discrete_edges(values)
            @test edges !== nothing
            @test all(isfinite, edges)
            @test issorted(edges, lt = <=)
            @test length(edges) == bars + 1

            weights = StatsPlots.Plots._make_hist((values,), edges).weights
            @test count(!iszero, weights) == bars
            @test sum(weights) == length(values)
        end

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

    @testset "every draw lands in a bin" begin
        # Edges are only worth returning when a histogram can be drawn over them, and bins
        # are closed on the left, so the last edge has to sit past the largest draw. Checked
        # over the awkward corners of each element type rather than a few chosen pairs,
        # because this is where overflow and rounding hide.
        function counted(values, bins)
            return sum(StatsPlots.Plots._make_hist((values,), bins).weights)
        end

        # Draws at the very end of the range need an edge past them that `Float64` has no
        # room for, so the plot has nowhere to put it and comes back empty. The draws are
        # still counted, and it is said out loud, for those draws and no others.
        for undrawable in (
            repeat([0.0, floatmax(Float64)], 100),
            fill(floatmax(Float64), 200),
            fill(-floatmax(Float64), 200),
            # Edges a canvas can hold one at a time, but not the width across them.
            repeat([-1.0e308, 1.0e308], 100),
        )
            edges = (@test_logs (:warn,) match_mode = :any MCMCChains._histogram_edges(
                [undrawable],
                :auto,
            ))
            @test counted(undrawable, edges) == length(undrawable)
        end
        for drawable in (
            randn(200),
            repeat([-1.0e305, 1.0e305], 100),
            repeat(Float16[-60000, 60000], 100),
            repeat([0.0, 1.0, nextfloat(1.0)], 100),
            repeat([0.0, nextfloat(0.0)], 100),
            fill(1.0e16, 200),
        )
            @test_logs MCMCChains._histogram_edges([drawable], :auto)
        end

        for T in (Float16, Float32, Float64, BigFloat)
            far = T === BigFloat ? BigFloat(floatmax(Float64)) : floatmax(T)
            corners = T[zero(T), one(T), -one(T), -far, far/2, -far/2, T(0.25)]
            for a in corners, b in corners
                a >= b && continue
                values = repeat(T[a, b], 100)
                for rule in (:auto, :fd, :sturges, :scott, :rice, :sqrt, 25)
                    bins = MCMCChains._histogram_edges([values], rule)
                    # Draws spread wider than their own type reaches leave no edges to
                    # give, and the count goes back to Plots to work out.
                    bins isa AbstractArray || continue
                    @test all(isfinite, bins)
                    @test issorted(bins, lt = <=)
                    @test first(bins) <= a
                    @test b < last(bins)
                    # Plots cannot place a draw in a bin when the span between the two ends
                    # is itself past what the type holds, whatever edges it is given.
                    isfinite(b - a) || continue
                    @test counted(values, bins) == length(values)
                end
            end
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
        # One bar holding the draw, rather than that bar beside empty ones. Half a unit is
        # no step at all once the value is large enough, so the bin widens until it holds.
        for (value, rule) in (
            (2.0, 5),
            (2.0, :auto),
            (1.0e16, :auto),
            (1.0e16, 25),
            (floatmax(Float64), :auto),
            (-floatmax(Float64), :auto),
            (1.0f16, :auto),
            (Float16(2), :auto),
            (nextfloat(0.0), :auto),
        )
            draws = fill(value, 200)
            edges = MCMCChains._histogram_edges([draws], rule)
            @test first(edges) <= value < last(edges)
            @test all(isfinite, edges)
            weights = StatsPlots.Plots._make_hist((draws,), edges).weights
            @test sum(weights) == 200
            @test count(iszero, weights) == 0
        end
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
