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

        for zero_pair in (
            repeat([-0.0, 0.0], 100),
            repeat([-0.0f0, 0.0f0], 100),
            repeat([BigFloat(-0.0), BigFloat(0.0)], 100),
        )
            @test MCMCChains._discrete_edges(zero_pair) === nothing
            fallback = MCMCChains._histogram_edges([zero_pair], :auto)
            @test all(isfinite, fallback)
            @test issorted(fallback, lt = <=)
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
        lopsided = vcat(fill(0.0, 13), [1.0, 2.0])
        @test MCMCChains._discrete_edges(lopsided) !== nothing
        @test StatsPlots.Plots._make_hist(
            (lopsided,),
            MCMCChains._histogram_edges([lopsided], :auto),
        ).weights == [13, 1, 1]
        @test MCMCChains._discrete_edges(vcat([0.0, 1.0], fill(2.0, 12))) === nothing

        function consecutive(start, n)
            values = [start]
            for _ = 2:n
                push!(values, nextfloat(last(values)))
            end
            return values
        end
        @test consecutive(1.0, 3) == [1.0, nextfloat(1.0), nextfloat(nextfloat(1.0))]

        @test MCMCChains._runs([0.0, 1.0, nextfloat(1.0)]) == [1:1, 2:3]
        @test MCMCChains._runs([1.0, 2.0, 3.0]) == [1:1, 2:2, 3:3]
        @test MCMCChains._runs([-1.0e308, 1.0e308]) == [1:1, 2:2]
        @test MCMCChains._runs(Float16[-60000, 60000]) == [1:1, 2:2]
        @test MCMCChains._runs(consecutive(1.0, 4)) == [1:4]
        @test MCMCChains._runs(vcat(consecutive(1.0, 4), 5.0)) == [1:4, 5:5]

        trio = consecutive(1.0f0, 5)
        for (values, expected) in (
            (repeat([1.0, 2.0], 200), [200, 200]),
            (vcat(fill(0.0, 50), fill(1.0, 150), fill(2.0, 100)), [50, 150, 100]),
            (repeat([0.0, 0.25, 1.0], 100), [100, 100, 100]),
            (repeat([-0.0, 0.0, 1.0], 100), [200, 100]),
            (repeat(Float32[1, 2, 3], 100), [100, 100, 100]),
            (vcat(fill(0.0, 100), fill(1.0, 100), fill(nextfloat(1.0), 200)), [100, 300]),
            (repeat([1.0, nextfloat(1.0)], 100), [200]),
            (repeat([0.0, nextfloat(0.0)], 100), [200]),
            (repeat(consecutive(1.0, 4), 100), [400]),
            (repeat(consecutive(0.0, 9), 50), [450]),
            (
                vcat(
                    fill(trio[1], 100),
                    fill(trio[2], 100),
                    fill(trio[3], 100),
                    fill(trio[5], 100),
                ),
                [300, 100],
            ),
            (vcat(fill(1.0, 90), fill(nextfloat(1.0), 10), fill(5.0, 50)), [100, 50]),
            (repeat([1.0e308, 1.1e308], 100), [100, 100]),
            (repeat([-1.0e308, 1.0e308], 100), [100, 100]),
            (repeat([0.0, floatmax(Float64)], 100), [100, 100]),
            (repeat(Float16[-60000, 60000], 100), [100, 100]),
            (repeat(Float16[-60000, 0, 60000], 100), [100, 100, 100]),
        )
            edges = MCMCChains._discrete_edges(values)
            @test edges !== nothing
            @test all(isfinite, edges)
            @test issorted(edges, lt = <=)
            @test length(edges) == length(expected) + 1
            @test StatsPlots.Plots._make_hist((values,), edges).weights == expected
        end

        for (crowded, bars) in (
            (repeat([prevfloat(1.0), 1.0], 100), 1),
            (repeat([1.0, nextfloat(1.0)], 100), 1),
            (repeat([-1.0e308, 1.0e308], 100), 2),
        )
            edges = MCMCChains._discrete_edges(crowded)
            @test edges !== nothing
            @test length(edges) - 1 == bars
            weights = StatsPlots.Plots._make_hist((crowded,), edges).weights
            @test count(!iszero, weights) == bars
            @test sum(weights) == length(crowded)
        end

        far_apart = repeat([-1.0e308, 1.0e308], 100)
        for scale in (5.0e307, 1.0e308)
            wide = repeat([-scale, scale], 100)
            @test MCMCChains._bin_count(wide, :fd) == 3
            @test MCMCChains._bin_count(wide, :auto) == 3
            @test MCMCChains._bin_count(wide, :scott) == 4
        end
        for rule in (:auto, :fd, :sturges, :scott, :rice, :sqrt, 25)
            edges = MCMCChains._histogram_edges([far_apart], rule)
            @test edges isa AbstractArray
            @test StatsPlots.Plots._make_hist((far_apart,), edges).weights |> sum ==
                  length(far_apart)
        end

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

    @testset "every draw lands in a bin" begin
        function counted(values, bins)
            return sum(StatsPlots.Plots._make_hist((values,), bins).weights)
        end

        for undrawable in (
            repeat([0.0, floatmax(Float64)], 100),
            fill(floatmax(Float64), 200),
            fill(-floatmax(Float64), 200),
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
                    bins isa AbstractArray || continue
                    @test all(isfinite, bins)
                    @test issorted(bins, lt = <=)
                    @test first(bins) <= a
                    @test b < last(bins)
                    @test counted(values, bins) == length(values)
                end
            end
        end
    end

    @testset "a rule we do not compute is left to Plots" begin
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
        for (value, rule, expected_bins) in (
            (2.0, :auto, 1),
            (2.0, :fd, 1),
            (2.0, 5, 5),
            (1.0e16, :auto, 1),
            (1.0e16, 25, 25),
            (floatmax(Float64), :auto, 1),
            (-floatmax(Float64), :auto, 1),
            (-floatmax(Float64), 4, 4),
            (1.0f16, :auto, 1),
            (Float16(2), :auto, 1),
            (Float16(2), 3, 3),
            (nextfloat(0.0), :auto, 1),
            (BigFloat(1.0e16), 25, 25),
            (floatmax(Float64), 25, 25),
            (-floatmax(Float64), 25, 25),
        )
            draws = fill(value, 200)
            edges = MCMCChains._histogram_edges([draws], rule)
            @test first(edges) <= value < last(edges)
            @test all(isfinite, edges)
            @test length(edges) - 1 == expected_bins
            @test issorted(collect(edges), lt = <=)
            @test length(unique(collect(edges))) == expected_bins + 1
            @test sum(StatsPlots.Plots._make_hist((draws,), edges).weights) == 200
        end
    end

    @testset "BigFloat draws keep their precision" begin
        lower = setprecision(BigFloat, 512) do
            BigFloat(1) - BigFloat(2)^(-300)
        end
        setprecision(BigFloat, 256) do
            far = repeat([lower, BigFloat(2)], 100)
            given = range(BigFloat(0), BigFloat(2); length = 3)
            @test MCMCChains._histogram_edges([far], given) === given
            edges = MCMCChains._histogram_edges([far], 25)
            @test first(edges) <= lower < last(edges)
            @test length(unique(collect(edges))) == 26
            @test sum(StatsPlots.Plots._make_hist((far,), edges).weights) == 200
            plot = histogram(Chains(reshape(far, :, 1, 1), [:x]); bins = 25)
            plot_edges = plot.series_list[1].plotattributes[:bins]
            @test sum(StatsPlots.Plots._make_hist((far,), plot_edges).weights) == 200

            close = repeat([lower, BigFloat(1)], 100)
            @test MCMCChains._bin_count(close, :fd) == 3
            @test MCMCChains._bin_count(close, :scott) == 4
            for rule in (:auto, :fd, :scott)
                bins = MCMCChains._histogram_edges([close], rule)
                @test first(bins) <= lower < last(bins)
                @test sum(StatsPlots.Plots._make_hist((close,), bins).weights) == 200
            end
            @test StatsPlots.Plots._make_hist(
                (close,),
                MCMCChains._discrete_edges(close),
            ).weights == [200]
            @test precision(BigFloat) == 256
        end
    end

    @testset "mixed BigFloat precisions share narrow gaps" begin
        low = setprecision(BigFloat, 256) do
            BigFloat(1)
        end
        high = setprecision(BigFloat, 512) do
            BigFloat(1) + BigFloat(2)^(-300)
        end
        lower = setprecision(BigFloat, 512) do
            BigFloat(1) - BigFloat(2)^(-300)
        end
        for pair in ((low, high), (lower, low))
            draws = repeat(collect(pair), 100)
            @test MCMCChains._runs(collect(pair)) == [1:2]
            bins = MCMCChains._histogram_edges([draws], :auto)
            @test StatsPlots.Plots._make_hist((draws,), bins).weights == [200]
            @test length(unique(Float64.(bins))) == 2
            plot = histogram(Chains(reshape(draws, :, 1, 1), [:x]))
            plot_bins = plot.series_list[1].plotattributes[:bins]
            @test StatsPlots.Plots._make_hist((draws,), plot_bins).weights == [200]
        end
        @test MCMCChains._runs([lower, BigFloat(2)]) == [1:1, 2:2]
        beyond = setprecision(BigFloat, 512) do
            BigFloat(nextfloat(low)) + BigFloat(2)^(-300)
        end
        @test MCMCChains._runs([low, beyond]) == [1:1, 2:2]
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
