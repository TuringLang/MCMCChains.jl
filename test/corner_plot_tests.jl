using Test, MCMCChains, StatsPlots, Random

@testset "StatsPlots corner plots" begin
    chn = Chains(
        randn(MersenneTwister(42), 30, 3, 3),
        [:a, :b, :lp],
        Dict(:internals => [:lp]),
    )
    function compare_plots(actual, expected)
        @test actual[:size] == expected[:size]
        @test length(actual.subplots) == length(expected.subplots)
        @test length(actual.series_list) == length(expected.series_list)
        for (a, b) in zip(actual.series_list, expected.series_list)
            @test a[:seriestype] == b[:seriestype]
            @test isequal(a[:x], b[:x])
            @test isequal(a[:y], b[:y])
            @test isequal(a[:label], b[:label])
        end
    end
    for (c, parameters, kwargs, expected_names) in (
        (chn, (), (;), [:a, :b]),
        (chn[:, :, 2:3], (), (; compact = false), [:a, :b]),
        (chn, ([:b, :lp],), (; size = (450, 450)), [:b, :lp]),
        (chn, (), (; sections = nothing), [:a, :b, :lp]),
        (chn, (:a,), (;), [:a]),
        (chn, ([2, 3],), (;), [:b, :lp]),
    )
        @testset "$(parameters), $(kwargs)" begin
            matrix = vcat([Array(c.value[:, expected_names, i]) for i = 1:size(c, 3)]...)
            expected = StatsPlots.cornerplot(
                matrix;
                label = permutedims(string.(expected_names)),
                compact = get(kwargs, :compact, length(expected_names) > 1),
                size = get(kwargs, :size, (600, 600)),
            )
            actual = StatsPlots.cornerplot(c, parameters...; kwargs...)
            compare_plots(actual, expected)
        end
    end
    @testset "append" begin
        matrix = vcat([Array(chn.value[:, [:a, :b], i]) for i = 1:size(chn, 3)]...)
        labels = permutedims(["a", "b"])
        expected =
            StatsPlots.cornerplot(matrix; label = labels, compact = true, size = (600, 600))
        StatsPlots.cornerplot!(expected, matrix; label = labels, compact = true)
        actual = StatsPlots.cornerplot(chn)
        @test StatsPlots.cornerplot!(actual, chn) === actual
        compare_plots(actual, expected)
        StatsPlots.cornerplot!(expected, matrix; label = labels, compact = true)
        plot!(actual)
        @test StatsPlots.cornerplot!(chn) === actual
        compare_plots(actual, expected)
    end
    @testset "render" begin
        previous_backend = StatsPlots.Plots.backend()
        try
            StatsPlots.gr()
            plot = StatsPlots.cornerplot(chn)
            @test occursin("<svg", sprint(show, MIME("image/svg+xml"), plot))
            @test StatsPlots.cornerplot!(plot, chn) === plot
            @test occursin("<svg", sprint(show, MIME("image/svg+xml"), plot))
            custom_margin = 5 * StatsPlots.Plots.mm
            custom = StatsPlots.cornerplot(chn; left_margin = custom_margin)
            @test custom.subplots[1][:left_margin] == custom_margin
            @test occursin("<svg", sprint(show, MIME("image/svg+xml"), custom))
        finally
            StatsPlots.Plots.backend(previous_backend)
        end
    end
end
