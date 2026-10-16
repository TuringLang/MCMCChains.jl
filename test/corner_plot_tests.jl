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
        end
    end
    for (c, parameters, kwargs) in (
        (chn, (), (;)),
        (chn[:, :, 2:3], (), (; compact = false)),
        (chn, ([:b, :lp],), (; size = (450, 450))),
        (chn, (), (; sections = nothing)),
    )
        @testset "$(parameters), $(kwargs)" begin
            expected = corner(c, parameters...; kwargs...)
            actual = StatsPlots.cornerplot(c, parameters...; kwargs...)
            compare_plots(actual, expected)
        end
    end
    @testset "append" begin
        expected = corner(chn)
        corner!(expected, chn)
        actual = StatsPlots.cornerplot(chn)
        @test StatsPlots.cornerplot!(actual, chn) === actual
        compare_plots(actual, expected)
        corner!(expected, chn)
        plot!(actual)
        @test StatsPlots.cornerplot!(chn) === actual
        compare_plots(actual, expected)
    end
end
