using Test
using MCMCChains
using StatsPlots

@testset "Interval plot helpers" begin
    chn = Chains(randn(200, 3, 2), [:a, :b, :c])

    @testset "_interval_ylims" begin
        # Rows sit at riser .+ (0:nparams-1) .* spacer, so the limits must clear both ends
        # by half a row. Without this the outermost row is drawn on the frame.
        lo, hi = MCMCChains._interval_ylims(0.2, 0.5, 3, -Inf)
        @test lo < 0.2
        @test hi > 0.2 + 2 * 0.5
        @test 0.2 - lo ≈ 0.25
        @test hi - (0.2 + 2 * 0.5) ≈ 0.25

        # A ridge taller than the last baseline must not be clipped.
        _, hi_tall = MCMCChains._interval_ylims(0.2, 0.5, 3, 10.0)
        @test hi_tall > 10.0
    end

    @testset "_interval_rows" begin
        rows = MCMCChains._interval_rows(chn, [:a, :b, :c])
        @test length(rows) == 3
        @test all(haskey(r, :val) && haskey(r, :h) for r in rows)
        # Baselines increase, which is what the axis limits are derived from.
        @test issorted([r.h for r in rows])
    end

    @testset "argument handling" begin
        # Omitting the parameter list previously threw a BoundsError from inside the recipe.
        @test ridgelineplot(chn) isa Plots.Plot
        @test forestplot(chn) isa Plots.Plot
        @test ridgelineplot(chn, [:a, :b]) isa Plots.Plot

        @test_throws ArgumentError ridgelineplot(chn, [:a], :extra)
        @test_throws ArgumentError forestplot(chn, [:a], :extra)
    end
end
