using Test
using MCMCChains
using StatsPlots
using Random

@testset "Plot style" begin
    Random.seed!(42)
    chn = Chains(randn(200, 2, 3), [:a, :b])

    try
        @test plot_style().grid === :none
        @test plot_style().legend_panel === :first

        @testset "settings reach the plot" begin
            plot_style!(grid = :dots, framestyle = :box)
            p = meanplot(chn)
            @test p[1][:framestyle] === :box
            # Dotted gridlines are gridlines.
            @test p[1][:xaxis][:grid]

            plot_style!(grid = :none)
            @test !meanplot(chn)[1][:xaxis][:grid]
        end

        @testset "a transparent background is transparent" begin
            plot_style!(background = :transparent)
            p = meanplot(chn)
            @test p[1][:background_color_inside] == RGBA(0, 0, 0, 0)
        end

        @testset "the palette is used for chains" begin
            reset_plot_style!()
            plot_style!(palette = ["#000000", "#FFFFFF"])
            @test MCMCChains._chain_colours(3) == ["#000000" "#FFFFFF" "#000000"]
        end

        @testset "one legend by default" begin
            reset_plot_style!()
            # Three parameters, two series types, so six panels and one legend.
            p = plot(chn)
            @test count(sp -> sp[:legend_position] !== :none, p.subplots) == 1

            plot_style!(legend_panel = :none)
            @test count(sp -> sp[:legend_position] !== :none, plot(chn).subplots) == 0

            plot_style!(legend_panel = :all)
            @test count(sp -> sp[:legend_position] !== :none, plot(chn).subplots) ==
                  length(plot(chn).subplots)
        end

        @testset "what the caller asks for wins" begin
            reset_plot_style!()
            # Plots splits an axis attribute into x, y and z before the recipe runs, so
            # these are the ones a careless default would silently overrule.
            p = plot(chn; grid = true, legend = :topright, tickfontsize = 14)
            @test p[1][:xaxis][:grid]
            @test p[1][:xaxis][:tickfontsize] == 14
            @test all(sp -> sp[:legend_position] === :topright, p.subplots)

            @test plot(chn; framestyle = :box)[1][:framestyle] === :box
            @test !plot(chn)[1][:xaxis][:grid]
        end

        @testset "bad settings are rejected" begin
            @test_throws ArgumentError plot_style!(grid = :squiggles)
            @test_throws ArgumentError plot_style!(legend_panel = :sometimes)
            @test_throws ArgumentError plot_style!(palette = [])
            @test_throws ArgumentError plot_style!(nonsense = 1)
        end
    finally
        reset_plot_style!()
    end

    @test plot_style() == MCMCChains._DEFAULT_STYLE
end
