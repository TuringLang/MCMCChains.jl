using Test
using MCMCChains
using StatsPlots

@testset "Count plots" begin
    @testset "shared categories" begin
        values = [2 1; 0 1; 0 3]
        levels, counts = MCMCChains._value_counts(values)
        @test levels == [0, 1, 2, 3]
        @test counts == [2 0; 0 2; 1 0; 0 1]
        @test MCMCChains._value_counts(reverse(values; dims = 1)) == (levels, counts)
        @test MCMCChains._value_counts(reshape(1:21, :, 1)) ==
              (collect(1:21), ones(Int, 21, 1))
    end

    @testset "values retain their identity" begin
        cases = (
            [0.0, 0.25, 1.0],
            [1.0, nextfloat(1.0)],
            [0.0, nextfloat(0.0)],
            Float16[-60000, 60000],
            [-1e308, 1e308],
            [0.0, floatmax(Float64)],
            [1e16],
            [floatmax(Float64)],
            [typemax(Int) - 1, typemax(Int)],
        )
        for levels in cases
            values = reshape(repeat(levels, 3), :, 1)
            @test MCMCChains._value_counts(values) == (levels, fill(3, length(levels), 1))
            p = countplot(Chains(reshape(values, :, 1, 1)))
            bars = only(filter(s -> s[:seriestype] === :shape, p.series_list))
            positions, labels = only(xticks(p))
            @test labels == string.(levels)
            @test all(==(1), diff(positions))
            for i in eachindex(levels)
                x = bars[:x][(6i-5):(6i-1)]
                @test all(isapprox.(extrema(x), (positions[i] - 0.4, positions[i] + 0.4)))
                @test maximum(bars[:y][(6i-5):(6i-1)]) == 3
            end
        end
        for T in (Float32, Float64, BigFloat)
            values = reshape(T[-0.0, 0.0, 1.0], :, 1)
            @test MCMCChains._value_counts(values) == (T[0, 1], reshape([2, 1], :, 1))
        end
    end

    @testset "mixed BigFloat precision" begin
        low = BigFloat(1; precision = 256)
        high, nearby = setprecision(BigFloat, 512) do
            BigFloat(1), 1 + BigFloat(2)^-300
        end
        before = precision(BigFloat)
        for values in ([low, high, nearby], [nearby, high, low])
            levels, counts = MCMCChains._value_counts(reshape(values, :, 1))
            @test levels == [low, nearby]
            @test counts == reshape([2, 1], :, 1)
            p = countplot(Chains(reshape(values, :, 1, 1)))
            @test only(xticks(p))[2] == string.([low, nearby])
            @test length(unique(only(xticks(p))[2])) == 2
            @test precision(BigFloat) == before
        end
    end

    @testset "public API" begin
        values = reshape([0.0, 0, 1, 1, 2, 2], 3, 1, 2)
        chn = Chains(values, [:a])
        for p in (countplot(chn, :a), plot(chn; seriestype = :countplot))
            bars = filter(s -> s[:seriestype] === :shape, p.series_list)
            @test length(bars) == 2
            @test [s[:label] for s in bars] == ["Chain 1", "Chain 2"]
            @test bars[1][:fillcolor] != bars[2][:fillcolor]
            @test bars[1][:y][1:6:end] == [2, 1, 0]
            @test bars[2][:y][1:6:end] == [0, 1, 2]
            @test bars[1][:x][3:6:end] ≈ bars[2][:x][1:6:end]
        end
        p = countplot(chn; append_chains = true)
        bars = only(filter(s -> s[:seriestype] === :shape, p.series_list))
        @test bars[:y][1:6:end] == [2, 2, 2]

        multi = Chains(cat(values, values .+ 10; dims = 2), [:a, :b])
        p = countplot(multi)
        @test length(p.subplots) == 2
        @test xticks(p[2])[2] == ["10.0", "11.0", "12.0"]
        p = countplot(multi; colordim = :parameter)
        @test length(p.subplots) == 2
        bars = filter(s -> s[:seriestype] === :shape, p.series_list)
        @test [s[:label] for s in bars] == ["a", "b", "a", "b"]
        @test bars[1][:y][1:6:end] == [2, 1, 0, 0]
        @test bars[2][:y][1:6:end] == [0, 0, 2, 1]

        p = countplot(chn; xticks = ([1, 2, 3], ["a", "b", "c"]), ylabel = "Draws")
        @test p[1][:xaxis][:ticks] == ([1, 2, 3], ["a", "b", "c"])
        @test p[1][:yaxis][:guide] == "Draws"
        @test :countplot in MCMCChains.supportedplots
    end

    @testset "missing values and proportions" begin
        values = Union{Missing,Float64}[0 1; 1 missing; 1 missing]
        @test MCMCChains._value_counts(values) == ([0.0, 1], [1 0; 2 1])
        chn = Chains(reshape(values, 3, 1, 2))
        p = countplot(chn; relative = true)
        bars = filter(s -> s[:seriestype] === :shape, p.series_list)
        @test bars[1][:y][1:6:end] ≈ [1 / 3, 2 / 3]
        @test bars[2][:y][1:6:end] == [0, 1]
        @test p[1][:yaxis][:guide] == "Proportion"
        p = countplot(chn; relative = true, append_chains = true)
        bars = only(filter(s -> s[:seriestype] === :shape, p.series_list))
        @test bars[:y][1:6:end] == [0.25, 0.75]

        values[:, 2] .= missing
        chn = Chains(reshape(values, 3, 1, 2))
        @test countplot(chn) isa Plots.Plot
        @test_throws ArgumentError countplot(chn; relative = true)
        @test_throws ArgumentError countplot(Chains(fill(missing, 3, 1, 1)))
        for v in (Inf, -Inf, NaN)
            @test_throws ArgumentError countplot(Chains(reshape([0.0, v], 2, 1, 1)))
        end
    end
end
