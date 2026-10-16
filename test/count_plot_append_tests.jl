using Test
using MCMCChains
using StatsPlots

@testset "Appending count plots" begin
    chain(values) = Chains(reshape(values, :, 1, 1))
    centers(p) = [s[:x] for s in p.series_list if s[:seriestype] === :scatter]

    @testset "overlapping and new values" begin
        p = countplot(chain([0, 0, 1]))
        before = copy(only(centers(p)))
        countplot!(p, chain([1.0, 2.0, 2.0]))
        positions, labels = deepcopy(xticks(p[1]))
        @test labels == ["0", "1", "2.0"]
        @test positions[1:2] == before
        @test centers(p)[1] == before
        @test centers(p)[2][1] == before[2]
        @test centers(p)[2][2] > before[2]
        @test length(unique(positions)) == 3

        countplot!(p, chain([0, 1]))
        @test centers(p)[3] == before
        @test xticks(p[1]) == (positions, labels)
        countplot!(p, chain([-1, 3]))
        @test xticks(p[1])[2] == ["0", "1", "2.0", "-1", "3"]
        @test xticks(p[1])[1][1:3] == positions
    end

    @testset "custom labels" begin
        p = countplot(chain([0, 1]))
        positions = copy(only(centers(p)))
        custom = (positions, ["zero", "one"])
        plot!(p; xticks = custom)
        countplot!(p, chain([1, 2]))
        @test xticks(p[1]) == custom
        @test centers(p)[1] == positions
        @test centers(p)[2][1] == positions[2]
    end

    @testset "separate subplots" begin
        p = plot(countplot(chain([0, 1])), countplot(chain([10, 11])); layout = (1, 2))
        first_ticks = xticks(p[1])
        before = deepcopy(centers(p))
        countplot!(p, chain([11, 12]); subplot = 2)
        @test xticks(p[1]) == first_ticks
        @test xticks(p[2])[2] == ["10", "11", "12"]
        @test centers(p)[1:2] == before
        @test centers(p)[3][1] == before[2][2]
    end

    @testset "exact values across calls" begin
        for values in (
            [typemax(Int) - 1, typemax(Int)],
            [1.0, nextfloat(1.0)],
            [BigFloat(1), nextfloat(BigFloat(1))],
        )
            p = countplot(chain(values[2:2]))
            before = only(only(centers(p)))
            countplot!(p, chain(values))
            @test xticks(p[1])[2] == string.(reverse(values))
            @test centers(p)[2][2] == before
            @test centers(p)[2][1] > before
        end
    end

    @testset "grouped chains" begin
        a = Chains(reshape([0, 0, 1, 0, 1, 1], 3, 1, 2))
        b = Chains(reshape([1, 1, 2, 1, 2, 2], 3, 1, 2))
        p = countplot(a; relative = true)
        before = deepcopy(centers(p))
        countplot!(p, b; relative = true)
        @test xticks(p[1])[2] == ["0", "1", "2"]
        @test centers(p)[1:2] == before
        @test centers(p)[3][1] == before[1][2]
        @test centers(p)[4][1] == before[2][2]
    end
end
