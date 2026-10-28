using Test
using MCMCChains
using Statistics
using LinearAlgebra
using AxisArrays: Axis, AxisArray

@testset "Summary statistics" begin
    values = reshape([1, 3, 5, 2, 4, 8, 11, 13, 15, 3, 6, 9], 3, 2, 2)
    chn = Chains(values, [:a, :b])
    pooled = vcat(values[:, :, 1], values[:, :, 2])

    @testset "var and std" begin
        for (f, column) in ((var, :var), (std, :std))
            @test f(chn) isa ChainDataFrame
            @test f(chn).nt[column] ≈ [f(pooled[:, i]) for i = 1:2]
            for corrected in (true, false)
                @test f(chn; corrected).nt[column] ≈
                      [f(pooled[:, i]; corrected) for i = 1:2]
                per_chain = f(chn; append_chains = false, corrected)
                @test length(per_chain) == 2
                for i = 1:2
                    @test per_chain[i].nt[column] ≈
                          [f(values[:, j, i]; corrected) for j = 1:2]
                end
            end
            reordered = f(Chains(values[:, [2, 1], :], [:b, :a]))
            @test reordered.nt.parameters == [:b, :a]
            @test reordered.nt[column] ≈ [f(pooled[:, i]) for i in (2, 1)]
        end
    end

    @testset "covariance" begin
        result = cov(chn)
        @test result isa AxisArray
        @test size(result) == (2, 2)
        @test result ≈ cov(pooled)
        @test diag(result) ≈ var(chn).nt.var
        @test cholesky(Symmetric(result)).L ≈ cholesky(Symmetric(cov(pooled))).L
        @test result * [1.0, 2.0] ≈ cov(pooled) * [1.0, 2.0]
        @test result[Axis{:row}(:a), Axis{:col}(:b)] ≈ cov(pooled)[1, 2]
        for corrected in (true, false)
            @test cov(chn; corrected) ≈ cov(pooled; corrected)
            per_chain = cov(chn; append_chains = false, corrected)
            @test length(per_chain) == 2
            for i = 1:2
                @test per_chain[i] ≈ cov(values[:, :, i]; corrected)
                @test per_chain[i][Axis{:row}(:b), Axis{:col}(:a)] ≈
                      cov(values[:, :, i]; corrected)[2, 1]
            end
        end
        reordered = cov(Chains(values[:, [2, 1], :], [:b, :a]))
        @test reordered ≈ cov(pooled[:, [2, 1]])
        @test reordered[Axis{:row}(:a), Axis{:col}(:b)] ≈ result[1, 2]
        named = Chains(values, [:parameters, :b])
        @test cov(named)[Axis{:row}(:parameters), Axis{:col}(:b)] ≈ result[1, 2]
        @test_throws MethodError cov(chn; correctedd = false)
    end

    @testset "sections" begin
        sectioned = Chains(values, [:a, :lp], Dict(:internals => [:lp]))
        for f in (var, std, cov)
            @test size(f(sectioned))[1] == 1
            @test size(f(sectioned; sections = nothing))[1] == 2
            @test size(f(sectioned; sections = :internals))[1] == 1
        end
        @test var(sectioned).nt.parameters == [:a]
        @test std(sectioned; sections = :internals).nt.parameters == [:lp]
        @test cov(sectioned; sections = :internals)[1, 1] ≈ var(pooled[:, 2])
    end

    @testset "missing and constant draws" begin
        data = Union{Missing,Float64}[1 2; missing 4; 3 6]
        missing_chain = Chains(data, [:a, :b])
        for f in (var, std)
            @test f(missing_chain).nt[Symbol(f)] ≈ [f([1.0, 3.0]), f([2.0, 4.0, 6.0])]
            @test f(missing_chain; corrected = false).nt[Symbol(f)] ≈
                  [f([1.0, 3.0]; corrected = false), f([2.0, 4.0, 6.0]; corrected = false)]
        end
        @test isequal(cov(missing_chain), cov(data))
        for data in (
            fill(2.0, 3, 1, 1),
            fill(2.0, 1, 1, 1),
            fill(2, 3, 1, 1),
            fill!(Array{Union{Missing,Float64}}(undef, 3, 1, 1), missing),
        )
            c = Chains(data, [:a])
            for corrected in (true, false)
                for f in (var, std)
                    @test isequal(
                        f(c; corrected).nt[Symbol(f)][1],
                        f(collect(skipmissing(vec(data))); corrected),
                    )
                end
                @test isequal(cov(c; corrected), cov(reshape(data, :, 1); corrected))
            end
        end
    end
end
