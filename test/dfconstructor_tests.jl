using MCMCChains, Test
using DataFrames

@testset "DataFrame constructor tests" begin
    val = rand(1000, 8, 4)
    chn = Chains(val,
                ["a", "b", "c", "d", "e", "f", "g", "h"],
                Dict(:internals => ["c", "d", "e", "f", "g", "h"]))

    df = DataFrame(chn, showall=true)
    @test size(df) == (4000, 8)

    df1 = DataFrame(chn, [:parameters])
    @test size(df1) == (4000, 2)

    df2 = DataFrame(chn, [:internals, :parameters])
    @test size(df2) == (4000, 8)

    df3 = DataFrame(chn[:a])
    @test size(df3) == (4000, 1)

    df4 = DataFrame(chn[:b])
    @test size(df4) == (4000, 1)

    df5 = DataFrame(chn, [:parameters], append_chains=false)
    @test size(df5) == (4, )
    @test size(df5[1]) == (1000, 2)

    df6 = DataFrame(chn, [:parameters, :internals], append_chains=false)
    @test size(df6) == (4, )
    @test size(df6[1]) == (1000, 8)

    df7 = DataFrame(chn, [:parameters, :internals],
        remove_missing_union=false)
    @test size(df7) == (4000, 8)

    df8 = DataFrame(chn, [:parameters, :internals],
        remove_missing_union=false, append_chains=false)
    @test size(df8) == (4, )
    @test size(df8[1]) == (1000, 8)
end

@testset "concretize of DataFrame" begin
    val = rand(10, 4, 4)
    chn = Chains(val, ["a", "b", "c", "d"])

    df = DataFrame(chn)
    @test MCMCChains.concretize(df) === df

    val2 = convert(Array{Union{Missing,Real}}, val)
    chn2 = Chains(val2, ["a", "b", "c", "d"])

    df2 = DataFrame(chn2)
    @test all(x -> eltype(x) === Float64, eachcol(df2))
    @test df2 == df

    df3 = DataFrame(chn2; remove_missing_union = false)
    @test all(x -> eltype(x) === Union{Missing,Real}, eachcol(df3))
    
    df4 = MCMCChains.concretize(df2)
    @test all(x -> eltype(x) === Float64, eachcol(df4))
    @test df4 == df
end

@testset "ChainDataFrame converter" begin
    val = rand(1000, 5, 4)
    chain = Chains(val)

    cdfs = describe(chain)

    dfs = DataFrame(cdfs)
    df1 = DataFrame(cdfs[1])
    df2 = DataFrame(cdfs[2])

    @test dfs[1] == df1
    @test dfs[2] == df2
end
