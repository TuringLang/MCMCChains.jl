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
