using MCMCChains, Test
using Statistics: std

@testset "Summarize to DataFrame tests" begin
    val = rand(1000, 8, 4)
    chns = Chains(
        val,
        ["a", "b", "c", "d", "e", "f", "g", "h"],
        Dict(:internals => ["c", "d", "e", "f", "g", "h"])
    )

    parm_df = summarize(chns, sections=[:parameters])

    @test 0.48 < parm_df[:a, :mean][1] < 0.52
    @test names(parm_df) == [:parameters, :mean, :std, :naive_se, :mcse, :ess, :rhat]
    @test parm_df[:a, :] == parm_df[:a]

    all_sections_df = summarize(chns, sections=[:parameters, :internals])
    @test all_sections_df[:,:parameters] == [:a, :b, :c, :d, :e, :f, :g, :h]
    @test size(all_sections_df) == (8, 7)

    two_parms_two_funs_df = summarize(chns[[:a, :b]], mean, std)
    @test two_parms_two_funs_df[:, :parameters] == [:a, :b]
    @test size(two_parms_two_funs_df) == (2, 3)

    three_parms_df = summarize(chns[[:a, :b, :c]], mean, std, sections=[:parameters, :internals])
    @test three_parms_df[:, :parameters] == [:a, :b, :c]
    @test size(three_parms_df) == (3, 3)

    three_parms_df_2 = summarize(chns[[:a, :b, :g]], mean, std,
    sections=[:parameters, :internals], func_names=[:mean, :sd])
    @test three_parms_df_2[:, :parameters] == [:a, :b, :g]
    @test size(three_parms_df_2) == (3, 3)
end
