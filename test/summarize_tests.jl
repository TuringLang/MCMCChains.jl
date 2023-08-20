using MCMCChains, Test
using Statistics: std

@testset "Summarize tests" begin
    val = rand(1000, 8, 4)
    parm_names = ["a", "b", "c", "d", "e", "f", "g", "h"]
    chns = Chains(
        val,
        parm_names,
        Dict(:internals => ["c", "d", "e", "f", "g", "h"])
    )

    parm_df = summarize(chns, sections=[:parameters])
    parm_array_df = summarize(PermutedDimsArray(val[:, 1:2, :], (1, 3, 2)); var_names=[:a, :b])

    # check that display of SummaryStats does not error
    println("compact display:")
    show(stdout, parm_df)
    println("\nverbose display:")
    show(stdout, "text/plain", parm_df)

    @test 0.48 < parm_df[:mean][1] < 0.52
    @test parm_df == parm_array_df

    # Indexing tests
    @test parm_df[:parameter] == [:a, :b]

    all_sections_df = summarize(chns, sections=[:parameters, :internals])
    all_sections_array_df = summarize(PermutedDimsArray(val, (1, 3, 2)); var_names=Symbol.(parm_names))
    @test all_sections_df isa SummaryStats
    @test all_sections_df[:parameter] == Symbol.(parm_names)
    @test all_sections_array_df == all_sections_df
    @test all_sections_df.name == "SummaryStats"

    all_sections_dfs = summarize(chns, sections=[:parameters, :internals], name = "Summary", append_chains = false)
    @test all_sections_dfs isa Vector{<:SummaryStats}
    for (i, all_sections_df) in enumerate(all_sections_dfs)
        @test all_sections_df[:parameter] == Symbol.(parm_names)
        @test length(keys(all_sections_df)) == length(keys(all_sections_array_df))
        @test all_sections_df.name == "Summary (Chain $i)"
    end

    two_parms_two_funs_df = summarize(chns[[:a, :b]], mean, std)
    @test two_parms_two_funs_df[:parameter] == [:a, :b]
    @test keys(two_parms_two_funs_df) == (:parameter, :mean, :std)

    three_parms_df = summarize(chns[[:a, :b, :c]], mean, std, sections=[:parameters, :internals])
    @test three_parms_df[:parameter] == [:a, :b, :c]
    @test keys(three_parms_df) == (:parameter, :mean, :std)

    three_parms_df_2 = summarize(chns[[:a, :b, :g]], :mymean => mean, :mystd => std,
    sections=[:parameters, :internals])
    @test three_parms_df_2[:parameter] == [:a, :b, :g]
    @test keys(three_parms_df_2) == (:parameter, :mymean, :mystd)
end
