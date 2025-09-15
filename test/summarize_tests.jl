using DataFrames
using MCMCChains, Test
using PosteriorStats: SummaryStats
using Statistics: std

@testset "Summarize tests" begin
    val = rand(1000, 8, 4)
    parm_names = ["a", "b", "c", "d", "e", "f", "g", "h"]
    chns = Chains(
        val,
        parm_names,
        Dict(:internals => ["c", "d", "e", "f", "g", "h"]),
    )

    parm_stats = summarize(chns, sections=[:parameters])
    @test parm_stats isa SummaryStats
    @test parm_stats.name == "SummaryStats"
    parm_array_stats = summarize(
        PermutedDimsArray(val[:, 1:2, :], (1, 3, 2));
        var_names =[:a, :b],
    )

    # check that display of SummaryStats does not error
    println("compact display:")
    show(stdout, parm_stats)
    println("\nverbose display:")
    show(stdout, "text/plain", parm_stats)

    parm_df = DataFrame(parm_stats)
    parm_array_df = DataFrame(parm_array_stats)

    @test 0.48 < parm_df[1, :mean] < 0.52
    @test parm_df == parm_array_df

    # Indexing tests
    @test parm_df[!, 1] == [:a, :b]

    all_sections_stats = summarize(chns; sections = [:parameters, :internals])
    all_sections_array_stats = summarize(
        PermutedDimsArray(val, (1, 3, 2)); var_names = Symbol.(parm_names),
    )
    @test all_sections_stats isa SummaryStats
    all_sections_df = DataFrame(all_sections_stats)
    all_sections_array_df = DataFrame(all_sections_array_stats)
    @test all_sections_df[!, 1] == Symbol.(parm_names)
    @test all_sections_array_df == all_sections_df
    @test all_sections_stats.name == "SummaryStats"

    all_sections_stats = summarize(
        chns;
        sections = [:parameters, :internals],
        name = "Summary",
        append_chains = false,
    )
    @test all_sections_stats isa Vector{<:SummaryStats}
    for (i, all_sections_stats_i) in enumerate(all_sections_stats)
        all_sections_df_i = DataFrame(all_sections_stats_i)
        @test all_sections_df_i[!, 1] == Symbol.(parm_names)
        @test size(all_sections_df_i, 2) == size(all_sections_df, 2)
        @test all_sections_stats_i.name == "Summary (Chain $i)"
    end

    two_parms_two_funs_df = DataFrame(summarize(chns[[:a, :b]], mean, std))
    @test two_parms_two_funs_df[!, 1] == [:a, :b]
    @test propertynames(two_parms_two_funs_df) == [:label, :mean, :std]

    three_parms_df = DataFrame(summarize(
        chns[[:a, :b, :c]],
        mean, std;
        sections = [:parameters, :internals],
    ))
    @test three_parms_df[!, 1] == [:a, :b, :c]
    @test propertynames(three_parms_df) == [:label, :mean, :std]

    three_parms_df_2 = DataFrame(summarize(
        chns[[:a, :b, :g]],
        :mymean => mean,
        :mystd => std;
        sections = [:parameters, :internals],
    ))
    @test three_parms_df_2[!, 1] == [:a, :b, :g]
    @test propertynames(three_parms_df_2) == [:label, :mymean, :mystd]
end
