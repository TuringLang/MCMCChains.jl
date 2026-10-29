using MCMCChains, Tables, Test
using Statistics: std

@testset "Summarize to DataFrame tests" begin
    val = rand(1000, 8, 4)
    chns = Chains(
        val,
        ["a", "b", "c", "d", "e", "f", "g", "h"],
        Dict(:internals => ["c", "d", "e", "f", "g", "h"]),
    )

    parm_df = summarize(chns, sections = [:parameters])

    # check that display of ChainDataFrame does not error
    println("compact display:")
    show(stdout, parm_df)
    println("\nverbose display:")
    show(stdout, "text/plain", parm_df)

    @test 0.48 < parm_df[:a, :mean][1] < 0.52
    @test names(parm_df) ==
          [:parameters, :mean, :std, :mcse, :ess_bulk, :ess_tail, :rhat, :ess_per_sec]

    # Indexing tests
    @test isequal(convert(Array, parm_df[:a, :]), convert(Array, parm_df[:a]))
    @test parm_df[:a, :][:, :parameters] == :a
    @test parm_df[[:a, :b], :][:, :parameters] == [:a, :b]

    all_sections_df = summarize(chns, sections = [:parameters, :internals])
    @test all_sections_df isa ChainDataFrame
    @test all_sections_df[:, :parameters] == [:a, :b, :c, :d, :e, :f, :g, :h]
    @test size(all_sections_df) == (8, 8)
    @test all_sections_df.name == ""

    all_sections_dfs = summarize(
        chns,
        sections = [:parameters, :internals],
        name = "Summary",
        append_chains = false,
    )
    @test all_sections_dfs isa Vector{<:ChainDataFrame}
    for (i, all_sections_df) in enumerate(all_sections_dfs)
        @test all_sections_df[:, :parameters] == [:a, :b, :c, :d, :e, :f, :g, :h]
        @test size(all_sections_df) == (8, 8)
        @test all_sections_df.name == "Summary (Chain $i)"
    end

    two_parms_two_funs_df = summarize(chns[[:a, :b]], mean, std)
    @test two_parms_two_funs_df[:, :parameters] == [:a, :b]
    @test size(two_parms_two_funs_df) == (2, 3)

    three_parms_df =
        summarize(chns[[:a, :b, :c]], mean, std, sections = [:parameters, :internals])
    @test three_parms_df[:, :parameters] == [:a, :b, :c]
    @test size(three_parms_df) == (3, 3)

    three_parms_df_2 = summarize(
        chns[[:a, :b, :g]],
        mean,
        std,
        sections = [:parameters, :internals],
        func_names = [:mean, :sd],
    )
    @test three_parms_df_2[:, :parameters] == [:a, :b, :g]
    @test size(three_parms_df_2) == (3, 3)
end

@testset "summary array conversion" begin
    chn = Chains(reshape(1:12, 3, 2, 2), ["x[2]", "x[1]"])
    summary = mean(group(chn, :x))
    expected = reshape([5.0, 8.0], 2, 1)
    @test summary.nt.parameters == [Symbol("x[2]"), Symbol("x[1]")]
    @test Array(summary) == expected
    @test vec(Array(summary)) == Tables.getcolumn(summary, :mean)
    @test Array{Float32}(summary) == Float32.(expected)
    @test convert(Array, summary) == expected
    @test eltype(convert(Array, summary)) === Float64
    single = mean(group(chn, Symbol("x[2]")))
    @test size(Array(single)) == (1, 1)
    @test vec(Array(single)) == Tables.getcolumn(single, :mean)
    per_chain = mean(chn; append_chains = false)
    @test convert(Array, per_chain) == reshape([2.0, 5.0, 8.0, 11.0], 2, 1, 2)
    for nt in (
        (parameters = [:a, :b], mean = Float32[1, 2], count = [3, 4]),
        (parameters = [:a, :b], mean = Union{Missing,Float64}[1, missing]),
        (parameters = Symbol[], mean = Float64[]),
    )
        df = ChainDataFrame(nt)
        expected = hcat(values(Base.tail(nt))...)
        @test isequal(Array(df), expected)
        @test eltype(Array(df)) === eltype(expected)
    end
end
