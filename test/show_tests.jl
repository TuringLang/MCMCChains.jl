using Test
using MCMCChains

@testset "Show tests" begin
    rng = MersenneTwister(1234)
    val = rand(rng, 100, 4, 4)
    parm_names = ["a", "b", "c", "d"]
    chns = Chains(val, parm_names, Dict(:internals => ["b", "d"]))[1:2:99, :, :]
    str = sprint(show, "text/plain", chns)
    expected_str = """
    Chains MCMC chain (50×4×4 Array{Float64, 3}):

    Iterations        = 1:2:99
    Number of chains  = 4
    Samples per chain = 50
    parameters        = a, c
    internals         = b, d


    Use `describe(chains)` for summary statistics and quantiles.
    """
    @test str == expected_str

    describe_str = sprint(describe, chns)
    @test occursin("Summary Statistics", describe_str)
    @test occursin("Quantiles", describe_str)
    @test occursin("parameters", describe_str)
    @test occursin("internals", describe_str)
end
