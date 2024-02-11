using Test
using MCMCChains

@testset "Show tests" begin
    rng = MersenneTwister(1234)
    val = rand(rng, 100, 4, 4)
    parm_names = ["a", "b", "c", "d"]
    chns = Chains(val, parm_names, Dict(:internals => ["b", "d"]))[1:2:99, :, :]
    str = sprint(show, "text/plain", chns)
    stats_str = sprint(show, "text/plain", summarystats(chns))
    quantile_str = sprint(show, "text/plain", quantile(chns))
    @test str == """Chains MCMC chain (50×4×4 Array{Float64, 3}):

    Iterations        = 1:2:99
    Number of chains  = 4
    Samples per chain = 50
    parameters        = a, c
    internals         = b, d

    $stats_str

    $quantile_str"""
end
