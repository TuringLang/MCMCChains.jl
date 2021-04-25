using Test, Random
using MCMCChains
using Distributions

@testset "deviance information criterion" begin
    # Ensure function at least runs
    chain = Chains(rand(100, 2, 1), [:a, :b])
    val = dic(chain, :a) 
    @test isa(val, Float64)
    # Should fail if variable does not exist
    @test_throws ArgumentError dic(chain, :c) 
end
