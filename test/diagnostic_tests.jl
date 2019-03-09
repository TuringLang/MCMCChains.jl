using MCMCChains
using Test

## CHAIN TESTS
# Define the experiment
n_iter = 500
n_name = 3
n_chain = 2

# some sample experiment results
val = randn(n_iter, n_name, n_chain) .+ [1, 2, 3]'
val = hcat(val, rand(1:2, n_iter, 1, n_chain))

# construct a Chains object
chn = Chains(val, start = 1, thin = 2)

@testset "basic chains functions" begin
    @test first(chn) == 1
    @test step(chn) == 2
    @test last(chn) == 999
    @test size(chn) == (999, 4, 2)
    @test keys(chn) == ["Param1", "Param2", "Param3", "Param4"]
    @test isa(chn[:,1,:], MCMCChains.AbstractChains)
    @test isa(chn[200:300,"Param1",:], MCMCChains.AbstractChains)
    @test isa(chn[200:300,["Param1", "Param3"],:], MCMCChains.AbstractChains)
    @test isa(chn[200:300,"Param1",1], MCMCChains.AbstractChains)
    @test length(vec(chn[:,1,:].value)) == n_chain * n_iter
    @test all(collect(skipmissing(chn[:,1,1].value)) .== val[:,1,1])
    @test all(chn[:,1,2].value .== val[:,1,2])
    @test all(MCMCChains.indiscretesupport(chn) .== [false, false, false, true])
    @test setinfo(chn, NamedTuple{(:A, :B)}((1,2))).info == NamedTuple{(:A, :B)}((1,2))
end

@testset "function tests" begin
    # the following tests only check if the function calls work!
    @test MCMCChains.diag_all(rand(100, 2), :weiss, 1, 1, 1) != nothing
    @test MCMCChains.diag_all(rand(100, 2), :hangartner, 1, 1, 1) != nothing
    @test MCMCChains.diag_all(rand(100, 2), :billingsley, 1, 1, 1) != nothing

    @test isa(discretediag(chn[:,4,:]), MCMCChains.ChainSummary)
    @test isa(gelmandiag(chn[:,1,:]), MCMCChains.ChainSummary)
    @test isa(gewekediag(chn[:,1,:]), MCMCChains.ChainSummary)
    @test isa(heideldiag(chn[:,1,:]), MCMCChains.ChainSummary)
    @test isa(rafterydiag(chn[:,1,:]), MCMCChains.ChainSummary)
end

@testset "concatenation tests" begin
    v1 = rand(500, 5, 1)
    v2 = rand(500, 5, 1)
    v3 = rand(500, 4, 1)
    v4 = rand(500, 5, 4)

    c1 = Chains(v1, [:a1, :a2, :a3, :a4, :a5])
    c2 = Chains(v2, [:a1, :a2, :a3, :a4, :a5], start = 501)
    c3 = Chains(v3, [:z1, :z2, :z3, :z4])
    c4 = Chains(v4, [:a1, :a2, :a3, :a4, :a5])

    # Test dim 1
    c1_2 = cat(c1, c2, dims = 1)
    @test c1_2.value.data == cat(v1, v2, dims=1)
    @test range(c1_2) == 1:1:1000
    @test names(c1_2) == names(c1) == names(c2)
    @test chains(c1_2) == chains(c1) == chains(c2)
    @test c1_2.value == vcat(c1, c2).value

    # Test dim 2
    c1_3 = cat(c1, c3, dims = 2)
    @test c1_3.value.data == cat(v1, v3, dims=2)
    @test range(c1_3) == 1:1:500
    @test names(c1_3) == cat(names(c1), names(c3), dims=1)
    @test chains(c1_3) == chains(c1) == chains(c3)
    @test c1_3.value == hcat(c1, c3).value

    # Test dim 3
    c1_4 = cat(c1, c4, dims = 3)
    @test c1_4.value.data == cat(v1, v4, dims=3)
    @test range(c1_4) == 1:1:500
    @test names(c1_4) == names(c1) == names(c4)
    @test length(chains(c1_4)) == length(chains(c1)) + length(chains(c4))
    @test c1_4.value == chainscat(c1, c4).value
end
