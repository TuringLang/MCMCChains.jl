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
    @test length(vec(chn[:,1,:].value)) == n_chain * n_iter
    @test all(collect(skipmissing(chn[:,1,1].value)) .== val[:,1,1])
    @test all(chn[:,1,2].value .== val[:,1,2])
    @test all(MCMCChains.indiscretesupport(chn) .== [false, false, false, true])
end

@testset "function tests" begin
    # do not test the following functions
    # - wrtsp
    # - window2inds (tested above, see getindex)
    # -

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
