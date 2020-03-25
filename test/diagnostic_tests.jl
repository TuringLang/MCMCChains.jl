using MCMCChains
using AbstractMCMC: AbstractChains
using Test

## CHAIN TESTS
# Define the experiment
n_iter = 4000
n_name = 3
n_chain = 2

# some sample experiment results
val = randn(n_iter, n_name, n_chain) .+ [1, 2, 3]'
val = hcat(val, rand(1:2, n_iter, 1, n_chain))

# construct a Chains object
chn = Chains(val, start = 1, thin = 2)

# Chains object for discretediag
val_disc = rand(Int16, 200, n_name, n_chain)
chn_disc = Chains(val_disc, start = 1, thin = 2)

@testset "basic chains functions" begin
    @test first(chn) == 1
    @test step(chn) == 2
    @test last(chn) == 7999
    @test size(chn) == (4000, 4, 2)
    @test size(chn[1:1000, :, :]) == (1000, 4, 2)
    @test keys(chn) == [:param1, :param2, :param3, :param4]
    @test isa(chn[:,1,:], AbstractChains)
    @test isa(chn[200:300, "param1", :], AbstractChains)
    @test isa(chn[200:300, ["param1", "param3"], :], AbstractChains)
    @test isa(chn[200:300, "param1", 1], AbstractChains)
    @test length(vec(chn[:,1,:].value)) == n_chain * n_iter
    @test all(collect(skipmissing(chn[:,1,1].value)) .== val[:,1,1])
    @test all(chn[:,1,2].value .== val[:,1,2])
    @test all(MCMCChains.indiscretesupport(chn) .== [false, false, false, true])
    @test setinfo(chn, NamedTuple{(:A, :B)}((1,2))).info == NamedTuple{(:A, :B)}((1,2))
    @test isa(set_section(chn, Dict(:internals => ["param1"])), AbstractChains)
    @test mean(chn) isa ChainDataFrame
    chn2 = @inferred replacenames(chn, "param1" => "PARAM1")
    @test names(chn2) == [:PARAM1, :param2, :param3, :param4]
end

@testset "function tests" begin
    # the following tests only check if the function calls work!
    @test MCMCChains.diag_all(rand(50, 2), :weiss, 1, 1, 1) != nothing
    @test MCMCChains.diag_all(rand(50, 2), :hangartner, 1, 1, 1) != nothing
    @test MCMCChains.diag_all(rand(50, 2), :billingsley, 1, 1, 1) != nothing

    @test eltype(discretediag(chn_disc[:,2,:])) <: ChainDataFrame
    @test typeof(gelmandiag(chn[:,1,:])) <: ChainDataFrame
    @test eltype(gewekediag(chn[:,1,:])) <: ChainDataFrame
    @test eltype(heideldiag(chn[:,1,:])) <: ChainDataFrame
    @test eltype(rafterydiag(chn[:,1,:])) <: ChainDataFrame
end

@testset "stats tests" begin
    @test autocor(chn) isa Vector{<:ChainDataFrame}
    @test autocor(chn; append_chains = false) isa Vector{<:ChainDataFrame}

    @test MCMCChains.cor(chn) isa ChainDataFrame
    @test MCMCChains.cor(chn; append_chains = false) isa Vector{<:ChainDataFrame}

    @test MCMCChains.changerate(chn) isa ChainDataFrame
    @test MCMCChains.changerate(chn; append_chains = false) isa Vector{<:ChainDataFrame}

    @test hpd(chn) isa ChainDataFrame
    @test hpd(chn; append_chains = false) isa Vector{<:ChainDataFrame}
end

@testset "vector of vectors" begin
    val = [rand(20) for _ in 1:10]

    chn = Chains(val)
    chn2 = Chains(reduce(hcat, val)')

    @test chains(chn) == chains(chn2)
    @test names(chn) == names(chn2)
end

@testset "sorting" begin
    chn_unsorted = Chains(rand(100,3,1), ["2", "1", "3"])
    chn_sorted = sort(chn_unsorted)

    @test names(chn_sorted) == Symbol.([1, 2, 3])
    @test names(chn_unsorted) == Symbol.([2, 1, 3])
end
