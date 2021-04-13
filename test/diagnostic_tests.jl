using MCMCChains
using AbstractMCMC: AbstractChains
using Dates
using Test

## CHAIN TESTS
# Define the experiment
niter = 4000
nparams = 3
nchains = 2

# some sample experiment results
val = randn(niter, nparams, nchains) .+ [1, 2, 3]'
val = hcat(val, rand(1:2, niter, 1, nchains))

# construct a Chains object
chn = Chains(val, start = 1, thin = 2)

# Chains object for discretediag
val_disc = rand(Int16, 200, nparams, nchains)
chn_disc = Chains(val_disc, start = 1, thin = 2)

@testset "basic chains functions" begin
    @test first(chn) == 1
    @test step(chn) == 2
    @test last(chn) == 7999
    @test size(chn) == (4000, 4, 2)
    @test size(chn[1:1000, :, :]) == (1000, 4, 2)
    @test keys(chn) == names(chn) == [:param_1, :param_2, :param_3, :param_4]

    @test range(chn) == range(1; step = 2, length = niter)

    @test_throws ErrorException setrange(chn, 1:10)
    @test_throws MethodError setrange(chn, float.(range(chn)))

    chn2 = setrange(chn, range(1; step = 10, length = niter))
    @test range(chn2) == range(1; step = 10, length = niter)
    @test names(chn2) === names(chn)
    @test chains(chn2) === chains(chn)
    @test chn2.value.data === chn.value.data
    @test chn2.logevidence === chn.logevidence
    @test chn2.name_map === chn.name_map
    @test chn2.info == chn.info

    chn3 = resetrange(chn)
    @test range(chn3) == 1:niter
    @test names(chn3) === names(chn)
    @test chains(chn3) === chains(chn)
    @test chn3.value.data === chn.value.data
    @test chn3.logevidence === chn.logevidence
    @test chn3.name_map === chn.name_map
    @test chn3.info == chn.info

    @test all(MCMCChains.indiscretesupport(chn) .== [false, false, false, true])
    @test setinfo(chn, NamedTuple{(:A, :B)}((1,2))).info == NamedTuple{(:A, :B)}((1,2))
    @test isa(set_section(chn, Dict(:internals => ["param_1"])), AbstractChains)
    @test mean(chn) isa ChainDataFrame
    @test mean(chn, ["param_1", "param_3"]) isa ChainDataFrame
    @test 0.95 ≤ mean(chn, "param_1") ≤ 1.05
end

@testset "Chain times" begin
    t1 = time()
    t2 = t1 + 1.5
    chn_timed = Chains(val, info = (start_time=t1, stop_time=t2))

    @test MCMCChains.max_stop(chn_timed) == unix2datetime(t2)
    @test MCMCChains.min_start(chn_timed) == unix2datetime(t1)
    @test MCMCChains.wall_duration(chn_timed) <= 1.6
end

@testset "indexing tests" begin
    @test chn[:,1,:] isa AbstractMatrix
    @test chn[200:300, "param_1", :] isa AbstractMatrix
    @test chn[200:300, ["param_1", "param_3"], :] isa Chains
    @test chn[200:300, "param_1", 1] isa AbstractVector
    @test size(chn[:,1,:]) == (niter, nchains)
    @test chn[:,1,1] == val[:,1,1]
    @test chn[:,1,2] == val[:,1,2]
end

@testset "names and groups tests" begin
    chn2 = @inferred replacenames(chn, "param_2" => "param[2]", "param_3" => "param[3]")
    @test chn2.value ==
        (@inferred replacenames(chn, Dict("param_2" => "param[2]",
                                          "param_3" => "param[3]"))).value
    @test names(chn2) == [:param_1, Symbol("param[2]"), Symbol("param[3]"), :param_4]
    @test namesingroup(chn2, "param") == Symbol.(["param[2]", "param[3]"])

    chn3 = group(chn2, "param")
    @test names(chn3) == Symbol.(["param[2]", "param[3]"])
    @test chn3.value == chn[:, [:param_2, :param_3], :].value
end

@testset "function tests" begin
    tchain = Chains(rand(niter, nparams, nchains), ["a", "b", "c"], Dict(:internals => ["c"]))

    # the following tests only check if the function calls work!
    @test MCMCChains.diag_all(rand(50, 2), :weiss, 1, 1, 1) != nothing
    @test MCMCChains.diag_all(rand(50, 2), :hangartner, 1, 1, 1) != nothing
    @test MCMCChains.diag_all(rand(50, 2), :billingsley, 1, 1, 1) != nothing

    @test eltype(discretediag(chn_disc[:,2:2,:])) <: ChainDataFrame

    gelman = gelmandiag(tchain)
    geweke = gewekediag(tchain)
    heidel = heideldiag(tchain)
    raferty = rafterydiag(tchain)

    # test raw return values
    @test typeof(gelman) <: ChainDataFrame
    @test typeof(geweke) <: Array{<:ChainDataFrame}
    @test typeof(heidel) <: Array{<:ChainDataFrame}
    @test typeof(raferty) <: Array{<:ChainDataFrame}

    # test ChainDataFrame sizes
    @test size(gelman) == (2,3)
    @test size(geweke[1]) == (2,3)
    @test size(heidel[1]) == (2,7)
    @test size(raferty[1]) == (2,6)
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

    result = hpd(chn)
    @test all(result.nt.upper .> result.nt.lower)
end

@testset "vector of vectors" begin
    val = [rand(20) for _ in 1:10]

    chn = Chains(val)
    chn2 = Chains(reduce(hcat, val)')

    @test chains(chn) == chains(chn2)
    @test names(chn) == names(chn2)
end

@testset "sorting" begin
    chn_unsorted = Chains(rand(100, nparams, 1), ["2", "1", "3"])
    chn_sorted = sort(chn_unsorted)

    @test names(chn_sorted) == Symbol.([1, 2, 3])
    @test names(chn_unsorted) == Symbol.([2, 1, 3])
end
