using MCMCChains
using AbstractMCMC: AbstractChains
using Dates
using DataFrames
using PosteriorStats: SummaryStats
using Tables
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
@test_throws ErrorException Chains(val; start=0, thin=2)
@test_throws ErrorException Chains(val; start=niter, thin=-1)
@test_throws ErrorException Chains(val; iterations=1:(niter - 1))
@test_throws ErrorException Chains(val; iterations=range(0; step=2, length=niter))
@test_throws ErrorException Chains(val; iterations=niter:-1:1)
@test_throws ErrorException Chains(val; iterations=ones(Int, niter))

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
    @test range(chn) == range(Chains(val; iterations=range(chn)))
    @test range(chn) == range(Chains(val; iterations=collect(range(chn))))

    @test_throws ErrorException setrange(chn, 1:10)
    @test_throws ErrorException setrange(chn, 0:(niter - 1))
    @test_throws ErrorException setrange(chn, niter:-1:1)
    @test_throws ErrorException setrange(chn, ones(Int, niter))
    @test_throws MethodError setrange(chn, float.(range(chn)))

    chn2a = setrange(chn, range(1; step = 10, length = niter))
    chn2b = setrange(chn, collect(range(1; step = 10, length = niter)))
    for chn2 in (chn2a, chn2b)
        @test range(chn2) == range(1; step = 10, length = niter)
        @test names(chn2) === names(chn)
        @test chains(chn2) === chains(chn)
        @test chn2.value.data === chn.value.data
        @test chn2.logevidence === chn.logevidence
        @test chn2.name_map === chn.name_map
        @test chn2.info == chn.info
    end

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
    @test mean(chn) isa SummaryStats
    @test mean(chn, ["param_1", "param_3"]) isa SummaryStats
    @test 0.95 ≤ mean(chn, "param_1") ≤ 1.05
end

@testset "Chain times" begin
    t1 = time()
    t2 = t1 + 1.5
    chn_timed = Chains(val, info = (start_time=t1, stop_time=t2))

    @test MCMCChains.max_stop(chn_timed) == unix2datetime(t2)
    @test MCMCChains.min_start(chn_timed) == unix2datetime(t1)
    @test MCMCChains.wall_duration(chn_timed) <= 1.6

    n1 = now()
    n2 = now() + Second(1)
    unix1 = datetime2unix(n1)
    unix2 = datetime2unix(n2)

    @test MCMCChains.min_datetime(n1) == n1
    @test MCMCChains.min_datetime([n1, n2]) == n1
    @test MCMCChains.min_datetime(unix1) == n1
    @test MCMCChains.min_datetime([unix1, unix2]) == n1
    @test MCMCChains.min_datetime("red") === missing

    @test MCMCChains.max_datetime(n1) == n1
    @test MCMCChains.max_datetime([n1, n2]) == n2
    @test MCMCChains.max_datetime(unix1) == n1
    @test MCMCChains.max_datetime([unix1, unix2]) == n2
    @test MCMCChains.max_datetime("red") === missing
end

@testset "indexing tests" begin
    c = chn[:, 1, :]
    @test c isa AbstractMatrix
    @test size(c) == (niter, nchains)
    @test c == val[:, 1, :]

    for i in 1:2
        c = chn[:, 1, i]
        @test c isa AbstractVector
        @test length(c) == niter
        @test c == val[:, 1, i]
    end

    for p in (:param_1, "param_1", SubString("param_1", 1))
        c = chn[200:300, p, :]
        @test c isa AbstractMatrix
        @test size(c) == (101, size(chn, 3))
        @test c == val[200:300, 1, :]

        c = chn[200:300, p, 1]
        @test c isa AbstractVector
        @test length(c) == 101
        @test c == val[200:300, 1, 1]
    end

    for ps in (
        [:param_1, :param_3],
        ["param_1", "param_3"],
        [SubString("param_1", 1), "param_3"],
        ["param_1", SubString("param_3", 1)],
        [SubString("param_1", 1), SubString("param_3", 1)],
    )
        c = chn[200:300, ps, :]
        @test c isa Chains
        @test size(c) == (101, 2, nchains)
        @test c.value.data == val[200:300, [1, 3], :]
    end
end

@testset "names and groups tests" begin
    chn2 = @inferred replacenames(chn, "param_2" => "param[2]", "param_3" => "param[3]")
    @test chn2.value ==
        (@inferred replacenames(chn, Dict("param_2" => "param[2]",
                                          "param_3" => "param[3]"))).value
    @test names(chn2) == [:param_1, Symbol("param[2]"), Symbol("param[3]"), :param_4]
    for p in (:param, "param", SubString("param", 1))
        @test namesingroup(chn2, p) == Symbol.(["param[2]", "param[3]"])
    end

    for p in (:param, "param", SubString("param", 1))
        chn3 = group(chn2, p)
        @test names(chn3) == Symbol.(["param[2]", "param[3]"])
        @test chn3.value == chn[:, [:param_2, :param_3], :].value
    end

    stan_chn = Chains(rand(100, 3, 1), ["a.1", "a[2]", "b"])
    for p in (:a, "a", SubString("a", 1))
        @test namesingroup(stan_chn, p; index_type=:dot) == [Symbol("a.1")]
        @test names(group(stan_chn, p; index_type=:dot)) == [Symbol("a.1")]
        @test_throws Exception namesingroup(stan_chn, p; index_type=:x)
        @test_throws Exception group(stan_chn, p; index_type=:x)    
    end
end

@testset "function tests" begin
    tchain = Chains(rand(niter, nparams, nchains), ["a", "b", "c"], Dict(:internals => ["c"]))

    @test eltype(discretediag(chn_disc[:, 2:2, :])) <: SummaryStats

    gelman = gelmandiag(tchain)
    gelmanmv = gelmandiag_multivariate(tchain)
    geweke = gewekediag(tchain)
    heidel = heideldiag(tchain)
    raferty = rafterydiag(tchain)

    # test raw return values
    @test typeof(gelman) <: SummaryStats
    @test typeof(gelmanmv) <: Tuple{SummaryStats,Float64}
    @test typeof(geweke) <: Array{<:SummaryStats}
    @test typeof(heidel) <: Array{<:SummaryStats}
    @test typeof(raferty) <: Array{<:SummaryStats}

    # test SummaryStats sizes
    for s in (gelman, gelmanmv[1], geweke[1], heidel[1], raferty[1])
        @test s isa SummaryStats
        df = DataFrame(s)
        @test size(df, 1) == 2
    end
    @test size(DataFrame(gelman), 2) == 3
    @test size(DataFrame(gelmanmv[1]), 2) == 3
    @test size(DataFrame(geweke[1]), 2) == 3
    @test size(DataFrame(heidel[1]), 2) == 7
    @test size(DataFrame(raferty[1]), 2) == 6
end

@testset "stats tests" begin
    # issue #363 (short chains)
    for c in (chn, chn[1:10, :, :])
        for append_chains in (true, false)
            # Number of samples for estimation
            carray = Array(c; append_chains=append_chains)
            n = append_chains ? size(carray, 1) : size(first(carray), 1)
            
            # Default lags
            lags = MCMCChains._default_lags(c, append_chains)
            @test lags == filter!(x -> x < n, [1, 5, 10, 50])

            acor = autocor(c; append_chains=append_chains)
            # Number of columns in the SummaryStats: lags + parameters
            ncols = length(lags) + 1
            if append_chains
                @test acor isa SummaryStats
                @test length(keys(acor)) == ncols
            else
                @test acor isa Vector{<:SummaryStats}
                @test all(length(keys(a)) == ncols for a in acor)
            end
        end
        @test autocor(c) isa SummaryStats
        @test autocor(c) == autocor(c; append_chains = true)
    end

    @test MCMCChains.cor(chn) isa SummaryStats
    @test MCMCChains.cor(chn; append_chains = false) isa Vector{<:SummaryStats}

    @test MCMCChains.changerate(chn) isa Tuple{SummaryStats,Float64}
    @test MCMCChains.changerate(chn; append_chains = false) isa
          Vector{<:Tuple{SummaryStats,Float64}}

    @test eti(chn) isa SummaryStats
    @test eti(chn; append_chains = false) isa Vector{<:SummaryStats}

    @test hdi(chn) isa SummaryStats
    @test hdi(chn; append_chains = false) isa Vector{<:SummaryStats}

    result = hdi(chn)
    @test :hdi94 in Tables.columnnames(result)

    @test_deprecated hpd(chn)
    @test hpd(chn) == hdi(chn; prob = 0.95)
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
