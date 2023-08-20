using MCMCChains
using StatsBase
using Test
using Random

# Tests for missing values.
function testdiff(cdf1, cdf2)
    return all(((x, y),) -> isapprox(x, y; atol=1e-2), Iterators.drop(zip(cdf1, cdf2), 1))
end

@testset "utils" begin
    x = [1, missing, 3, 2]
    @test MCMCChains.cummean(x) == [1., 1., 2., 2.]
end

@testset "diagnostic functions" begin
    Random.seed!(1234)
    nchains = 2
    chn = Chains(randn(5000, 2, nchains))

    # Add missing values.
    chn_m = Chains(cat(chn.value, ones(1, 2, nchains) .* missing, dims = 1))

    # Currently we only check if diagnostic functions throw a `MethodError` if the
    # element type of the chain is not a subtype of `Real`.
    @test_throws MethodError gelmandiag(chn_m)

    gw_1 = gewekediag(chn)
    gw_2 = gewekediag(chn_m)
    hd_1 = heideldiag(chn)
    hd_2 = heideldiag(chn_m)
    rf_1 = rafterydiag(chn)
    rf_2 = rafterydiag(chn_m)

    @testset "diagnostics missing tests" for i in 1:nchains
        @test all(Base.splat(testdiff), zip(gw_1, gw_2))
        @test all(Base.splat(testdiff), zip(hd_1, hd_2))
        @test all(Base.splat(testdiff), zip(rf_1, rf_2))
    end

    @test_throws MethodError discretediag(chn_m)
end
