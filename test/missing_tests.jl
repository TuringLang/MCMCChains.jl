using MCMCChains
using StatsBase
using Test
using Random

# Tests for missing values.
function testdiff(cdf1, cdf2)
    m1 = convert(Array, cdf1)
    m2 = convert(Array, cdf2)
    return all(((x, y),) -> isapprox(x, y; atol=1e-2), zip(m1, m2))
end

@testset "utils" begin
    x = [1, missing, 3, 2]
    @test MCMCChains.cummean(x) == [1., 1., 2., 2.]
end

# for now, summarystats will return missing whenever the chain contains one or more missing values
# @testset "stats" begin
#     Random.seed!(1234)
#     chn = Chains(randn(1000, 2, 2))

#     # Call describe without missing values.
#     describe(devnull, chn; sections = nothing)
#     s1, s2 = summarystats(chn), quantile(chn)

#     # Add missing values.
#     chn_m = Chains(cat(chn.value, ones(1, 2, 2) .* missing, dims = 1))

#     # Call describe with missing values.
#     describe(devnull, chn_m; sections = nothing)
#     m1, m2 = summarystats(chn_m), quantile(chn_m)

#     @test testdiff(s1, m1)
# end

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
        @test testdiff(gw_1, gw_2)
        @test testdiff(hd_1, hd_2)
        @test testdiff(rf_1, rf_2)
    end

    @test_throws MethodError discretediag(chn_m)
end
