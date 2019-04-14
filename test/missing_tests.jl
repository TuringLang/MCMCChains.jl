using MCMCChains
using Test

# Tests for missing values.

@testset "utils" begin
    x = [1, missing, 3, 2]
    @test MCMCChains.cummean(x) == [1., 1., 2., 2.]
end

@testset "stats" begin
    chn = Chains(randn(1000, 2, 2))

    # Call describe without missing values.
    describe(devnull, chn; showall=true)
    s1, s2 = MCMCChains.summarystats(chn),
      MCMCChains.quantile(chn)

    # Add missing values.
    chn_m = Chains(cat(chn.value, ones(1, 2, 2) .* missing, dims = 1))

    # Call describe with missing values.
    describe(devnull, chn_m; showall=true)
    m1, m2 = MCMCChains.summarystats(chn),
      MCMCChains.quantile(chn)

    @test s1[:,2:4] == m1[:,2:4]
end

@testset "diagnostic functions" begin
    nchains = 2
    chn = Chains(randn(5000, 2, nchains))

    # Add missing values.
    chn_m = Chains(cat(chn.value, ones(1, 2, nchains) .* missing, dims = 1))

    # Currently we only check if diag. functions throw an assertion if a value is missing.
    @test_throws AssertionError gelmandiag(chn_m)

    gw_1 = gewekediag(chn)
    gw_2 = gewekediag(chn_m)
    hd_1 = heideldiag(chn)
    hd_2 = heideldiag(chn_m)
    rf_1 = rafterydiag(chn)
    rf_2 = rafterydiag(chn_m)

    @testset "diagnostics missing tests" for i in 1:nchains
        @test all(gw_1[i].df == gw_2[i].df)
        @test all(hd_1[i].df == hd_2[i].df)
        @test all(rf_1[i].df == rf_2[i].df)
    end

    @test_throws AssertionError discretediag(chn_m)
end
