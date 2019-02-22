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
    s1, s2 = MCMCChains.summarystats(chn).summaries[1], MCMCChains.quantile(chn).summaries[1]

    # Add missing values.
    chn = Chains(cat(chn.value, ones(1, 2, 2) .* missing, dims = 1))

    # Call describe with missing values.
    describe(devnull, chn; showall=true)
    m1, m2 = MCMCChains.summarystats(chn).summaries[1], MCMCChains.quantile(chn).summaries[1]

    @test all(s1.value[:,1:4,:] .== m1.value[:,1:4,:])
    @test all(s1.value[:,5,:] .+ 1 .== m1.value[:,5,:])
end

@testset "diagnostic functions" begin
    chn = Chains(randn(5000, 2, 2))

    # Add missing values.
    chn_m = Chains(cat(chn.value, ones(1, 2, 2) .* missing, dims = 1))

    # Currently we only check if diag. functions throw an assertion if a value is missing.
    @test_throws AssertionError gelmandiag(chn_m)

    @test all(gewekediag(chn).value .== gewekediag(chn_m).value)
    @test all(heideldiag(chn).value .== heideldiag(chn_m).value)
    @test all(rafterydiag(chn).value .== rafterydiag(chn_m).value)

    @test_throws AssertionError discretediag(chn_m)
end
