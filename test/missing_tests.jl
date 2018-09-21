using MCMCChain
using Test

# tests for missing values

@testset "utils" begin
    x = [1, missing, 3, 2]
    @test MCMCChain.cummean(x) == [1., 1., 2., 2.]
end

@testset "stats" begin
    chn = Chains(randn(1000, 2, 2))

    # call describe without missing values
    (s1, s2) = describe(devnull, chn)

    # add missing values
    #chn.value[2,:,:] .= missing
    chn = Chains(cat(chn.value, ones(1, 2, 2) .* missing, dims = 1))

    # call describe with missing values
    (m1, m2) = describe(devnull, chn)

    @test all(s1.value[:,1:4,:] .== m1.value[:,1:4,:])
    @test all(s1.value[:,5,:] .+ 1 .== m1.value[:,5,:])
end
