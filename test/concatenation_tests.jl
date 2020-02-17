using MCMCChains
using Test

@testset "vcat" begin
    chn = Chains(rand(10, 5, 2), ["a", "b", "c", "d", "e"], Dict(:internal => ["d", "e"]))
    chn1 = Chains(rand(5, 5, 2), ["a", "b", "c", "d", "e"], Dict(:internal => ["a", "b"]))

    # incorrect thinning
    @test_throws ArgumentError vcat(chn, Chains(rand(2, 5, 2); thin = 2))

    # incorrect names
    @test_throws ArgumentError vcat(chn, Chains(rand(10, 5, 2), ["a", "b", "c", "d", "f"]))

    # incorrect number of chains
    @test_throws ArgumentError vcat(chn, Chains(rand(10, 5, 3), ["a", "b", "c", "d", "e"]))

    # concate the same chain
    chn2 = vcat(chn, chn)
    @test chn2.value.data == vcat(chn.value.data, chn.value.data)
    @test size(chn2) == (20, 5, 2)
    @test names(chn2) == names(chn)
    @test range(chn2) == 1:20
    @test chn2.name_map == (internal = ["d", "e"], parameters = ["c", "b", "a"])
    
    chn2a = cat(chn, chn)
    @test chn2a.value == chn2.value
    @test chn2a.name_map == chn2.name_map
    @test chn2a.info == chn2.info

    chn2b = cat(chn, chn; dims = Val(1))
    @test chn2b.value == chn2.value
    @test chn2b.name_map == chn2.name_map
    @test chn2b.info == chn2.info

    chn2c = cat(chn, chn; dims = 1)
    @test chn2c.value == chn2.value
    @test chn2c.name_map == chn2.name_map
    @test chn2c.info == chn2.info

    # concatenate a different chain
    chn3 = vcat(chn, chn1)
    @test chn3.value.data == vcat(chn.value.data, chn1.value.data)
    @test size(chn3) == (15, 5, 2)
    @test names(chn3) == names(chn)
    @test range(chn3) == 1:15
    @test chn3.name_map == (internal = ["a", "b"], parameters = ["c", "e", "d"])
    
    chn3a = cat(chn, chn1)
    @test chn3a.value == chn3.value
    @test chn3a.name_map == chn3.name_map
    @test chn3a.info == chn3.info

    chn3b = cat(chn, chn1; dims = Val(1))
    @test chn3b.value == chn3.value
    @test chn3b.name_map == chn3.name_map
    @test chn3b.info == chn3.info

    chn3c = cat(chn, chn1; dims = 1)
    @test chn3c.value == chn3.value
    @test chn3c.name_map == chn3.name_map
    @test chn3c.info == chn3.info
end

@testset "hcat" begin
    chn = Chains(rand(10, 3, 2), ["a", "b", "c"], Dict(:internal => ["c"]))
    chn1 = Chains(rand(10, 2, 2), ["d", "e"], Dict(:internal => ["d"]))

    # incorrect ranges
    @test_throws ArgumentError hcat(chn, Chains(rand(10, 2, 2); start = 2))
    @test_throws ArgumentError hcat(chn, Chains(rand(10, 2, 2); thin = 2))

    # non-unique names
    @test_throws ArgumentError hcat(chn, chn)

    # incorrect number of chains
    @test_throws ArgumentError hcat(chn, Chains(rand(10, 2, 3)))

    # concatenate two chains
    chn2 = hcat(chn, chn1)
    @test chn2.value.data == hcat(chn.value.data, chn1.value.data)
    @test size(chn2) == (10, 5, 2)
    @test names(chn2) == vcat(names(chn), names(chn1))
    @test range(chn2) == 1:10
    @test chn2.name_map == (internal = ["d"], parameters = ["e", "c", "b", "a"])
    
    chn2a = cat(chn, chn1; dims = Val(2))
    @test chn2a.value == chn2.value
    @test chn2a.name_map == chn2.name_map
    @test chn2a.info == chn2.info

    chn2b = cat(chn, chn1; dims = 2)
    @test chn2b.value == chn2.value
    @test chn2b.name_map == chn2.name_map
    @test chn2b.info == chn2.info
end

@testset "chainscat" begin
    chn = Chains(rand(10, 3, 3), ["a", "b", "c"], Dict(:internal => ["c"]))
    chn1 = Chains(rand(10, 3, 2), ["a", "b", "c"], Dict(:internal => ["a"]))

    # incorrect ranges
    @test_throws ArgumentError chainscat(chn, Chains(rand(10, 3, 2); start = 2))
    @test_throws ArgumentError chainscat(chn, Chains(rand(10, 3, 2); thin = 2))

    # incorrect names
    @test_throws ArgumentError chainscat(chn, Chains(rand(10, 3, 2), ["d", "e", "f"]))

    # concatenate two chains
    chn2 = chainscat(chn, chn1)
    @test chn2.value.data == cat(chn.value.data, chn1.value.data; dims = 3)
    @test size(chn2) == (10, 3, 5)
    @test names(chn2) == names(chn)
    @test range(chn2) == 1:10
    @test chn2.name_map == (internal = ["a"], parameters = ["c", "b"])
    
    chn2a = cat(chn, chn1; dims = Val(3))
    @test chn2a.value == chn2.value
    @test chn2a.name_map == chn2.name_map
    @test chn2a.info == chn2.info

    chn2b = cat(chn, chn1; dims = 3)
    @test chn2b.value == chn2.value
    @test chn2b.name_map == chn2.name_map
    @test chn2b.info == chn2.info
end