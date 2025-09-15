@testset "describe" begin
    @test describe(Chains(rand(10, 3, 1))) isa Any
end
