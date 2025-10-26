@testset "describe / show" begin
    c = Chains(rand(10, 3, 2))
    @test describe(c) isa Any
    @test show(c) isa Any
    @test show(mean(c)) isa Any
end
