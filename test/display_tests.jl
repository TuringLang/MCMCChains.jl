@testset "describe" begin
    @test describe(Chains(rand(10, 3, 1))) isa Any

    chn = Chains(rand(100, 2, 2), [:a, :b])
    function render(; context_digits = nothing, kwargs...)
        io = IOBuffer()
        output = context_digits === nothing ? io : IOContext(io, :digits => context_digits)
        describe(output, chn; kwargs...)
        return String(take!(io))
    end
    @test render(digits = 2) != render()
    @test render(context_digits = 2) == render(digits = 2)
    @test render(context_digits = 3, digits = 1) == render(context_digits = 1)
    @test render(digits = 0) == render(context_digits = 0)
    @test_throws ArgumentError describe(IOBuffer(), chn; digits = -1)
    for invalid in (true, false)
        io = IOBuffer()
        @test_throws ArgumentError describe(io, chn; digits = invalid)
        @test isempty(take!(io))
    end
end
