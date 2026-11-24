import AbstractMCMC

@testset "SamplingChains" begin
    if pkgversion(AbstractMCMC) >= v"5.18.0"
        draws = [(x = [i, j],) for i = 1:5, j = 1:2]
        output = AbstractMCMC.SamplingOutput(
            draws;
            iterations = 3:2:11,
            sampler_states = [[1], [2]],
        )
        chain = SamplingChains(output)
        selected = chain[[begin, 2, end], 2]

        @test chain.output === output
        @test size(selected) == (3, 1)
        @test range(selected) == [3, 5, 11]
        @test selected[2, 1] === draws[2, 2]
        @test parent(selected.output.samples) === draws
        @test selected.output.sampler_states[1] === output.sampler_states[2]
    else
        @test !isdefined(MCMCChains, :SamplingChains)
    end
end
