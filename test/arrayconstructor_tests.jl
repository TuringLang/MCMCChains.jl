using MCMCChains, Test

@testset "Array constructor tests" begin
    @testset "Smoke tests" begin
        val = convert(Array{Union{Missing, Float64},3}, rand(1000, 5, 5))
        chns = Chains(val, ["a", "b", "c", "d", "e"], [:internals => ["d", "e"]])

        # Config.
        main_params = 3
        internal_params = 2
        total_params = main_params + internal_params

        d, p, c = size(chns)

        @test size(Array(chns)) == (d*c, main_params)
        @test size(Array(chns, [:parameters])) == (d*c, main_params)
        @test size(Array(chns, [:parameters])) == size(Array(chns))
        @test size(Array(chns, [:parameters, :internals])) == (d*c, total_params)
        @test size(Array(chns, [:parameters, :internals])) ==
            size(Array(chns, showall=true))
        @test size(Array(chns, [:internals])) == (d*c, internal_params)
        @test size(Array(chns, append_chains=true)) == (d*c, main_params)
        @test size(Array(chns, append_chains=false)) == (5,)
        @test size(Array(chns, append_chains=false)[1]) == (d, main_params)
        @test typeof(Array(chns, append_chains=true)) == Array{Float64, 2}
        @test size(Array(chns)) == (d*c, main_params)
        @test size(Array(chns, append_chains=false)) == (5,)
        @test size(Array(chns, append_chains=false)[1]) == (d, main_params)
        @test size(Array(chns[:b])) == (d*c,)

        Array(chns)
        Array(chns[:a])
        Array(chns, [:parameters])
        Array(chns, [:parameters, :internals])
    end
    @testset "Accuracy" begin
        nchains = 5
        val = rand(500, 5, nchains)
        long_val = vcat([val[:,:,i] for i in 1:nchains]...)
        squished_val = [val[:,:,i] for i in 1:nchains]
        chn = Chains(val)

        # Test no_append
        arr1 = Array(chn, append_chains=true)
        @test all(arr1 .== long_val)

        # Test append
        arr2 = Array(chn, append_chains=false)
        @test all(vcat([arr2[i] .== squished_val[i] for i in 1:nchains]...))
    end
end
