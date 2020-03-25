using MCMCChains, Test

@testset "Array constructor tests" begin
    @testset "Smoke tests" begin
        val = convert(Array{Union{Missing,Real},3}, rand(1000, 5, 5))
        chns = Chains(val, ["a", "b", "c", "d", "e"], [:internals => ["d", "e"]])
        chns_a = chns[:a]

        # Config.
        main_params = 3
        internal_params = 2
        total_params = main_params + internal_params

        d, p, c = size(chns)

        # return type
        @test Array(chns) isa Matrix{Float64}
        @test Array(MCMCChains.concretize(chns)) isa Matrix{Float64}
        @test MCMCChains.concretize(Array(chns)) isa Matrix{Float64}

        @test Array(chns; append_chains = true) isa Matrix{Float64}
        @test Array(MCMCChains.concretize(chns); append_chains = true) isa Matrix{Float64}
        @test MCMCChains.concretize(Array(chns; append_chains = true)) isa Matrix{Float64}

        @test Array(chns; append_chains = false) isa Vector{Matrix{Float64}}
        @test Array(MCMCChains.concretize(chns); append_chains = false) isa
            Vector{Matrix{Float64}}
        @test MCMCChains.concretize(Array(chns; append_chains = false)) isa
            Vector{Matrix{Float64}}

        @test Array(chns_a; append_chains = false) isa Vector{Vector{Float64}}
        @test Array(MCMCChains.concretize(chns_a); append_chains = false) isa
            Vector{Vector{Float64}}
        @test MCMCChains.concretize(Array(chns_a; append_chains = false)) isa
            Vector{Vector{Float64}}

        @test Array(chns; remove_missing_union = false) isa Matrix{Union{Missing,Real}}
        @test Array(chns; append_chains = true, remove_missing_union = false) isa
            Matrix{Union{Missing,Real}}
        @test Array(chns; append_chains = false, remove_missing_union = false) isa
            Vector{Matrix{Union{Missing,Real}}}
        @test Array(chns_a; append_chains = false, remove_missing_union = false) isa
            Vector{Vector{Union{Missing,Real}}}

        # type inference (needs concretely typed Chains)
        @inferred MCMCChains.to_matrix(chns)
        @inferred MCMCChains.to_vector(chns)
        @inferred MCMCChains.to_vector_of_vectors(chns_a)
        @inferred MCMCChains.to_vector_of_matrices(chns)

        # sizes
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

    @testset "concretize" begin
        z = Any[0.5 1; 0.5 missing]
        @test eltype(MCMCChains.concretize(z)) === Union{Float64,Missing}

        zz = [Any[0.5, 1, missing] for _ in 1:3]
        @test eltype(MCMCChains.concretize(zz)) === Vector{Union{Float64,Missing}}

        zzz = Any[Any[0.5, 1, missing] for _ in 1:3]
        @test eltype(MCMCChains.concretize(zzz)) === Vector{Union{Float64,Missing}}
    end
end
