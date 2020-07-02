using MCMCChains
using Distributions

using Random
using Serialization
using Test

Random.seed!(20)

ProjDir = mktempdir()

@testset "serialization read and write test" begin
    # Consider the model
    # ```math
    # s ~ IG(2, 3)
    # m ~ Normal(0, √s)
    # x ~ Normal(m, √s)
    # ```
    # for variable ``x``. Given observations ``x₁ = 1.5`` and ``x₂ = 2``, sample one chain
    # with 500 samples from the posterior
    # ```math
    # m, s ~ Normal-IG(7/6, 3, 3, 49/12)
    # ```
    vals = Matrix{Float64}(undef, 500, 2)
    rand!(InverseGamma(3, 49/12), view(vals, :, 2))
    for i in 1:size(vals, 1)
        vals[i, 1] = rand(Normal(7/6, sqrt(vals[i, 2] / 3)))
    end
    chn1 = Chains(vals, ["m", "s"])

    # Julia 1.0 doesn't support `serialize(::AbstractString, value)`
    # and `deserialize(::AbstractString)`
    open(joinpath(ProjDir, "chn1.jls"), "w") do io
        serialize(io, chn1)
    end
    chn2 = open(deserialize, joinpath(ProjDir, "chn1.jls"), "r")

    open(joinpath(ProjDir, "chn1.txt"), "w") do io
        describe(io, chn1);
    end

    open(joinpath(ProjDir, "chn2.txt"), "w") do io
        describe(io, chn2);
    end

    @test open(f->read(f, String), joinpath(ProjDir, "chn1.txt")) ==
        open(f->read(f, String), joinpath(ProjDir, "chn2.txt"))
end

rm(ProjDir, force=true, recursive=true)
