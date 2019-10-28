using Turing, Test

ProjDir = mktempdir()

@testset "serialization read and write test" begin
    @model gdemo(x) = begin
        s ~ InverseGamma(2,3)
        m ~ Normal(0, sqrt(s))
        for i in eachindex(x)
            x[i] ~ Normal(m, sqrt(s))
        end
    end

    model = gdemo([1.5, 2.0])
    sampler = HMC(0.01, 5)
    chn1 = sample(model, sampler, 500, save_state=true)

    write(joinpath(ProjDir, "chn1.jls"), chn1)
    chn2 = read(joinpath(ProjDir, "chn1.jls"), MCMCChains.Chains)

    open(joinpath(ProjDir, "chn1.txt"), "w") do io
        describe(io, chn1);
    end

    open(joinpath(ProjDir, "chn2.txt"), "w") do io
        describe(io, chn2);
    end

    @test open(f->read(f, String), joinpath(ProjDir, "chn1.txt")) ==
        open(f->read(f, String), joinpath(ProjDir, "chn2.txt"))

    # Test whether sampler state made it through serialization/deserialization.
    chn3 = Turing.Utilities.resume(chn2, 100)
    @test range(chn3) == 1:1:600
end

rm(ProjDir, force=true, recursive=true)
