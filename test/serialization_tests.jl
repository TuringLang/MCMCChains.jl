using  MCMCChains, Test

ProjDir = mktempdir()

@testset "serialization read and write test" begin
    val = rand(500, 5, 1)
    chn1 = Chains(val)

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
end

rm(ProjDir, force=true, recursive=true)
