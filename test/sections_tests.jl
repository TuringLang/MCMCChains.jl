using  MCMCChains, Test

@testset "describe sections" begin

    a3d = rand(500, 9, 4)

    cnames = ["lp__"  , "accept_stat__", "stepsize__" , "treedepth__" ,
        "n_leapfrog__" , "divergent__"  , "energy__", "sigma", "mu" ]

    pi = filter(p -> length(p) > 2 && p[end-1:end] == "__", cnames)
    p = filter(p -> !(p in  pi), cnames)

    c = Chains(a3d,
        cnames,
            Dict(
                :parameters => p,
                :internals => pi
                )
        )

    c2 = Chains(a3d,
        cnames,
            Dict(
                :parameters => p,
                :internals => pi
            )
        )

    tmpdir = tempdir()
    open(joinpath(tmpdir, "sections_test1.txt"), "w") do io
        describe(io, c, sections=:parameters);
        describe(io, c, sections=:internals);
        describe(io, c);
    end

    open(joinpath(tmpdir, "sections_test2.txt"), "w") do io
        describe(io, c2, sections=:parameters);
        describe(io, c2, sections=:internals);
        describe(io, c2);
    end


    @test open(f->read(f, String), joinpath(tmpdir, "sections_test1.txt")) ==
        open(f->read(f, String), joinpath(tmpdir, "sections_test2.txt"))
end
