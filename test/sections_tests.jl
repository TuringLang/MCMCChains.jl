using  MCMCChains, Test

# https://github.com/TuringLang/AdvancedMH.jl/pull/63
# https://github.com/TuringLang/MCMCChains.jl/issues/443
@testset "order of parameters" begin
    params = vcat(["μ[$i]" for i in 1:9], :p1, :p2, :p3)
    chains = Chains(rand(2, length(params)+1, 4), vcat(params, :lp), (internals=[:lp],))
    @test names(chains, :parameters) == map(Symbol, params)
    @test collect(keys(get(chains; section=:parameters))) == [:μ, :p1, :p2, :p3]
end

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
