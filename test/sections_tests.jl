using  MCMCChains, Test

@testset "describe sections" begin
  
  a3d = reshape(1:3600, 100, 9, 4)

  cnames = ["lp__"  , "accept_stat__", "stepsize__" , "treedepth__" ,
    "n_leapfrog__" , "divergent__"  , "energy__", "sigma", "mu" ]

  pi = filter(p -> length(p) > 2 && p[end-1:end] == "__", cnames)
  p = filter(p -> !(p in  pi), cnames)

  c = Chains(a3d,
    Symbol.(cnames),
    Dict(
      :parameters => Symbol.(p),
      :internals => Symbol.(pi)
    )
  )
  
  tmpdir = tempdir()
  open(joinpath(tmpdir, "sections_test.txt"), "w") do io
    describe(io, c, section=:parameters);
    describe(io, c, section=:internals);
    describe(io, c);
  end

  
  @test open(f->read(f, String), joinpath(tmpdir, "sections_test.txt")) ==
      open(f->read(f, String), joinpath(@__DIR__, "sections_test.txt"))
    
  #println(open(f->read(f, String), joinpath(tmpdir, "sections_test.txt")))
  
end
