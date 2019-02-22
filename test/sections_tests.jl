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

@testset "describe multilevel model sections" begin
  
  cnames =  ["lp__", "accept_stat__", "stepsize__" , "treedepth__",
  "n_leapfrog__" , "divergent__", "energy__",
  "a_1", "a_2", "a_3", "a_4", "a_5", "a_6", "a_7",        
  "bp", "bpC"]        
  
  # Below section is only used to create a valid Chains object.
  # The object is written (serialized) into `test/sections_m10.4s.jls`
  # Note that if the input chain data (a3d) changes, the file
  # `sections_m10.4s.txt` needs to be updated.
  
  #=
  pi = filter(p -> length(p) > 2 && p[end-1:end] == "__", cnames)
  p = filter(p -> !(p in  pi), cnames)

  c = Chains(a3d,
    Symbol.(cnames),
    Dict(
      :parameters => Symbol.(p),
      :internals => Symbol.(pi)
    )
  )
  
  write(joinpath(@__DIR__, "sections_m10.4s.jls"), c)
  =#
  
  global m10_4s
  m10_4s = read(joinpath(@__DIR__, "sections_m10.4s.jls"), MCMCChains.Chains)
  
  #describe(m10_4s)
  
  pi = filter(p -> length(p) > 2 && p[end-1:end] == "__", cnames)
  p = filter(p -> !(p in  pi), cnames)
  p1 = ["bp", "bpC"]
  p2 = filter(p -> !(p in p1), p)

  c1 = Chains(m10_4s.value[:,:,:],
    Symbol.(cnames),
    Dict(
      :parameters => Symbol.(p1),
      :pooled => Symbol.(p2),
      :internals => Symbol.(pi)
    )
  )
  
  tmpdir = tempdir()
  open(joinpath(tmpdir, "sections_m10.4s.txt"), "w") do io
    describe(io, c1, section=:pooled);
  end

  @test open(f->read(f, String), joinpath(tmpdir, "sections_m10.4s.txt")) ==
      open(f->read(f, String), joinpath(@__DIR__, "sections_m10.4s.txt"))
    
  #println(open(f->read(f, String), joinpath(tmpdir,   "sections_m10.4s.txt"))
  
end
