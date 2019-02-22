using  MCMCChains, Test

ProjDir = @__DIR__

@testset "serialization write test" begin
  
  m10_4s = read(joinpath(@__DIR__, "sections_m10.4s.jls"), MCMCChains.Chains)

  cnames =  ["lp__", "accept_stat__", "stepsize__" , "treedepth__",
  "n_leapfrog__" , "divergent__", "energy__",
  "a_1", "a_2", "a_3", "a_4", "a_5", "a_6", "a_7",        
  "bp", "bpC"]        
  
  pi = filter(p -> length(p) > 2 && p[end-1:end] == "__", cnames)
  p = filter(p -> !(p in  pi), cnames)
  p1 = ["bp", "bpC"]
  p2 = filter(p -> !(p in p1), p)

  c = Chains(m10_4s.value[:,:,:],
    Symbol.(cnames),
    Dict(
      :parameters => Symbol.(p1),
      :pooled => Symbol.(p2),
      :internals => Symbol.(pi)
    )
  )
  
  write(joinpath(@__DIR__, "sections_m10.4.1s.jls"), c)
end


@testset "serialization read test" begin
  
  c = read(joinpath(ProjDir, "sections_m10.4.1s.jls"), MCMCChains.Chains)

  tmpdir = tempdir()
  open(joinpath(tmpdir, "sections_m10.4.1s.txt"), "w") do io
    describe(io, c, section=:pooled);
  end

  @test open(f->read(f, String), joinpath(tmpdir, "sections_m10.4.1s.txt")) ==
      open(f->read(f, String), joinpath(ProjDir, "sections_m10.4s.txt"))
    
end
