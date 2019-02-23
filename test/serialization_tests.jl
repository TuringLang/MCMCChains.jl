using  MCMCChains, Test

ProjDir = @__DIR__

@testset "serialization read and write test" begin
  
  m10_4s = read(joinpath(ProjDir, "sections_m10.4s.jls"), MCMCChains.Chains)

  write(joinpath(ProjDir, "sections_m10.4.1s.jls"), m10_4s)
  
  open(joinpath(ProjDir, "sections_m10.4s.txt"), "w") do io
    describe(io, m10_4s);
  end
  
  c = read(joinpath(ProjDir, "sections_m10.4.1s.jls"), MCMCChains.Chains)

  tmpdir = tempdir()
  open(joinpath(ProjDir, "sections_m10.4.1s.txt"), "w") do io
    describe(io, c);
  end

  @test open(f->read(f, String), joinpath(ProjDir, "sections_m10.4.1s.txt")) ==
      open(f->read(f, String), joinpath(ProjDir, "sections_m10.4s.txt"))
    
end
