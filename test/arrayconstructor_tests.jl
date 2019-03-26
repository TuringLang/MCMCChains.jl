using Turing, MCMCChains, Test

@testset "Array constructor tests" begin
  
  @model gdemo(x) = begin
      m ~ Normal(1, 0.01)
      s ~ Normal(5, 0.01)
  end

  model = gdemo([1.5, 2.0])
  sampler = HMC(1000, 0.01, 5)

  chns = [sample(model, sampler, save_state=true) for i in 1:4]
  chns = chainscat(chns...) 
  
  d, p, c = size(chns.value.data)

  @test size(Array(chns)) == (d*c, p)
  @test size(Array(chns, [:parameters])) == (d*c, 2)
  @test size(Array(chns, [:parameters, :internals])) == (d*c, p)
  @test size(Array(chns, [:internals])) == (d*c, 6)
  @test size(Array(chns, append_chains=true)) == (d*c, p)
  @test size(Array(chns, append_chains=false)) == (4,)
  @test size(Array(chns, append_chains=false)[1]) == (d, p)
  @test typeof(Array(chns, append_chains=true)) == Array{Float64, 2}
  @test size(Array(chns, remove_missing_union=false)) == (d*c, p)
  @test size(Array(chns, append_chains=false, remove_missing_union=false)) == (4,)
  @test size(Array(chns, append_chains=false, remove_missing_union=false)[1]) == (d, p)
  @test typeof(Array(chns, append_chains=true, remove_missing_union=false)) == 
    Array{Union{Missing, Float64}, 2}
  @test size(Array(chns[:m])) == (d*c,)

  Array(chns)
  Array(chns[:s])
  Array(chns, [:parameters])
  Array(chns, [:parameters, :internals])
  
end