using Turing, MCMCChains, KernelDensity, StatsBase, Test, Statistics

#@testset "sampling api" begin

  @model gdemo(x) = begin
      m ~ Normal(1, 0.01)
      s ~ Normal(5, 0.01)
  end

  model = gdemo([1.5, 2.0])
  sampler = HMC(500, 0.01, 5)
  chn = sample(model, sampler, save_state=true);

  chn_sample = sample(chn, 5)
  @test range(chn_sample) == 1:1:5

  c = kde(Array(chn[:s]))
  chn_weighted_sample = sample(c.x, Weights(c.density), 100000)

  @test mean(Array(chn[:s])) ≈ 5.0 atol=0.1

  #end
