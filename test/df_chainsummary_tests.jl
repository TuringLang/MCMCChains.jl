using Turing, MCMCChains, StatsBase

#@testset "Dataframe summary" begin
  
  @model gdemo(x) = begin
      m ~ Normal(1, 0.01)
      s ~ Normal(5, 0.01)
  end

  model = gdemo([1.5, 2.0])
  sampler = HMC(1000, 0.01, 5)

  chns = [sample(model, sampler, save_state=true) for i in 1:4]
  chns = chainscat(chns...) 
  
  d, p, c = size(chns.value.data)


  describe(chns)

  display(StatsBase.summarystats(chns).summaries[1])
  pars = StatsBase.summarystats(chns).summaries[1].value
  println()

  display(StatsBase.summarystats(chns, section=:internals).summaries[1])
  internals = StatsBase.summarystats(chns, section=:internals).summaries[1].value

  convert(::Type{Array}, chn::Chains) = chn.value.data
  
  chns.name_map
  
  arry = convert(Array, chns)
  size(arry)
  arry[1:5, :, 1:2]
  
  #end