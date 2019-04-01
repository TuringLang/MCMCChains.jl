using Turing, MCMCChains, Test

#@testset "DataFrame chain summary tests" begin
  
  @model gdemo(x) = begin
      m ~ Normal(1, 0.01)
      s ~ Normal(5, 0.01)
  end

  model = gdemo([1.5, 2.0])
  sampler = HMC(1000, 0.01, 5)

  chns = [sample(model, sampler, save_state=true) for i in 1:4]
  chns = chainscat(chns...) 
  
  parm_df = summarize(chns, sections=[:parameters])

  @test 0.9 < parm_df[:m, :mean][1] < 1.1
  @test names(parm_df) == [:parameters, :mean, :std, :naive_se, :mcse, :ess]
  
  all_sections_df = summarize(chns, sections=[:parameters, :internals])
  @test all_sections_df[:parameters] == [:m, :s, :elapsed, :epsilon, :eval_num, :lf_eps, :lf_num, :lp]
  @test size(all_sections_df) == (8, 6)
  #end