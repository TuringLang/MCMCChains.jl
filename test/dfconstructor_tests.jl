using Turing, MCMCChains, Test

@testset "DataFrame constructor tests" begin
  
  @model gdemo(x) = begin
      m ~ Normal(1, 0.01)
      s ~ Normal(5, 0.01)
  end

  model = gdemo([1.5, 2.0])
  sampler = HMC(1000, 0.01, 5)

  chns = [sample(model, sampler, save_state=true) for i in 1:4]
  chn = chainscat(chns...)
  
  df = DataFrame(chn)
  @test size(df) == (4000, 8)
  
  df1 = DataFrame(chn, [:parameters])
  @test size(df1) == (4000, 2)
  
  df2 = DataFrame(chn, [:internals, :parameters])
  @test size(df2) == (4000, 8)
  
  df3 = DataFrame(chn[:s]) 
  @test size(df3) == (4000, 1)
  
  df4 = DataFrame(chn[:lp])  
  @test size(df4) == (4000, 1)
  
  df5 = DataFrame(chn, [:parameters], append_chains=false)  
  @test size(df5) == (4, )
  @test size(df5[1]) == (1000, 2)

  df6 = DataFrame(chn, [:parameters, :internals], append_chains=false)  
  @test size(df6) == (4, )
  @test size(df6[1]) == (1000, 8)

  df7 = DataFrame(chn, [:parameters, :internals], remove_missing_union=false)  
  @test size(df7) == (4000, 8)

  df8 = DataFrame(chn, [:parameters, :internals], remove_missing_union=false,
    append_chains=false)  
  @test size(df8) == (4, )
  @test size(df8[1]) == (1000, 8)

end