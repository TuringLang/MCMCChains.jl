using MCMCChains

using Random
using StatsPlots
using Test

@testset "PPC Plot Tests" begin
    Random.seed!(123)
    n_iter = 1000
    n_params = 2
    n_chains = 2
    
    posterior_data = randn(n_iter, n_params, n_chains)
    posterior_chains = Chains(posterior_data, [:μ, :σ])
    
    pp_data = zeros(n_iter, 10, n_chains)
    for i in 1:n_iter, j in 1:n_chains
        μ = posterior_data[i, 1, j]
        σ = abs(posterior_data[i, 2, j]) + 0.1
        pp_data[i, :, j] = randn(10) * σ .+ μ
    end
    pp_chains = Chains(pp_data)
    
    Random.seed!(456)
    observed = randn(10) * 1.2 .+ 0.5
    
    @testset "Basic PPC Plot Creation" begin
        # Test that plots can be created without errors
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed)
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; kind=:density)
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; kind=:histogram)
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; kind=:scatter)
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; kind=:cumulative)
    end
    
    @testset "PPC Plot Arguments Validation" begin
        # Test error handling
        @test_throws ErrorException ppcplot()  # No arguments
        @test_throws ErrorException ppcplot(posterior_chains)  # Too few arguments
        @test_throws ErrorException ppcplot(posterior_chains, pp_chains)  # Too few arguments
        @test_throws ErrorException ppcplot("not_chains", pp_chains, observed)  # Wrong type
        @test_throws ErrorException ppcplot(posterior_chains, "not_chains", observed)  # Wrong type
        @test_throws ErrorException ppcplot(posterior_chains, pp_chains, "not_vector")  # Wrong type
        @test_throws ErrorException ppcplot(posterior_chains, pp_chains, observed; kind=:invalid)
        @test_throws ErrorException ppcplot(posterior_chains, pp_chains, observed; colors=["red"])  # Wrong length
    end
    
    @testset "PPC Plot Options" begin
        # Test various options work without error
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; alpha=0.5)
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; num_pp_samples=10)
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; mean_pp=false)
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; observed_rug=true)
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; colors=["blue", "red", "green"])
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; jitter=0.3)
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; legend=false)
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; random_seed=42)
        
        # Test scatter-specific defaults
        @test_nowarn ppcplot(posterior_chains, pp_chains, observed; kind=:scatter, num_pp_samples=3)
    end
    
    @testset "Edge cases" begin
        # Test with single chain
        single_chain_post = Chains(randn(100, 2, 1), [:μ, :σ])
        single_chain_pp = Chains(randn(100, 5, 1))
        obs_single = randn(5)
        @test_nowarn ppcplot(single_chain_post, single_chain_pp, obs_single)
        
        # Test with small sample size
        small_post = Chains(randn(10, 2, 1), [:μ, :σ])
        small_pp = Chains(randn(10, 3, 1))
        obs_small = randn(3)
        @test_nowarn ppcplot(small_post, small_pp, obs_small)
        
        @test_logs (:warn, r"Requested .* samples but only .* available") ppcplot(small_post, small_pp, obs_small; num_pp_samples=100)  # Should warn and use all
    end
    
    @testset "Group Types and Observed Data Control" begin
        # Test explicit group types
        mcmc_chains = Chains(randn(100, 2, 2), [:μ, :σ])  # Multiple chains, MCMC names
        pp_chains_test = Chains(randn(100, 5, 2))
        obs_test = randn(5)
        @test_nowarn ppcplot(mcmc_chains, pp_chains_test, obs_test; ppc_group=:posterior)
        
        # Test prior group (should hide observed data by default)
        prior_chains = Chains(randn(50, 2, 1), [:param_1, :param_2])  # Single chain, generic names
        prior_pp_chains = Chains(randn(50, 5, 1))
        @test_nowarn ppcplot(prior_chains, prior_pp_chains, obs_test; ppc_group=:prior)
        
        # Test explicit observed control
        @test_nowarn ppcplot(mcmc_chains, pp_chains_test, obs_test; ppc_group=:prior, observed=true)
        @test_nowarn ppcplot(mcmc_chains, pp_chains_test, obs_test; ppc_group=:posterior, observed=false)
        
        # Test invalid group type
        @test_throws ErrorException ppcplot(mcmc_chains, pp_chains_test, obs_test; ppc_group=:invalid)
    end
end
