using MCMCChains

function simulate_ar1(ρ, N; σ = 1.0)
    x = Vector{Float64}(undef, N)
    z = σ * randn() / √(1 - ρ^2)
    for i in 1:N
        z = ρ*z + randn()*σ
        x[i] = z
    end
    x
end

ρ = 0.5
N = 100
nparams = 3
ar1_ess_factor(ρ) = 1/(1 + 2*ρ/(1-ρ))

vals = reshape(hcat([simulate_ar1(ρ, N) for _ in 1:nparams]...), (N, nparams, 1))
chn = Chains(vals)
e = ess(chn)[2] ./ (ar1_ess_factor(ρ)*N)
display(ess(chn))


# actual = ar1_ess_factor(ρ)*N
# ######################################
# chn2 = Chains(randn(1000,10,1));
# ess(chn2)
# # @test 9900 ≤ effective_sample_size(vcat(chains...)) ≤ 10100
# # @test 1 ≤ potential_scale_reduction(chains...) ≤ 1.001
#
# # Tamas:
# # │ Row │ parameters │ getfield(MCMCChains, Symbol("#ess_f#516"))() │
# # │     │ Symbol     │ Float64                                      │
# # ├─────┼────────────┼──────────────────────────────────────────────┤
# # │ 1   │ Param1     │ 25.6979                                      │
# # │ 2   │ Param2     │ 38.2962                                      │
# # │ 3   │ Param3     │ 53.374                                       │
#
# # Cameron:
# # │ Row │ parameters │ ess     │ r_hat    │
# # │     │ Symbol     │ Any     │ Any      │
# # ├────┼──────────┼───────┼─────────┤
# # │ 1   │ Param1     │ 9.80762 │ 1.04235  │
# # │ 2   │ Param2     │ 15.1556 │ 0.997153 │
# # │ 3   │ Param3     │ 19.7391 │ 1.06333  │
#
# MCMCChains.autocorrelation([1,2,3], 1)
