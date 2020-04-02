using Test, Random
using MCMCChains
using Distributions

## Test Chain
# Define the experiment
n_iter = 4000
n_name = 3
n_chain = 2

# observations
Random.seed!(1234)
x = map(i -> i + 2e-1*randn(), 1:n_name)

# some sample experiment results
Random.seed!(1234)
val1 = 0.5*randn(n_iter, n_name, n_chain) .+ x'
val2 = 2*randn(n_iter, n_name, n_chain) .+ x'

# construct a Chains object
chn1 = Chains(val1)
chn2 = Chains(val2)

lpfun = function f(chain::Chains)

    p1 = Array(chain[:,1,:])
    p2 = Array(chain[:,2,:])
    p3 = Array(chain[:,3,:])

    lp = map(p -> logpdf(Normal(p), x[1]), p1)
    lp += map(p -> logpdf(Normal(p), x[2]), p2)
    lp += map(p -> logpdf(Normal(p), x[3]), p3)

    return lp
end

@testset "deviance information criterion" begin

    DIC1, pD1 = dic(chn1, lpfun)
    DIC2, pD2 = dic(chn2, lpfun)

    @test DIC1 < DIC2
    @test pD1 < pD2
end

@testset "hpd tests" begin
    result = hpd(chn1)
    @test all(result.nt.upper .> result.nt.lower)
end
