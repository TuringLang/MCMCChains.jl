using Test
using StatsPlots
using MCMCChains

n_iter = 500
n_name = 3
n_chain = 2

val = randn(n_iter, n_name, n_chain) .+ [1, 2, 3]'
val = hcat(val, rand(1:2, n_iter, 1, n_chain))

chn = Chains(val)

@testset "Plotting tests" begin

    # plotting singe plotting types
    println("traceplot")
    ps_trace = traceplot(chn, 1)
    @test isa(ps_trace, Plots.Plot)

    println("meanplot")
    ps_mean = meanplot(chn, 1)
    @test isa(ps_mean, Plots.Plot)

    println("density")
    ps_density = density(chn, 1)
    @test isa(ps_density, Plots.Plot)

    ps_density = density(chn, 1, append_chains=true)
    @test isa(ps_density, Plots.Plot)

    println("autocorplot")
    ps_autocor = autocorplot(chn, 1)
    @test isa(ps_autocor, Plots.Plot)

    #ps_contour = plot(chn, :contour)

    println("histogram")
    ps_hist = histogram(chn, 1)
    @test isa(ps_hist, Plots.Plot)

    println("mixeddensity")
    ps_mixed = mixeddensity(chn, 1)
    @test isa(ps_mixed, Plots.Plot)

    # plotting combinations
    ps_trace_mean = plot(chn)
    @test isa(ps_trace_mean, Plots.Plot)

    ps_trace_mean = plot(chn, append_chains=true)
    @test isa(ps_trace_mean, Plots.Plot)

    savefig("demo-plot.png")

    ps_mixed_auto = plot(chn, seriestype = (:mixeddensity, :autocorplot))
    @test isa(ps_mixed_auto, Plots.Plot)

    # Test plotting using colordim keyword
    p_colordim = plot(chn, colordim = :parameter)
    @test isa(p_colordim, Plots.Plot)

    # Test if plotting a sub-set work.s
    p_subset = plot(chn, 2)
    @test isa(p_subset, Plots.Plot)

    p_subset_colordim = plot(chn, 2, colordim = :parameter)
    @test isa(p_subset_colordim, Plots.Plot)

    rm("demo-plot.png")
end
