println("imports")
@time using Test
@time using StatsPlots
@time using MCMCChains

n_iter = 500
n_name = 3
n_chain = 2

println("setup")
@time val = randn(n_iter, n_name, n_chain) .+ [1, 2, 3]'
@time val = hcat(val, rand(1:2, n_iter, 1, n_chain))

@time chn = Chains(val)

@testset "Plotting tests" begin

    # plotting singe plotting types
    println("trace")
    @time ps_trace = traceplot(chn, 1)
    @test isa(ps_trace, Plots.Plot)

    println("mean")
    @time ps_mean = meanplot(chn, 1)
    @test isa(ps_mean, Plots.Plot)

    println("density")
    @time ps_density = density(chn, 1)
    @test isa(ps_density, Plots.Plot)

    @time ps_density = density(chn, 1, append_chains=true)
    @test isa(ps_density, Plots.Plot)

    println("autocor")
    @time ps_autocor = autocorplot(chn, 1)
    @test isa(ps_autocor, Plots.Plot)

    #ps_contour = plot(chn, :contour)

    println("histogram")
    @time ps_hist = histogram(chn, 1)
    @test isa(ps_hist, Plots.Plot)

    println("mixeddensity")
    @time ps_mixed = mixeddensity(chn, 1)
    @test isa(ps_mixed, Plots.Plot)

    # plotting combinations
    println("trace_mean")
    @time ps_trace_mean = plot(chn)
    @test isa(ps_trace_mean, Plots.Plot)

    @time ps_trace_mean = plot(chn, append_chains=true)
    @test isa(ps_trace_mean, Plots.Plot)

    savefig("demo-plot.png")

    println("mixedauto")
    @time ps_mixed_auto = plot(chn, seriestype = (:mixeddensity, :autocorplot))
    @test isa(ps_mixed_auto, Plots.Plot)

    # Test plotting using colordim keyword
    println("colordim")
    @time p_colordim = plot(chn, colordim = :parameter)
    @test isa(p_colordim, Plots.Plot)

    # Test if plotting a sub-set work.s
    println("subset")
    @time p_subset = plot(chn, 2)
    @test isa(p_subset, Plots.Plot)

    println("subsetcolordim")
    @time p_subset_colordim = plot(chn, 2, colordim = :parameter)
    @test isa(p_subset_colordim, Plots.Plot)

    rm("demo-plot.png")
end
