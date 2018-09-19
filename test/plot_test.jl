using Test
using Plots
using StatPlots
using MCMCChain

n_iter = 500
n_name = 3
n_chain = 2

val = randn(n_iter, n_name, n_chain) .+ [1, 2, 3]'
val = hcat(val, rand(1:2, n_iter, 1, n_chain))

chn = Chains(val)

# plotting singe plotting types
ps_trace = traceplot(chn, 1)
@test isa(ps_trace, Plots.Plot)

ps_mean = meanplot(chn, 1)
@test isa(ps_mean, Plots.Plot)

ps_density = densityplot(chn, 1)
@test isa(ps_density, Plots.Plot)

ps_autocor = autocorplot(chn, 1)
@test isa(ps_autocor, Plots.Plot)

#ps_contour = plot(chn, :contour)

ps_hist = histogramplot(chn, 1)
@test isa(ps_hist, Plots.Plot)

ps_mixed = mixeddensityplot(chn, 1)
@test isa(ps_mixed, Plots.Plot)

# plotting combinations
ps_trace_mean = plot(chn)
@test isa(ps_trace_mean, Plots.Plot)

#savefig("demo-plot.png")

ps_mixed_auto = plot(chn, ptypes = [MixedDensityPlot, AutocorPlot])
@test isa(ps_mixed_auto, Plots.Plot)
