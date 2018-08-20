# run this file in interactive mode to pop-up plots (`julia -i plot_test.jl`)
# NOTE: I don't understand how Plots evals plotting call - I
#       tried to add `unicodeplots()` to plot in terminal but doesn't work.
#       However, the same code works in interactive mode.
using Chain

n_iter = 500
n_name = 3
n_chain = 1

val = randn(n_iter, n_name, n_chain)
chn = Chains(val)

ps_trace = plot(chn, :trace)
ps_mean = plot(chn, :mean)
ps_density = plot(chn, :density)
ps_autocor = plot(chn, :autocor)
ps_contour = plot(chn, :contour)
ps_bar = plot(chn, :bar)
