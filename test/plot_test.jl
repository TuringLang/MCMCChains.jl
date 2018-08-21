using Chain

n_iter = 500
n_name = 3
n_chain = 2

val = randn(n_iter, n_name, n_chain) .+ [1, 2, 3]'
val = hcat(val, rand(1:2, n_iter, 1, n_chain))

chn = Chains(val)

# plotting singe plotting types
ps_trace = plot(chn, :trace)
ps_mean = plot(chn, :mean)
ps_density = plot(chn, :density)
ps_autocor = plot(chn, :autocor)
#ps_contour = plot(chn, :contour)
ps_hist = plot(chn, :histogram)
ps_mixed = plot(chn, :mixeddensity)

# plotting combinations
ps_trace_mean = plot(chn, [:trace, :mean])
ps_mixed_auto = plot(chn, [:mixeddensity, :autocor])
