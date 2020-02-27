using Test
using StatsPlots
using MCMCChains

import Logging

unicodeplots()

n_iter = 500
n_name = 3
n_chain = 2

val = randn(n_iter, n_name, n_chain) .+ [1, 2, 3]'
val = hcat(val, rand(1:2, n_iter, 1, n_chain))

chn = Chains(val)

# Silence all warnings.
level = Logging.min_enabled_level(Logging.current_logger())
Logging.disable_logging(Logging.Warn)

@testset "Plotting tests" begin
    # plotting singe plotting types
    println("traceplot")
    display(traceplot(chn, 1))
    println()
    
    println("meanplot")
    display(meanplot(chn, 1))
    println()
    
    println("density")
    display(density(chn, 1))
    display(density(chn, 1, append_chains=true))
    println()
    
    println("autocorplot")
    display(autocorplot(chn, 1))
    println()
    
    #ps_contour = plot(chn, :contour)

    println("histogram")
    display(histogram(chn, 1))
    println()
    
    println("\nmixeddensity")
    display(mixeddensity(chn, 1))
    
    # plotting combinations
    display(plot(chn))
    display(plot(chn, append_chains=true))
    display(plot(chn, seriestype = (:mixeddensity, :autocorplot)))

    # Test plotting using colordim keyword
    display(plot(chn, colordim = :parameter))
    
    # Test if plotting a sub-set work.s
    display(plot(chn, 2))
    display(plot(chn, 2, colordim = :parameter))
    println()
end

# Reset log level.
Logging.disable_logging(level)