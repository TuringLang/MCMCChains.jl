using Test
using StatsPlots
using MCMCChains

import Logging

unicodeplots()

n_iter = 500
n_name = 3
n_chain = 3

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
    display(density(chn, 1, append_chains = true))
    println()

    println("autocorplot")
    display(autocorplot(chn, 1))
    println()

    println("ridgelineplot")
    display(ridgelineplot(chn, chn.name_map[:parameters]))
    println()

    println("forestplot")
    display(forestplot(chn, chn.name_map[:parameters]))
    println()

    #ps_contour = plot(chn, :contour)

    println("histogram")
    display(histogram(chn, 1))
    println()

    println("\nmixeddensity")
    display(mixeddensity(chn, 1))

    println("corner")
    display(corner(chn[:, 1:2, :]))
    # https://github.com/TuringLang/MCMCChains.jl/issues/413
    display(corner(chn[:, 1:2, 2:3]))

    # Violinplot tests
    println("\nviolinplot")
    display(violinplot(chn)) # All parameters, default colordim (:chain)
    display(violinplot(chn, colordim = :parameter)) # All chains, colordim = :parameter
    display(violinplot(chn, 1)) # Single parameter, default colordim (:chain)
    display(violinplot(chn, 1, colordim = :parameter)) # Single chain, colordim = :parameter
    display(violinplot(chn, 1, show_boxplot = false)) # Single parameter, no boxplot
    display(violinplot(chn, 1, append_chains = true)) # Single parameter, chains appended
    display(violinplot(chn, append_chains = true)) # All parameters, chains appended
    println()

    # Plot() with violinplot seriestype
    println("\nplot() with violinplot seriestype")
    # "seriestype = :violin will" also work fine
    display(plot(chn, seriestype = :violinplot)) # All parameters with violinplot
    display(plot(chn, 1, seriestype = :violinplot)) # Specific parameter(s) with violinplot
    display(plot(chn, 1, seriestype = :violinplot, show_boxplot = false)) # Specific parameter(s) with violinplot and no boxplot
    display(plot(chn, seriestype = :violinplot, append_chains = true)) # All parameters, chains appended
    println()

    # plotting combinations
    display(plot(chn))
    display(plot(chn, append_chains = true))
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
