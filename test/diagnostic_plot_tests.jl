using Test
using MCMCChains
using StatsPlots
using Random

@testset "Diagnostic plots" begin
    Random.seed!(42)
    niter, nchains = 400, 4

    chn = Chains(randn(niter, 3, nchains), [:a, :b, :c])

    nuts_val = hcat(
        randn(niter, 2, nchains),
        clamp.(0.9 .+ 0.05 .* randn(niter, 1, nchains), 0, 1),
        0.3 .+ 0.01 .* randn(niter, 1, nchains),
        float.(rand(2:5, niter, 1, nchains)),
        float.(rand(niter, 1, nchains) .< 0.05),
    )
    nuts_names = [:x, :y, :acceptance_rate, :step_size, :tree_depth, :numerical_error]
    chn_nuts = Chains(
        nuts_val,
        nuts_names,
        (
            parameters = [:x, :y],
            internals = [:acceptance_rate, :step_size, :tree_depth, :numerical_error],
        ),
    )

    @testset "per-parameter diagnostics" begin
        @test essplot(chn) isa Plots.Plot
        @test essplot(chn, kind = :tail) isa Plots.Plot
        @test essplot(chn, relative = true) isa Plots.Plot
        @test rhatplot(chn) isa Plots.Plot
        @test mcseplot(chn) isa Plots.Plot
        @test mcseplot(chn, relative = true) isa Plots.Plot

        # One dot per parameter, plus the dashed reference line.
        p = rhatplot(chn)
        @test length(p.series_list) == 2
        @test length(p.series_list[2][:x]) == 3

        # No reference line to draw for a Monte Carlo standard error.
        @test length(mcseplot(chn).series_list) == 1
    end

    @testset "diagnostics are the ones reported" begin
        p = rhatplot(chn)
        @test p.series_list[2][:x] ≈ collect(rhat(chn).nt.rhat)

        p = essplot(chn, kind = :tail)
        @test p.series_list[2][:x] ≈ collect(ess(chn, kind = :tail).nt.ess)

        # The relative version is the same numbers over the draw count.
        p = essplot(chn, relative = true)
        @test p.series_list[2][:x] ≈ collect(ess(chn).nt.ess) ./ (niter * nchains)
    end

    @testset "only a chain is accepted" begin
        @test_throws ArgumentError rhatplot(randn(10, 2))
        @test_throws ArgumentError essplot(chn, chn)
    end

    @testset "the whole diagnostics table at once" begin
        @test diagnosticsplot(chn) isa Plots.Plot

        par_names, columns, values, scores = MCMCChains._diagnostics_table(chn)
        @test par_names == names(chn)
        @test length(columns) == 5
        @test size(values) == (3, 5)

        # The first column is R-hat, and it is the R-hat that is reported.
        @test values[:, 1] ≈ collect(rhat(chn).nt.rhat)
        # Chains from one distribution score well on everything.
        @test all(scores .== 2)

        # A diagnostic that cannot be computed is marked apart from a bad one.
        flat = Chains(fill(1.0, niter, 1, nchains), [:flat])
        _, _, flat_values, flat_scores = MCMCChains._diagnostics_table(flat)
        @test any(!isfinite, flat_values)
        @test any(flat_scores .== -1)

        # A chain that has not mixed scores badly on R-hat.
        split =
            Chains(cat(randn(niter, 1, 2), randn(niter, 1, 2) .+ 10; dims = 3), [:split])
        _, _, _, split_scores = MCMCChains._diagnostics_table(split)
        @test split_scores[1, 1] == 0
    end

    @testset "diagnostics against draws" begin
        @test evolutionplot(chn) isa Plots.Plot
        @test evolutionplot(chn, diagnostic = :rhat) isa Plots.Plot
        @test_throws ArgumentError evolutionplot(chn, diagnostic = :nonsense)

        points, par_names, values = MCMCChains._evolution(chn, :ess, 5)
        @test last(points) == niter          # always ends at the whole chain
        @test issorted(points)
        @test size(values) == (length(points), length(par_names))

        # The last point is the diagnostic for the whole chain.
        @test values[end, :] ≈ collect(ess(chn).nt.ess)
        # More draws, more information.
        @test values[end, 1] > values[1, 1]
    end

    @testset "empirical CDF" begin
        p = plot(chn, seriestype = :ecdfplot)
        @test p isa Plots.Plot

        xs, ys = MCMCChains._ecdf_series(randn(100, 3))
        @test length(xs) == 3
        @test all(issorted, xs)
        @test all(y -> last(y) ≈ 1, ys)
        @test all(y -> issorted(y), ys)
    end

    @testset "parallel coordinates" begin
        @test parallelplot(chn) isa Plots.Plot
        @test parallelplot(chn, standardise = false) isa Plots.Plot

        # Without divergences there is one bundle, with them there are two.
        @test length(parallelplot(chn).series_list) == 1
        @test length(parallelplot(chn_nuts).series_list) == 2

        # All the draws go into one series, broken between lines.
        p = parallelplot(chn)
        @test length(p.series_list[1][:y]) == niter * nchains * 4

        # Thinning keeps the requested number of lines and is reproducible.
        p = parallelplot(chn, num_draws = 50, random_seed = 7)
        q = parallelplot(chn, num_draws = 50, random_seed = 7)
        @test length(p.series_list[1][:y]) == 50 * 4
        @test isequal(p.series_list[1][:y], q.series_list[1][:y])

        @test_throws ArgumentError parallelplot(Chains(randn(niter, 1, nchains), [:a]))
    end

    @testset "standardising" begin
        m = [1.0 100.0; 3.0 300.0; 5.0 500.0]
        s = MCMCChains._standardise(m)
        @test all(abs.(sum(s; dims = 1)) .< 1e-10)
        @test s[:, 1] ≈ s[:, 2]
        # A parameter that never moved has nowhere to be placed but the middle.
        @test all(MCMCChains._standardise(fill(2.0, 4, 1)) .== 0)
    end

    @testset "NUTS statistics" begin
        for kind in (:acceptance, :stepsize, :treedepth, :divergence)
            @test nutsplot(chn_nuts, kind = kind) isa Plots.Plot
        end

        # The divergence count is per chain, and is the number recorded.
        p = nutsplot(chn_nuts, kind = :divergence)
        drawn = filter(!isnan, p.series_list[1][:y])
        expected = vec(sum(Array(chn_nuts[:, :numerical_error, :]) .> 0; dims = 1))
        @test sort(unique(filter(>(0), drawn))) == sort(unique(expected))

        @test_throws ArgumentError nutsplot(chn_nuts, kind = :nonsense)
        # A chain without sampler statistics says so rather than failing deep inside.
        @test_throws ArgumentError nutsplot(chn)
    end
end
