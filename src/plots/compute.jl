function _compute_plot_data(
    i::Integer,
    chains::Chains,
    par_names::AbstractVector{Symbol};
    hpd_val = [0.05, 0.2],
    q = [0.1, 0.9],
    spacer = 0.4,
    _riser = 0.2,
    barbounds = (-Inf, Inf),
    show_mean = true,
    show_median = true,
    show_qi = false,
    show_hpdi = true,
    fill_q = true,
    fill_hpd = false,
    ordered = false,
)

    chain_dic = Dict(zip(quantile(chains)[:, 1], quantile(chains)[:, 4]))
    sorted_chain = sort(collect(zip(values(chain_dic), keys(chain_dic))))
    sorted_par = [sorted_chain[i][2] for i = 1:length(par_names)]
    par = (ordered ? sorted_par : par_names)
    hpdi = sort(hpd_val)

    chain_sections = MCMCChains.group(chains, Symbol(par[i]))
    chain_vec = vec(chain_sections.value.data)
    lower_hpd =
        [MCMCChains.hpd(chain_sections, alpha = hpdi[j]).nt.lower for j = 1:length(hpdi)]
    upper_hpd =
        [MCMCChains.hpd(chain_sections, alpha = hpdi[j]).nt.upper for j = 1:length(hpdi)]
    h = _riser + spacer * (i - 1)
    qs = quantile(chain_vec, q)
    k_density = kde(chain_vec)
    if fill_hpd
        x_int = filter(x -> lower_hpd[1][1] <= x <= upper_hpd[1][1], k_density.x)
        val = pdf(k_density, x_int) .+ h
    elseif fill_q
        x_int = filter(x -> qs[1] <= x <= qs[2], k_density.x)
        val = pdf(k_density, x_int) .+ h
    else
        x_int = k_density.x
        val = k_density.density .+ h
    end
    chain_med = median(chain_vec)
    chain_mean = mean(chain_vec)
    min = minimum(k_density.density .+ h)
    q_int = (show_qi ? [qs[1], chain_med, qs[2]] : [chain_med])

    return (;
        par,
        hpdi,
        lower_hpd,
        upper_hpd,
        h,
        qs,
        k_density,
        x_int,
        val,
        chain_med,
        chain_mean,
        min,
        q_int,
    )
end

"""
    _interval_rows(chains, par_names; kwargs...)

Per-parameter plot data for the ridgeline and forest recipes.

Computing every row up front lets the recipes size the axis to the tallest ridge, which a
per-row loop cannot do because the first row does not know about the others.
"""
function _interval_rows(chains::Chains, par_names::AbstractVector{Symbol}; kwargs...)
    return [_compute_plot_data(i, chains, par_names; kwargs...) for i = 1:length(par_names)]
end

"""
    _rank_bin_counts(val, nbins)

Rank every draw against every other draw, then count each chain's ranks into `nbins`
equal-width bins.

`val` is iterations by chains. Ranking pools all chains, so under good mixing each chain
holds an equal share of every rank range and its counts are flat. Returns the bin edges and
a bins by chains matrix of counts.
"""
function _rank_bin_counts(val::AbstractMatrix, nbins::Integer)
    niter, nchains = size(val)
    ntotal = niter * nchains
    ranks = reshape(ordinalrank(vec(val)), niter, nchains)
    edges = range(0.5, ntotal + 0.5; length = nbins + 1)
    counts = zeros(Int, nbins, nchains)
    width = (ntotal) / nbins
    for c = 1:nchains, r in view(ranks, :, c)
        b = min(nbins, Int(fld(r - 1, width)) + 1)
        counts[b, c] += 1
    end
    return edges, counts
end

"""
    _bin_count(values, rule)

How many bins to cut `values` into.

`:sturges` reads the count off the sample size alone, which is the oldest rule and gives too
few bins on a large sample. `:fd` is Freedman and Diaconis (1981), width `2 IQR / cbrt(n)`,
rank based so that one far-out draw does not widen every bin. `:auto` takes whichever of the
two asks for more, so that a sample with no spread between its quartiles still gets a usable
plot.

# References
Freedman and Diaconis (1981). On the histogram as a density estimator: L2 theory. Zeitschrift für Wahrscheinlichkeitstheorie und verwandte Gebiete 57(4). https://doi.org/10.1007/BF01025868
"""
function _bin_count(values, rule::Symbol)
    rule in (:auto, :fd, :sturges) ||
        throw(ArgumentError("`bins` must be `:auto`, `:fd`, `:sturges`, a count or edges"))

    n = length(values)
    n > 1 || return 1
    lo, hi = extrema(values)
    hi > lo || return 1
    span = hi - lo

    sturges = ceil(Int, log2(n)) + 1
    rule === :sturges && return sturges

    iqr = quantile(values, 0.75) - quantile(values, 0.25)
    fd = iqr > 0 ? ceil(Int, span * cbrt(n) / (2 * iqr)) : 0
    rule === :fd && return clamp(fd, 1, 1000)
    return clamp(max(fd, sturges), 1, 1000)
end

"""
    _histogram_edges(series, bins)

Bin edges covering every chain of one parameter.

Binning each chain over its own range gives bars of different widths, which cannot be
compared, so the edges come from the chains pooled. `bins` is a rule, a count, or edges to
use as they are. The extra bin past the top is there because the last edge is open, and the
largest draw would otherwise be counted nowhere.
"""
function _histogram_edges(series, bins)
    pooled = collect(Iterators.flatten(series))
    bins isa Symbol && (bins = _bin_count(pooled, bins))
    bins isa Integer || return bins

    lo, hi = extrema(pooled)
    hi > lo || return range(lo - 0.5, lo + 0.5; length = bins + 1)
    width = (hi - lo) / bins
    return range(lo, hi + width; length = bins + 2)
end

"""
    _ecdf_series(val)

Sorted draws and the empirical cumulative probability at each, one column per chain.

Chains that agree lie on top of each other, which is easier to judge than two density
curves because the eye compares positions rather than areas.
"""
function _ecdf_series(val::AbstractMatrix)
    nchains = size(val, 2)
    xs = [sort!(collect(skipmissing(val[:, k]))) for k = 1:nchains]
    ys = [collect(1:length(x)) ./ length(x) for x in xs]
    return xs, ys
end

"""
    _evolution_points(ndraws, npoints)

Draw counts at which to recompute a diagnostic, spread evenly and always ending at `ndraws`.
"""
function _evolution_points(ndraws::Integer, npoints::Integer)
    npoints = max(2, npoints)
    step = max(1, ndraws ÷ npoints)
    points = collect(step:step:ndraws)
    isempty(points) && return [ndraws]
    last(points) == ndraws || push!(points, ndraws)
    return points
end

"""
    _evolution(chains, diagnostic, npoints)

A diagnostic recomputed on growing prefixes of the chain.

A single number says how well the sampler did in the end. Watching it against draw count
says whether it is still improving, which is what tells you to run for longer.
"""
function _evolution(chains::Chains, diagnostic::Symbol, npoints::Integer)
    sub = Chains(chains, :parameters)
    par_names = names(sub)
    ndraws = size(sub, 1)
    points = _evolution_points(ndraws, npoints)
    values = Matrix{Float64}(undef, length(points), length(par_names))

    for (i, n) in enumerate(points)
        prefix = sub[1:n, :, :]
        column = if diagnostic === :ess
            MCMCDiagnosticTools.ess(prefix).nt.ess
        elseif diagnostic === :rhat
            MCMCDiagnosticTools.rhat(prefix).nt.rhat
        else
            throw(ArgumentError("`diagnostic` must be `:ess` or `:rhat`, got `$diagnostic`"))
        end
        values[i, :] = collect(column)
    end
    return points, par_names, values
end

# The diagnostics worth seeing side by side, how to read each, and what counts as good.
# The effective sample size and R-hat thresholds are the ones recommended in Vehtari et al.
# (2021), https://doi.org/10.1214/20-BA1221.
const _DIAGNOSTIC_COLUMNS = (
    (name = "R-hat", good = v -> v <= 1.01, warn = v -> v <= 1.05, digits = 3),
    (name = "Bulk ESS", good = v -> v >= 400, warn = v -> v >= 100, digits = 0),
    (name = "Tail ESS", good = v -> v >= 400, warn = v -> v >= 100, digits = 0),
    (name = "ESS / draw", good = v -> v >= 0.5, warn = v -> v >= 0.1, digits = 2),
    (name = "MCSE / SD", good = v -> v <= 0.05, warn = v -> v <= 0.1, digits = 3),
)

"""
    _diagnostics_table(chains)

Parameter names, column names, the value in each cell, and how each cell scores.

A score is 2 for good, 1 for worth a look, 0 for bad and -1 where the diagnostic does not
apply, such as a tail effective sample size for a parameter that takes two values. Reading
one table beats reading five plots, because the question is almost always whether anything
at all is wrong.
"""
function _diagnostics_table(chains::Chains)
    sub = Chains(chains, :parameters)
    par_names = names(sub)
    ndraws, _, nchains = size(sub)
    total = ndraws * nchains

    bulk = collect(MCMCDiagnosticTools.ess(sub; kind = :bulk).nt.ess)
    tail = collect(MCMCDiagnosticTools.ess(sub; kind = :tail).nt.ess)
    rhats = collect(MCMCDiagnosticTools.rhat(sub).nt.rhat)
    mcses = collect(MCMCDiagnosticTools.mcse(sub).nt.mcse)
    _, draws = _draws_by_parameter(sub)
    spread = vec(std(draws; dims = 1))

    values = hcat(rhats, bulk, tail, bulk ./ total, mcses ./ spread)
    scores = similar(values)
    for (j, column) in enumerate(_DIAGNOSTIC_COLUMNS)
        for i in Base.axes(values, 1)
            v = values[i, j]
            scores[i, j] = if !isfinite(v)
                -1.0
            elseif column.good(v)
                2.0
            elseif column.warn(v)
                1.0
            else
                0.0
            end
        end
    end
    return par_names, [c.name for c in _DIAGNOSTIC_COLUMNS], values, scores
end

"""
    _format_cell(value, digits)

A diagnostic rounded for reading rather than for precision.
"""
function _format_cell(value, digits::Integer)
    isfinite(value) || return "n/a"
    digits == 0 && return string(round(Int, value))
    return string(round(value; digits = digits))
end

"""
    _draws_by_parameter(chains)

Every draw of every parameter as a draws by parameters matrix, chains stacked.

Row order is iteration fastest, so it lines up with `vec` of any iterations by chains
statistic taken from the same chain.
"""
function _draws_by_parameter(chains::Chains)
    sub = Chains(chains, :parameters)
    vals = Array(sub.value.data)
    niter, npar, nchains = size(vals)
    return names(sub), reshape(permutedims(vals, (1, 3, 2)), niter * nchains, npar)
end

"""
    _internal_names(chains)

The sampler statistics a chain carries, empty when it has no `:internals` section.
"""
function _internal_names(chains::Chains)
    haskey(chains.name_map, :internals) || return Symbol[]
    return names(chains, :internals)
end

"""
    _divergences(chains)

A draw-length `Bool` vector marking divergent transitions, or `nothing` when the chain does
not record them.
"""
function _divergences(chains::Chains)
    internals = _internal_names(chains)
    for key in (:numerical_error, :divergent)
        if key in internals
            return vec(Array(chains[:, key, :])) .> 0
        end
    end
    return nothing
end

"""
    _standardise(m)

Centre and scale every column of `m`, leaving constant columns at zero.

Parallel coordinates put every parameter on one axis, so parameters measured on different
scales have to be made comparable first.
"""
function _standardise(m::AbstractMatrix)
    out = float(copy(m))
    for j in Base.axes(out, 2)
        col = view(out, :, j)
        centre = mean(col)
        spread = std(col)
        col .= spread > 0 ? (col .- centre) ./ spread : zero(eltype(out))
    end
    return out
end

"""
    _thin_draws(ndraws, keep, rng)

Row indices of at most `keep` draws, chosen uniformly.

Thinning is uniform over all draws, divergent ones included, so that the share of the lines
that are divergent stays the share of the draws that are divergent.
"""
function _thin_draws(ndraws::Integer, keep, rng)
    (keep === nothing || ndraws <= keep) && return collect(1:ndraws)
    return sort!(Random.randperm(rng, ndraws)[1:keep])
end

# Statistics an HMC sampler records, and what to call them on an axis.
const _NUTS_STATISTICS = (
    acceptance = (:acceptance_rate, "Acceptance rate"),
    stepsize = (:step_size, "Step size"),
    treedepth = (:tree_depth, "Tree depth"),
    divergence = (:numerical_error, "Divergent transitions"),
)

"""
    _nuts_statistic(chains, kind)

An iterations by chains matrix of the requested sampler statistic, and its axis label.
"""
function _nuts_statistic(chains::Chains, kind::Symbol)
    haskey(_NUTS_STATISTICS, kind) || throw(
        ArgumentError(
            "`kind` must be one of $(join(keys(_NUTS_STATISTICS), ", ")), got `$kind`",
        ),
    )
    key, label = _NUTS_STATISTICS[kind]
    if key ∉ _internal_names(chains)
        throw(
            ArgumentError(
                "`$key` is not in the chain's internals, so `kind = :$kind` is not available. " *
                "Sampler statistics like this one come from HMC samplers such as NUTS.",
            ),
        )
    end
    return Array(chains[:, key, :]), label
end
