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

    # Only worked out when it is asked for. Sorting runs `quantile` over every parameter
    # twice, once per row of the plot, and on a chain holding one parameter the table it
    # reads back is a scalar rather than a column.
    par = if ordered
        summary = quantile(chains)
        medians = Dict(zip(summary[:, 1], summary[:, 4]))
        sorted = sort(collect(zip(values(medians), keys(medians))))
        [sorted[j][2] for j = 1:length(par_names)]
    else
        par_names
    end
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
    niter, _ = size(val)
    flat = vec(val)
    keep = .!ismissing.(flat)
    ntotal = count(keep)
    # Tied draws take the average of the ranks they span. Breaking ties by position would
    # order them by which chain they came from, which draws a well mixed discrete parameter
    # as though every chain held its own stretch of the range.
    ranks = fill(NaN, length(flat))
    ranks[keep] = tiedrank(collect(skipmissing(flat)))

    edges = range(0.5, ntotal + 0.5; length = nbins + 1)
    counts = zeros(Int, nbins, size(val, 2))
    width = ntotal / nbins
    for (k, r) in enumerate(ranks)
        isnan(r) && continue
        b = clamp(Int(fld(r - 1, width)) + 1, 1, nbins)
        counts[b, cld(k, niter)] += 1
    end
    return edges, counts
end

# The binning rules this computes itself, so that every chain of a parameter shares one set
# of edges. Plots names these and `:wand`, which is left to Plots.
const BIN_RULES = (:auto, :fd, :sturges, :sqrt, :rice, :scott)

const MAX_DISCRETE_LEVELS = 20

"""
    _discrete_edges(values)

Return bins for a low-cardinality parameter, grouping adjacent values when needed, or `nothing`.
"""
function _discrete_edges(values)
    levels = _few_levels(values)
    levels === nothing && return nothing
    lo, hi = first(levels), last(levels)

    runs = _runs(levels)
    length(runs) == 1 && return _one_bin(lo, hi)
    return _bins_around(levels, runs, lo, hi)
end

"""
    _runs(levels)

Group consecutive representable values into runs that share one bar.
"""
function _runs(levels)
    runs = UnitRange{Int}[]
    start = firstindex(levels)
    for i = (start+1):lastindex(levels)
        nextfloat(levels[i-1]) < levels[i] || continue
        push!(runs, start:(i-1))
        start = i
    end
    push!(runs, start:lastindex(levels))
    return runs
end

"""
    _halfway(a, b)

Return a midpoint strictly between `a` and `b`, or `nothing` if arithmetic cannot produce one.
"""
function _halfway(a, b)
    mid = a + (b - a) / 2
    return a < mid < b ? mid : nothing
end

"""
    _wider_than(T)

Types to try for edges around draws of type `T`, widest last.
"""
_wider_than(::Type{T}) where {T} = (promote_type(T, Float64), BigFloat)

"""
    _one_bin(lo, hi, bins=1)

Return edges for `bins` bins covering `[lo, hi]`, or `nothing` if they do not fit.
"""
function _one_bin(lo, hi, bins = 1)
    for T in _wider_than(typeof(lo))
        a, b = T(lo), T(hi)
        # Allow enough room for distinct edges when a count was requested.
        room = bins * eps(max(abs(a), abs(b), one(T)))
        reach = max(T(0.5), room)
        low = min(a - reach, prevfloat(a))
        high = max(b + reach, nextfloat(b))
        (isfinite(low) && isfinite(high)) || continue
        edges = range(low, high; length = bins + 1)
        _edges_hold(edges, a, b) && return edges
    end
    return nothing
end

"""
    _bins_around(levels, runs, lo, hi)

Return one bin per run when suitable edges are found, or `nothing`.
"""
function _bins_around(levels, runs, lo, hi)
    n = length(runs)
    for T in _wider_than(eltype(levels))
        at = convert(Vector{T}, levels)
        # Use run boundaries so an edge cannot cut through a run.
        mids = [_halfway(at[last(runs[i])], at[first(runs[i+1])]) for i = 1:(n-1)]
        any(isnothing, mids) && continue

        # Fall back to the next representable edge at the type's limits.
        low = at[first(runs[1])] - (mids[1] - at[last(runs[1])])
        high = at[last(runs[n])] + (at[first(runs[n])] - mids[n-1])
        isfinite(low) || (low = prevfloat(at[first(runs[1])]))
        isfinite(high) || (high = nextfloat(at[last(runs[n])]))

        edges = vcat(low, mids, high)
        # Checking only the extrema could miss an edge inside a run.
        all(edges[i] <= at[first(runs[i])] && at[last(runs[i])] < edges[i+1] for i = 1:n) ||
            continue
        _edges_hold(edges, T(lo), T(hi)) && return edges
    end
    return nothing
end

"""
    _edges_hold(edges, lo, hi)

Check finite edges that hold `[lo, hi]` under left-closed bins.
"""
function _edges_hold(edges, lo, hi)
    all(isfinite, edges) || return false
    # A finite step can still span a width that breaks range bin lookup.
    if edges isa AbstractRange
        (isfinite(step(edges)) && isfinite(last(edges) - first(edges))) || return false
    end
    issorted(edges, lt = <=) || return false
    return first(edges) <= lo && hi < last(edges)
end

"""
    _few_levels(values)

Return sorted distinct values if the level and draw counts qualify, or `nothing`.
"""
function _few_levels(values)
    levels = Set{eltype(values)}()
    for v in values
        # Set uses `isequal`, which treats signed zeros as distinct values.
        push!(levels, iszero(v) ? zero(v) : v)
        length(levels) > MAX_DISCRETE_LEVELS && return nothing
    end
    # Count draws across all levels so rare levels still qualify.
    (1 < length(levels) && 5 * length(levels) <= length(values)) || return nothing
    return sort!(collect(levels))
end

"""
    _bin_count(values, rule)

How many bins to cut `values` into, for the binning rules Plots names.

`:sqrt`, `:sturges`, `:rice`, `:scott` and `:fd` are the rules Plots implements, computed the
same way here so that passing `bins` keeps meaning what it always did. `:auto` is `:fd`.
"""
function _bin_count(values, rule::Symbol)
    rule in BIN_RULES || throw(
        ArgumentError(
            "`bins` must be a count, a set of edges, or one of $(join(BIN_RULES, ", "))",
        ),
    )

    n = length(values)
    n > 1 || return 1
    lo, hi = extrema(values)
    hi > lo || return 1
    # A Float16 or Float32 span can overflow in its own type.
    W = promote_type(eltype(values), Float64)
    span = W(hi) - W(lo)
    # Clamp before converting a near-zero bin width to an integer count.
    count(x) = isnan(x) ? 1 : clamp(ceil(Int, clamp(x, 1.0, 10_000.0)), 1, 10_000)

    rule === :auto && (rule = :fd)
    nd = cbrt(n)
    if rule === :sqrt
        return count(sqrt(n))
    elseif rule === :sturges
        return count(log2(n) + 1)
    elseif rule === :rice
        return count(2 * nd)
    elseif rule === :scott
        # Squaring narrow floats can overflow before their standard deviation does.
        spread = std(W(v) for v in values)
        iszero(spread) && return 1
        denominator = 3.5 * spread
        estimate = span * nd / denominator
        if all(isfinite, (span, spread, denominator, estimate))
            return count(estimate)
        end
        wide_spread = std(BigFloat(v) for v in values)
        return count((BigFloat(hi) - BigFloat(lo)) * nd / (BigFloat(3.5) * wide_spread))
    end
    lower, upper = quantile(values, 0.25), quantile(values, 0.75)
    iqr = W(upper) - W(lower)
    denominator = 2 * (iqr > 0 ? iqr : one(iqr))
    estimate = span * nd / denominator
    if all(isfinite, (span, iqr, denominator, estimate))
        return count(estimate)
    end
    wide_iqr = BigFloat(upper) - BigFloat(lower)
    wide_denominator = 2 * (wide_iqr > 0 ? wide_iqr : one(wide_iqr))
    return count((BigFloat(hi) - BigFloat(lo)) * nd / wide_denominator)
end

"""
    _histogram_edges(series, bins)

Bin edges covering every chain of one parameter.

Binning each chain over its own range gives bars of different widths, which cannot be
compared, so the edges come from the chains pooled. `bins` is a rule, a count, or edges to
use as they are. The count asked for is the count returned: the top edge is nudged past the
largest draw, which the last bin would otherwise leave out for being open on the right.

A rule this does not compute, such as `:wand`, is handed back for Plots to apply. Those
chains are then binned one at a time, which is what Plots did with all of them before.

Draws that are not finite are left out of the range, the way Plots leaves them out of its
own.
"""
function _histogram_edges(series, bins)
    bins isa Symbol && bins ∉ BIN_RULES && return bins
    pooled = [float(v) for v in Iterators.flatten(series) if !ismissing(v) && isfinite(v)]
    if bins isa Symbol
        if bins === :auto
            discrete = _discrete_edges(pooled)
            discrete === nothing || return _announced(discrete, pooled)
        end
        bins = _bin_count(pooled, bins)
    end
    bins isa Integer || return bins

    isempty(pooled) && return range(-0.5, 0.5; length = bins + 1)
    first_draw, last_draw = extrema(pooled)
    # Half a unit cannot widen a bin around a sufficiently large draw.
    if last_draw <= first_draw
        single = _one_bin(first_draw, last_draw, bins)
        single === nothing && return bins
        return _announced(single, pooled)
    end

    for T in _wider_than(eltype(pooled))
        lo, hi = T(first_draw), T(last_draw)
        # The last bin is open on the right.
        top = nextfloat(hi)
        isfinite(top) || continue
        edges = range(lo, top; length = bins + 1)
        _edges_hold(edges, lo, hi) && return _announced(edges, pooled)
        # A wide range can round away the extra step at the last edge.
        nudged = collect(range(lo, hi; length = bins + 1))
        nudged[end] = top
        _edges_hold(nudged, lo, hi) && return _announced(nudged, pooled)
    end
    return bins
end

"""
    _announced(edges, pooled)

Return `edges` and warn when their `Float64` range cannot hold the draws or has an infinite span.
"""
function _announced(edges, pooled)
    lo, hi = extrema(pooled)
    at_canvas = map(Float64, edges)
    usable =
        _edges_hold(at_canvas, Float64(lo), Float64(hi)) &&
        isfinite(last(at_canvas) - first(at_canvas))
    usable || @warn(
        "Histogram range may exceed the plot axis. Draws were counted, but the plot may be empty.",
        extrema = (lo, hi),
        maxlog = 1,
    )
    return edges
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

# MCMCDiagnosticTools splits a chain in two and needs more than four draws in each half, so
# a shorter prefix than this has no effective sample size to report.
const MIN_EVOLUTION_DRAWS = 10

"""
    _evolution_points(ndraws, npoints)

Draw counts at which to recompute a diagnostic, spread evenly and always ending at `ndraws`.

Nothing below `MIN_EVOLUTION_DRAWS`, because the effective sample size is not defined on a
prefix that short and asking for it only returns a warning and a `NaN`.
"""
function _evolution_points(ndraws::Integer, npoints::Integer)
    npoints = max(2, npoints)
    first_point = min(MIN_EVOLUTION_DRAWS, ndraws)
    step = max(1, ndraws ÷ npoints)
    points = filter(>=(first_point), collect(step:step:ndraws))
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
        present = filter(isfinite, col)
        centre = isempty(present) ? zero(eltype(out)) : mean(present)
        spread = length(present) > 1 ? std(present) : zero(eltype(out))
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
