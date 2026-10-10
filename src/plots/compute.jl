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
