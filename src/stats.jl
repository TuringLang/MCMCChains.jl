#################### Posterior Statistics ####################
function autocor(chn::AbstractChains;
                 lags::Vector=[1, 5, 10, 50],
                 relative::Bool=true,
                 showall=false,
                 suppress_header=false,
                 section=:parameters)

    # Allocation of summary vector.
    summaries = Vector{ChainSummary}([])

    # Separate the chain into sections if showall=false.
    chns = showall ? [chn] : get_sections(chn, section)

    # Set showlabels bool.
    show_labels = true

    # Iterate through each chain.
    for c in chns
        # If we're only going through one section at a time,
        # set a section name.
        section_name = showall ? "" : string.(first(keys(c.name_map)))

        if relative
            lags *= step(c)
        elseif any(lags .% step(c) .!= 0)
            throw(ArgumentError("lags do not correspond to thinning interval"))
        end
        labels = map(x -> "Lag " * string(x), lags)

        (niter, nvar, nchain) = size(c.value)
        vals = zeros(nvar, length(lags), nchain)
        for k in 1:nchain
            for v in 1:nvar
                # skipping missing values
                # skipping should be done inside of autocor to ensure correct alignment
                # TODO: modify autocor to deal with missing values
                x = convert(Vector{Float64}, collect(skipmissing(c.value[:,v,k])))
                vals[v, :, k] = autocor(x, lags)
            end
        end

        new_summary = ChainSummary(vals,
            string.(names(c)),
            labels,
            "",
            true)

        push!(summaries, new_summary)
    end

    h = suppress_header ? "" : header(chn)
    return ChainSummaries(h, summaries)
end

function cor(chn::AbstractChains;
    showall=false, suppress_header=false, section=:parameters)
    # Allocation of summary vector.
    summaries = Vector{ChainSummary}([])

    # Separate the chain into sections if showall=false.
    chns = showall ? [chn] : get_sections(chn, section)

    # Set showlabels bool.
    show_labels = true

    # Iterate through each chain.
    for c in chns
        # If we're only going through one section at a time,
        # set a section name.
        section_name = showall ? "" : string.(first(keys(c.name_map)))
        new_summary = ChainSummary(cor(combine(c)),
            string.(names(c)),
            string.(names(c)),
            section_name)
        push!(summaries, new_summary)
    end

    h = suppress_header ? "" : header(chn)
    return ChainSummaries(h, summaries)
end

function changerate(chn::AbstractChains;
    showall=false, suppress_header=false, section=:parameters)
    # Check for missing values.
    @assert !any(ismissing.(chn.value)) "Change rate comp. doesn't support missing values."

    # Allocation of summary vector.
    summaries = Vector{ChainSummary}([])

    # Separate the chain into sections if showall=false.
    chns = showall ? [chn] : get_sections(chn, section)

    # Set showlabels bool.
    show_labels = true

    # Iterate through each chain.
    for c in chns
        # If we're only going through one section at a time,
        # set a section name.
        section_name = showall ? "" : string.(first(keys(c.name_map)))

        n, p, m = size(c.value)
        r = zeros(Float64, p, 1, 1)
        r_mv = 0.0
        delta = Array{Bool}(undef, p)
        for k in 1:m
            prev = c.value[1, :, k]
            for i in 2:n
                for j in 1:p
                    x = c.value[i, j, k]
                    dx = x != prev[j]
                    r[j] += dx
                    delta[j] = dx
                    prev[j] = x
                end
                r_mv += any(delta)
            end
        end
        vals = round.([r; r_mv] / (m * (n - 1)), digits = 3)
        new_summary = ChainSummary(vals,
            [string.(names(c)); "Multivariate"],
            ["Change Rate"],
            section_name)
        push!(summaries, new_summary)
    end

    h = suppress_header ? "" : header(chn)
    return ChainSummaries(h, summaries)
end

describe(c::AbstractChains; args...) = describe(stdout, c; args...)

function describe(io::IO,
                  c::AbstractChains;
                  q = [0.025, 0.25, 0.5, 0.75, 0.975],
                  etype=:bm,
                  showall=false,
                  section=:parameters,
                  args...
                 )
    # Print the chain header.
    println(io, header(c))

    # Generate summary statistics.
    ps_stats = summarystats(c; etype=etype,
        showall=showall,
        suppress_header=true,
        section=section,
        args...)
    ps_quantiles = quantile(c,
        q=q,
        showall=showall,
        suppress_header=true,
        section=section)

    # Get linewidths.
    stats_linewidth = maximum(map(max_width, ps_stats.summaries))
    quantiles_linewidth = maximum(map(max_width, ps_quantiles.summaries))
    linewidth = max(quantiles_linewidth, stats_linewidth)

    # Print stats.
    print(io, "Empirical Posterior Estimates:\n")
    println(repeat("=", linewidth))
    show(io, ps_stats)

    print(io, "Quantiles:\n")
    println(repeat("=", linewidth))
    show(io, ps_quantiles)
end

function hpd(x::Vector{T}; alpha::Real=0.05) where {T<:Real}
    n = length(x)
    m = max(1, ceil(Int, alpha * n))

    y = sort(x)
    a = y[1:m]
    b = y[(n - m + 1):n]
    _, i = findmin(b - a)

    return [a[i], b[i]]
end

function hpd(chn::AbstractChains; alpha::Real=0.05,
    showall=false, suppress_header=false, section=:parameters)
    # Allocation summary vector.
    summaries = Vector{ChainSummary}([])

    # Separate the chain into sections if showall=false.
    chns = showall ? [chn] : get_sections(chn, section)

    # Set showlabels bool.
    show_labels = true

    # Iterate through each chain.
    for c in chns
        # If we're only going through one section at a time,
        # set a section name.
        section_name = showall ? "" : string.(first(keys(c.name_map)))

        pct = first(showoff([100.0 * (1.0 - alpha)]))
        labels = ["$(pct)% Lower", "$(pct)% Upper"]
        vals = permutedims(
            mapslices(x -> hpd(collect(skipmissing(x)), alpha=alpha), c.value, dims = [1, 3]),
            [2, 1, 3]
        )
        new_summary = ChainSummary(convert(Array{Float64,3}, vals),
            string.(names(c)), labels, section_name)
        push!(summaries, new_summary)
    end

    h = suppress_header ? "" : header(chn)
    return ChainSummaries(h, summaries)
end

function quantile(chn::AbstractChains; q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975],
    showall=false, suppress_header=false, section=:parameters)
    # Allocation summary vector.
    summaries = Vector{ChainSummary}([])

    # Separate the chain into sections if showall=false.
    chns = showall ? [chn] : get_sections(chn, section)

    # Set showlabels bool.
    show_labels = true

    # Iterate through each chain.
    for c in chns
        # If we're only going through one section at a time,
        # set a section name.
        section_name = showall ? "" : string.(first(keys(c.name_map)))

        # Make labels.
        labels = map(x -> string(100 * x) * "%", q)
        vals = Array(hcat([quantile(collect(skipmissing(c[i]))) for i in names(c)]...)')
        new_summary = ChainSummary(round.(vals, digits=4), string.(names(c)), labels, section_name, true)
        push!(summaries, new_summary)
    end

    h = suppress_header ? "" : header(chn)
    return ChainSummaries(h, summaries)
end

function summarystats(chn::AbstractChains; etype=:bm,
    showall=false, suppress_header=false, section=:parameters, args...)
    # Allocation summary vector.
    summaries = Vector{ChainSummary}([])

    # Summary statistics function array.
    f(x) = [mean(x),
            std(x),
            sem(x),
            mcse(x, etype; args...)]

    # Separate the chain into sections if showall=false.
    chns = showall ? [chn] : get_sections(chn, section)

    # Set showlabels bool.
    show_labels = true

    # Iterate through each chain.
    for c in chns
        # If we're only going through one section at a time,
        # set a section name.
        section_name = showall ? "" : string.(first(keys(c.name_map)))

        # Make labels.
        labels = if show_labels
            show_labels = false
            ["Mean", "SD", "Naive SE", "MCSE", "ESS"]
        else
            repeat([""], 5)
        end

        # Make summary statistics and ChainSummary.
        vals = hcat([f(collect(skipmissing(c[i]))) for i in names(c)]...)'
        stats = [vals  min.((vals[:, 2] ./ vals[:, 4]).^2, size(c.value, 1))]
        new_summary = ChainSummary(round.(stats, digits=4), string.(names(c)), labels, section_name, true)
        push!(summaries, new_summary)
    end

    h = suppress_header ? "" : header(chn)
    return ChainSummaries(h, summaries)
end
