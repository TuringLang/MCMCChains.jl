@recipe function f(
    p::PPCPlot;
    kind = :density,
    alpha = nothing,
    num_pp_samples = nothing,
    mean_pp = true,
    observed = nothing,
    observed_rug = false,
    colors = [:steelblue, :black, :orange],
    jitter = nothing,
    legend = true,
    random_seed = nothing,
    ppc_group = :posterior,
)
    _apply_chrome!(plotattributes)

    if length(p.args) < 3
        error(
            "ppcplot requires at least 3 arguments: (posterior_chains, posterior_predictive_chains, observed_data)",
        )
    end

    posterior_chains = p.args[1]
    pp_chains = p.args[2]
    observed_data = p.args[3]

    if !(posterior_chains isa Chains)
        error("First argument must be a Chains object (posterior chains)")
    end
    if !(pp_chains isa Chains)
        error("Second argument must be a Chains object (posterior predictive chains)")
    end
    if !(observed_data isa AbstractVector)
        error("Third argument must be a vector (observed data)")
    end

    if kind ∉ (:density, :histogram, :scatter, :cumulative)
        error(
            "`kind` must be one of `:density`, `:histogram`, `:scatter`, or `:cumulative`",
        )
    end

    if ppc_group ∉ (:posterior, :prior)
        error("`ppc_group` must be one of `:posterior` or `:prior`")
    end

    if observed === nothing
        observed = (ppc_group == :posterior)
    end

    if length(colors) != 3
        error(
            "`colors` must be a vector of length 3: [predictive_color, observed_color, mean_color]",
        )
    end

    if alpha === nothing
        alpha = (kind == :scatter) ? 0.7 : 0.2
    end

    if jitter === nothing
        jitter = 0.0
    end

    pp_pooled = pool_chain(pp_chains)
    pp_data = Array(pp_pooled.value.data)
    pp_data = reshape(pp_data, size(pp_data, 1) * size(pp_data, 3), size(pp_data, 2))

    if random_seed !== nothing
        Random.seed!(random_seed)
    end

    total_pp_samples = size(pp_data, 1)
    if num_pp_samples === nothing
        if kind == :scatter
            num_pp_samples = min(5, total_pp_samples)
            if jitter == 0.0 && num_pp_samples <= 5
                jitter = 0.7
            end
        else
            num_pp_samples = total_pp_samples
        end
    end

    if num_pp_samples > total_pp_samples
        @warn "Requested $num_pp_samples samples but only $total_pp_samples available. Using all samples."
        num_pp_samples = total_pp_samples
    end

    if num_pp_samples < total_pp_samples
        sample_indices = Random.randperm(total_pp_samples)[1:num_pp_samples]
        pp_data = pp_data[sample_indices, :]
    end

    if ppc_group == :prior
        title := "Prior Predictive Check"
        predictive_label = "Prior Predictive"
        mean_label = "Prior Predictive Mean"
    else
        title := "Posterior Predictive Check"
        predictive_label = "Posterior Predictive"
        mean_label = "Posterior Predictive Mean"
    end
    legend := legend

    if kind == :density
        xaxis := "Value"
        yaxis := "Density"

        for i = 1:size(pp_data, 1)
            @series begin
                seriestype := :density
                label := i == 1 ? predictive_label : ""
                color := colors[1]
                alpha := alpha
                linewidth := 1
                pp_data[i, :]
            end
        end

        if observed
            @series begin
                seriestype := :density
                label := "Observed"
                color := colors[2]
                linewidth := 2
                observed_data
            end
        end

        if mean_pp
            pp_mean = vec(mean(pp_data, dims = 1))
            @series begin
                seriestype := :density
                label := mean_label
                color := colors[3]
                linewidth := 2
                linestyle := :dash
                pp_mean
            end
        end

        if observed_rug && observed
            y_min = 0
            @series begin
                seriestype := :scatter
                label := ""
                color := colors[2]
                markershape := :vline
                markersize := 2
                y := fill(y_min, length(observed_data))
                observed_data
            end
        end

    elseif kind == :histogram
        xaxis := "Value"
        yaxis := "Frequency"

        for i = 1:size(pp_data, 1)
            @series begin
                seriestype := :histogram
                label := i == 1 ? predictive_label : ""
                color := colors[1]
                alpha := alpha
                normalize := true
                bins := 30
                pp_data[i, :]
            end
        end

        if observed
            @series begin
                seriestype := :histogram
                label := "Observed"
                color := colors[2]
                alpha := 0.7
                normalize := true
                bins := 30
                observed_data
            end
        end

        if mean_pp
            pp_mean = vec(mean(pp_data, dims = 1))
            @series begin
                seriestype := :histogram
                label := mean_label
                color := colors[3]
                alpha := 0.7
                normalize := true
                bins := 30
                pp_mean
            end
        end

    elseif kind == :cumulative
        xaxis := "Value"
        yaxis := "Cumulative Probability"

        all_data = observed ? vcat(vec(pp_data), observed_data) : vec(pp_data)
        x_range = range(minimum(all_data), maximum(all_data), length = 200)

        for i = 1:size(pp_data, 1)
            pp_ecdf = ecdf(pp_data[i, :])
            y_vals = pp_ecdf.(x_range)
            @series begin
                seriestype := :path
                label := i == 1 ? predictive_label : ""
                color := colors[1]
                alpha := alpha
                linewidth := 1
                x := x_range
                y := y_vals
                ()
            end
        end

        if observed
            obs_ecdf = ecdf(observed_data)
            obs_y_vals = obs_ecdf.(x_range)
            @series begin
                seriestype := :path
                label := "Observed"
                color := colors[2]
                linewidth := 2
                x := x_range
                y := obs_y_vals
                ()
            end
        end

        if mean_pp
            pp_mean = vec(mean(pp_data, dims = 1))
            pp_mean_ecdf = ecdf(pp_mean)
            mean_y_vals = pp_mean_ecdf.(x_range)
            @series begin
                seriestype := :path
                label := mean_label
                color := colors[3]
                linewidth := 2
                linestyle := :dash
                x := x_range
                y := mean_y_vals
                ()
            end
        end

        if observed_rug && observed
            @series begin
                seriestype := :scatter
                label := ""
                color := colors[2]
                markershape := :vline
                markersize := 2
                x := observed_data
                y := fill(0, length(observed_data))
                ()
            end
        end

    elseif kind == :scatter
        xaxis := "Index"
        yaxis := "Value"

        if jitter > 0
            jitter_vals = jitter * (rand(size(pp_data, 2)) .- 0.5)
        else
            jitter_vals = zeros(size(pp_data, 2))
        end

        for i = 1:size(pp_data, 1)
            y_vals = pp_data[i, :] .+ jitter_vals[1:length(pp_data[i, :])]
            @series begin
                seriestype := :scatter
                label := i == 1 ? predictive_label : ""
                color := colors[1]
                alpha := alpha
                markersize := 3
                x := 1:length(y_vals)
                y := y_vals
                ()
            end
        end

        if observed
            obs_y = observed_data .+ jitter_vals[1:length(observed_data)]
            @series begin
                seriestype := :scatter
                label := "Observed"
                color := colors[2]
                markersize := 4
                markershape := :diamond
                x := 1:length(obs_y)
                y := obs_y
                ()
            end
        end

        if mean_pp
            pp_mean = vec(mean(pp_data, dims = 1))
            mean_y = pp_mean .+ jitter_vals[1:length(pp_mean)]
            @series begin
                seriestype := :scatter
                label := mean_label
                color := colors[3]
                markersize := 4
                markershape := :square
                x := 1:length(mean_y)
                y := mean_y
                ()
            end
        end
    end
end
