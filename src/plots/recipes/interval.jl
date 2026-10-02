@recipe function f(
    p::RidgelinePlot;
    hpd_val = [0.05, 0.2],
    q = [0.1, 0.9],
    spacer = 0.5,
    _riser = 0.2,
    show_mean = true,
    show_median = true,
    show_qi = false,
    show_hpdi = true,
    fill_q = true,
    fill_hpd = false,
    ordered = false,
)

    chn = p.args[1]
    par_names = p.args[2]

    for i = 1:length(par_names)
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
        q_int = _compute_plot_data(
            i,
            chn,
            par_names;
            hpd_val = hpd_val,
            q = q,
            spacer = spacer,
            _riser = _riser,
            show_mean = show_mean,
            show_median = show_median,
            show_qi = show_qi,
            show_hpdi = show_hpdi,
            fill_q = fill_q,
            fill_hpd = fill_hpd,
            ordered = ordered,
        )

        yticks --> (
            length(par_names) > 1 ?
            (_riser .+ ((1:length(par_names)) .- 1) .* spacer, string.(par)) : :default
        )
        yaxis --> (length(par_names) > 1 ? "Parameters" : "Density")
        @series begin
            seriestype := :hline
            label := nothing
            linecolor := "#BBBBBB"
            linewidth --> 1.2
            [h]
        end
        @series begin
            seriestype := :path
            label := nothing
            fillrange --> min
            fillalpha --> 0.8
            x_int, val
        end
        @series begin
            seriestype := :path
            label := nothing
            linecolor --> "#000000"
            k_density.x, k_density.density .+ h
        end
        @series begin
            seriestype := :path
            label --> (show_mean ? (i == 1 ? "Mean" : nothing) : nothing)
            linecolor --> "dark red"
            linewidth --> (show_mean ? 1.2 : 0)
            [chain_mean, chain_mean], [min, min + pdf(k_density, chain_mean)]
        end
        @series begin
            seriestype := :path
            label --> (show_median ? (i == 1 ? "Median" : nothing) : nothing)
            linecolor --> "#000000"
            linewidth --> (show_median ? 1.2 : 0)
            [chain_med, chain_med], [min, min + pdf(k_density, chain_med)]
        end
        @series begin
            seriestype := :scatter
            label := (show_qi ? (i == 1 ? "Q$(q[1]), Q$(q[2])" : nothing) : nothing)
            markershape --> (show_qi ? :diamond : :circle)
            markercolor --> "#000000"
            markersize --> (show_qi ? 2 : 0)
            q_int, [h]
        end
        @series begin
            seriestype := :path
            label := nothing
            linecolor := "#000000"
            linewidth --> (show_qi ? 1.2 : 0)
            [qs[1], qs[2]], [h, h]
        end
        @series begin
            seriestype := :path
            label := (
                show_hpdi ? (i == 1 ? "$(Integer((1-hpdi[1])*100))% HPDI" : nothing) :
                nothing
            )
            linewidth --> (show_hpdi ? 2 : 0)
            seriesalpha --> 0.80
            linecolor --> :darkblue
            [lower_hpd[1][1], upper_hpd[1][1]], [h, h]
        end
    end
end

@recipe function f(
    p::ForestPlot;
    hpd_val = [0.05, 0.2],
    q = [0.1, 0.9],
    spacer = 0.5,
    _riser = 0.2,
    show_mean = true,
    show_median = true,
    show_qi = false,
    show_hpdi = true,
    fill_q = true,
    fill_hpd = false,
    ordered = false,
)

    chn = p.args[1]
    par_names = p.args[2]

    for i = 1:length(par_names)
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
        q_int = _compute_plot_data(
            i,
            chn,
            par_names;
            hpd_val = hpd_val,
            q = q,
            spacer = spacer,
            _riser = _riser,
            show_mean = show_mean,
            show_median = show_median,
            show_qi = show_qi,
            show_hpdi = show_hpdi,
            fill_q = fill_q,
            fill_hpd = fill_hpd,
            ordered = ordered,
        )

        yticks --> (
            length(par_names) > 1 ?
            (_riser .+ ((1:length(par_names)) .- 1) .* spacer, string.(par)) : :default
        )
        yaxis --> (length(par_names) > 1 ? "Parameters" : "Density")

        for j = 1:length(hpdi)
            @series begin
                seriestype := :path
                label := (
                    show_hpdi ?
                    (i == 1 ? "$(Integer((1-hpdi[j])*100))% HPDI" : nothing) : nothing
                )
                linecolor --> j
                linewidth --> (show_hpdi ? 1.5 * j : 0)
                seriesalpha --> 0.80
                [lower_hpd[j][1], upper_hpd[j][1]], [h, h]
            end
        end
        @series begin
            seriestype := :scatter
            label := (show_median ? (i == 1 ? "Median" : nothing) : nothing)
            markershape --> :diamond
            markercolor --> "#000000"
            markersize --> (show_median ? length(hpdi) : 0)
            [chain_med], [h]
        end
        @series begin
            seriestype := :scatter
            label := (show_mean ? (i == 1 ? "Mean" : nothing) : nothing)
            markershape --> :circle
            markercolor --> :gray
            markersize --> (show_mean ? length(hpdi) : 0)
            [chain_mean], [h]
        end
        @series begin
            seriestype := :scatter
            label := (show_qi ? (i == 1 ? "Q1 = $(q[1]), Q3 = $(q[2])" : nothing) : nothing)
            markershape --> (show_qi ? :diamond : :circle)
            markercolor --> "#000000"
            markersize --> (show_qi ? 2 : 0)
            q_int, [h]
        end
        @series begin
            seriestype := :path
            label := nothing
            linecolor := "#000000"
            linewidth --> (show_qi ? 1.2 : 0.0)
            [qs[1], qs[2]], [h, h]
        end
    end
end
