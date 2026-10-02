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

    return par,
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
    q_int
end
