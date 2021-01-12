function check_div(chn; per_chain=false)
    divergent = Array(chn[:numerical_error])
    
    n_for_chains = sum(divergent, dims=1)
    
    n_divergent = sum(n_for_chains)
    if n_divergent > 0
        N = length(divergent)
        @info "Warning: $(n_divergent) of $(N) iterations ended with a divergence ($(100 * n_divergent / N) %)"
        if per_chain
            chain_len, num_chains = size(divergent)
            for chain_num in 1:num_chains
                if n_for_chains[chain_num] > 0
                    @info "Chain $(chain_num): $(n_for_chains[chain_num]) of $(chain_len) iterations ended with a divergence ($(100 * n_for_chains[chain_num] / chain_len) %)"
                end
            end
        end
        
        @info "Try running with larger target acceptance rate to remove the divergences."
        return false
    else
        return true
    end
end

function check_n_eff(chn)
    no_warning = true
    
    n_steps = size(chn, 1)
    ess_df = ess(chn)
    for ix in size(ess_df)[1]
        param_name = string(ess_df[ix,:parameters])
        n_eff = ess_df[ix,:ess]
        ratio = n_eff / n_steps
        if ratio < 0.1 || isnan(ratio) || isinf(ratio)
            no_warning = false
            @info "n_eff / n_steps for parameter $(param_name) is $(ratio)"
        end
    end
    
    if no_warning
        @info "n_eff / n_iter looks reasonable for all parameters"
    else
        @info "n_eff / iter below 0.001 indicates that the effective sample size has likely been overestimated"
    end
    
    return no_warning
end

function check_rhat(chn)
    no_warning = true
    
    ess_df = ess(chn)
    for ix in size(ess_df)[1]
        param_name = string(ess_df[ix,:parameters])
        rhat = ess_df[ix,:rhat]
        if isnan(rhat) || isinf(rhat) || rhat > 1.1 || rhat < 0.9
            no_warning = false
            @info "Rhat for parameter $(param_name) is $(rhat)"
        end
    end
    
    if no_warning
        @info "Rhat looks reasonable for all parameters"
    end
    
    return no_warning
end

function check_diagnostics(chn)
    div = check_div(chn, per_chain=true)
    n_eff = check_n_eff(chn)
    rhat = check_rhat(chn)
    return div && n_eff && rhat
end