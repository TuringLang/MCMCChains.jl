function p_direction(x::Vector{Float64})
    return maximum([sum(x .> 0) ./ length(x), sum(x .< 0) ./ length(x)])
end

function p_direction(chains::Chains; kwargs...)
    # Store everything.
    funs = [p_direction]
    func_names = [:p_direction]

    # Summarize.
    summary_df = summarize(
        chains, funs...;
        func_names=func_names,
        name="Probability of Direction (pd)",
        kwargs...
    )

    return summary_df
end
