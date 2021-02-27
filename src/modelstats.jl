export dic

#################### Posterior Statistics ####################

"""
    dic(chain::Chains)

Computes Deviance Information Criterion (DIC).

# Arguments
* `chain::Chains: an MCMC `Chain` object containing posterior log likelihood values

# Returns
* `dic::Real`: DIC value
"""
function dic(chain::Chains)
    lps = get(chain, :lp).lp[:] |> Array
    return dic(lps)
end
