export dic

#################### Posterior Statistics ####################

"""
    dic(chain::Chains)

Computes Deviance Information Criterion (DIC).

# Arguments
* `chain::Chains: an MCMC `Chain` object containing posterior log likelihood values
* `loglik::Symbol`: variable name associated with posterior log likelihood values

# Returns
* `dic::Real`: DIC value
"""
function dic(chain::Chains, loglik::Symbol)
    lps = Array(chain[:, loglik, :])
    return dic(vec(lps))
end
