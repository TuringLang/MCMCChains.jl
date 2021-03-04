export dic

#################### Posterior Statistics ####################

"""
    dic(chain::Chains, loglik::Symbol)

Compute the deviance information criterion (DIC).

Note: DIC assumes that the posterior distribution is approx. multivariate Gaussian and tends to select overfitted models.

# Arguments
* `chain::Chains: an MCMC `Chain` object containing posterior log likelihood values
* `loglik::Symbol`: variable name associated with posterior log likelihood values

"""
function StatsModelComparisons.dic(chain::Chains, loglik::Symbol)
    lps = Array(chain[:, loglik, :])
    return dic(vec(lps))
end
