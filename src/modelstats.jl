export dic

#################### Posterior Statistics ####################

"""
    dic(chain::Chains, loglik::Symbol)

Compute the deviance information criterion (DIC) from `chain` on posterior log likelihood samples specified by parameter name `loglik`.

Note: DIC assumes that the posterior distribution is approx. multivariate Gaussian and tends to select overfitted models.
"""
function StatsModelComparisons.dic(chain::Chains, loglik::Symbol)
    lps = Array(chain[:, loglik, :])
    return dic(vec(lps))
end
