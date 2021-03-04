export dic

#################### Posterior Statistics ####################

"""
    dic(chain::Chains, loglik::Symbol)

Compute the deviance information criterion (DIC) from Chain object on posterior log likelihood samples specified by varible name loglik.

Note: DIC assumes that the posterior distribution is approx. multivariate Gaussian and tends to select overfitted models.
"""
function StatsModelComparisons.dic(chain::Chains, loglik::Symbol)
    lps = Array(chain[:, loglik, :])
    return dic(vec(lps))
end
