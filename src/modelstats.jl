export dic

#################### Posterior Statistics ####################

"""
    dic(chain::Chains, logpdf::Function) -> (DIC, pD)

Compute the deviance information criterion (DIC).
(Smaller is better)

Note: DIC assumes that the posterior distribution is approx. multivariate Gaussian and tends to select overfitted models.

## Returns:

* `DIC`: The calculated deviance information criterion
* `pD`: The effective number of parameters

## Usage:

```
chn ... # sampling results
lpfun = function f(chain::Chains) # function to compute the logpdf values
    niter, nparams, nchains = size(chain)
    lp = zeros(niter + nchains) # resulting logpdf values
    for i = 1:nparams
        lp += map(p -> logpdf( ... , x), Array(chain[:,i,:]))
    end
    return lp
end

DIC, pD = dic(chn, lpfun)
```
"""
function dic(chain::Chains, logpdf::Function)

    # expectation of each parameter
    Eθ = reshape(mean(Array(chain), dims = [1,3]), 1,:,1)
    Echain = Chains(Eθ)
    EθD = -2*mean(logpdf(Echain))

    D = -2*logpdf(chain)
    ED = mean(D)

    pD = 2*(ED - EθD)

    DIC = EθD + pD

    return DIC, pD
end
