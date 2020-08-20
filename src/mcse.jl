#################### Monte Carlo Standard Errors ####################

function mcse(x::Vector{<:Real}, method::Symbol=:imse; args...)
  method == :bm ? mcse_bm(x; args...) :
  method == :imse ? mcse_imse(x) :
  method == :ipse ? mcse_ipse(x) :
    throw(ArgumentError("unsupported mcse method $method"))
end

function mcse_bm(x::Vector{<:Real}; size::Integer=100)
    n = length(x)
    m = div(n, size)
    if m < 2
        @debug "iterations are < $(2 * size) and batch size is > $(div(n, 2))"
        return missing
    end
    mbar = [mean(x[i * size .+ (1:size)]) for i in 0:(m - 1)]
    return sem(mbar)
end

function mcse_imse(x::Vector{<:Real})
    n = length(x)
    m = div(n - 2, 2)
    x_ = map(Float64, x)
    ghat = autocov(x_, [0, 1])
    Ghat = sum(ghat)
    value = -ghat[1] + 2 * Ghat
    for i in 1:m
        Ghat = min(Ghat, sum(autocov(x_, [2 * i, 2 * i + 1])))
        Ghat > 0 || break
        value += 2 * Ghat
    end
    return sqrt(value / n)
end

function mcse_ipse(x::Vector{<:Real})
    n = length(x)
    m = div(n - 2, 2)
    x_ = map(Float64, x)
    ghat = autocov(x_, [0, 1])
    value = ghat[1] + 2 * ghat[2]
    for i in 1:m
        Ghat = sum(autocov(x_, [2 * i, 2 * i + 1]))
        Ghat > 0 || break
        value += 2 * Ghat
    end
    return sqrt(value / n)
end
