using Chain
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

include("plot_test.jl")
