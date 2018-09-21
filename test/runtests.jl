using Test
using Pkg

# add packages
to_add = [
    PackageSpec(name="Plots"),
    PackageSpec(name="StatPlots"),
]

Pkg.add(to_add)

# run plotting tests
include("plot_test.jl")

# run function tests
include("diagnostic_tests.jl")

# run tests for missing values
include("missing_tests.jl")
