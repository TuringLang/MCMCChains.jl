using Test
using Pkg

# add packages
to_add = [
    PackageSpec(name="Plots"),
    PackageSpec(name="StatsPlots"),
    PackageSpec(name="Turing"),
]

Pkg.add(to_add)

@testset "MCMCChains" begin
    # run plotting tests
    include("plot_test.jl")

    # run function tests
    include("diagnostic_tests.jl")

    # run tests for missing values
    include("missing_tests.jl")

    # run tests for missing values
    include("sections_tests.jl")

    # run tests for missing values
    include("serialization_tests.jl")

    # run tests for sampling api
    include("sampling_tests.jl")

    # run tests for array constructor
    include("arrayconstructor_tests.jl")

    # run tests for dataframe constructor
    include("dfconstructor_tests.jl")
end
