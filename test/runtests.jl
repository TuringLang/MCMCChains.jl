using Test

@testset "MCMCChains" begin
    # run plotting tests
    include("plot_test.jl")

    # run function tests
    include("diagnostic_tests.jl")

    # run tests for missing values
    include("missing_tests.jl")

    # run tests for describing sections
    include("sections_tests.jl")

    # run tests for accessing parameters
    include("get_tests.jl")

    # run tests for serialization
    include("serialization_tests.jl")

    # run tests for sampling api
    include("sampling_tests.jl")

    # run tests for array constructor
    include("arrayconstructor_tests.jl")

    # run tests for dataframe constructor
    include("dfconstructor_tests.jl")

    # run tests for dataframe summary
    include("summarize_tests.jl")

    # run tests for posterior stats
    # include("modelstats_test.jl")
end
