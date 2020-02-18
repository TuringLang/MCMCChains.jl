using Test

@testset "MCMCChains" begin
    # run function tests
    @time include("diagnostic_tests.jl")

    # run plotting tests
    @time include("plot_test.jl")

    # run tests for missing values
    @time include("missing_tests.jl")

    # run tests for describing sections
    @time include("sections_tests.jl")

    # run tests for accessing parameters
    @time include("get_tests.jl")

    # run tests for serialization
    @time include("serialization_tests.jl")

    # run tests for sampling api
    @time include("sampling_tests.jl")

    # run tests for array constructor
    @time include("arrayconstructor_tests.jl")

    # run tests for dataframe constructor
    @time include("dfconstructor_tests.jl")

    # run tests for dataframe summary
    @time include("summarize_tests.jl")

    # run tests for posterior stats
    @time include("modelstats_test.jl")
end
