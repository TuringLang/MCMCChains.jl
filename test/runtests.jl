using Test

@testset "MCMCChains" begin
    # run tests for effective sample size
    include("ess_tests.jl")

    # run plotting tests
    println("Plotting")
    @time include("plot_test.jl")

    # run function tests
    println("Diagnostics")
    @time include("diagnostic_tests.jl")

    # run tests for missing values
    println("Missing values")
    @time include("missing_tests.jl")

    # run tests for describing sections
    println("Sections")
    @time include("sections_tests.jl")

    # run tests for accessing parameters
    println("Accessing parameters")
    @time include("get_tests.jl")

    # run tests for serialization
    println("Serialization")
    @time include("serialization_tests.jl")

<<<<<<< HEAD
    # run tests for array constructor
    include("arrayconstructor_tests.jl")
=======
    # run tests for sampling api
    println("Sampling")
    @time include("sampling_tests.jl")
>>>>>>> 338784ee4d17478724037ab151e4ec1cf437899c

    # run tests for array constructor
    println("Array")
    @time include("arrayconstructor_tests.jl")

    # run tests for tables interfaces
    println("Tables interfaces")
    @time include("tables_tests.jl")

<<<<<<< HEAD
    # run tests for summaries
    include("summarize_tests.jl")
=======
    # run tests for dataframe summary
    println("Summary")
    @time include("summarize_tests.jl")
>>>>>>> 338784ee4d17478724037ab151e4ec1cf437899c

    # run tests for posterior stats
    println("Model statistics")
    @time include("modelstats_test.jl")

    # run tests for concatenation
    println("Concatenation")
    @time include("concatenation_tests.jl")
end
