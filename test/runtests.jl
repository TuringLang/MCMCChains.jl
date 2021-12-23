using Pkg

# Activate test environment on older Julia versions
if VERSION < v"1.2"
    Pkg.activate(@__DIR__)
    Pkg.develop(PackageSpec(path=dirname(@__DIR__)))
    Pkg.instantiate()
end

using MCMCChains
using Documenter

using Test
using Random

# set seed for all testsets
Random.seed!(0)

@testset "MCMCChains" begin
    # MLJXGBoostInterface requires Julia >= 1.3
    # Array printing depends on Julia version, therefore we only test Julia >= 1.7
    # XGBoost errors on 32bit systems: https://github.com/dmlc/XGBoost.jl/issues/92
    if VERSION >= v"1.7" && Sys.WORD_SIZE == 64
        # run tests related to rstar statistic
        println("Rstar")
        Pkg.add([
            Pkg.PackageSpec(; name="MLJBase", version="0.18"),
            Pkg.PackageSpec(; name="MLJXGBoostInterface", version="0.1.5"),
        ])
        @time include("rstar_tests.jl")

        DocMeta.setdocmeta!(
            MCMCChains,
            :DocTestSetup,
            :(using MCMCChains);
            recursive=true
        )

        doctest(MCMCChains)
    end

    # run tests for effective sample size
    println("ESS")
    @time include("ess_tests.jl")

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

    # run tests for sampling api
    println("Sampling")
    @time include("sampling_tests.jl")

    # run tests for array constructor
    println("Array")
    @time include("arrayconstructor_tests.jl")

    # run tests for tables interfaces
    println("Tables interfaces")
    @time include("tables_tests.jl")

    # run tests for dataframe summary
    println("Summary")
    @time include("summarize_tests.jl")

    # run tests for posterior stats
    println("Model statistics")
    @time include("modelstats_test.jl")

    # run tests for concatenation
    println("Concatenation")
    @time include("concatenation_tests.jl")

end
