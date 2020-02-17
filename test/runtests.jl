using Test

@testset "Plotting" begin include("plot_test.jl") end
@testset "Diagnostics" begin include("diagnostic_tests.jl") end
@testset "Missing values" begin include("missing_tests.jl") end
@testset "Sections" begin include("sections_tests.jl") end
@testset "Accessing parameters" begin include("get_tests.jl") end
@testset "Serialization" begin include("serialization_tests.jl") end
@testset "Sampling" begin include("sampling_tests.jl") end
@testset "Array constructor" begin include("arrayconstructor_tests.jl") end
@testset "Dataframe constructor" begin include("dfconstructor_tests.jl") end
@testset "Dataframe summary" begin include("summarize_tests.jl") end
@testset "Posterior" begin include("modelstats_test.jl") end
@testset "Concatenation" begin include("concatenation_tests.jl") end