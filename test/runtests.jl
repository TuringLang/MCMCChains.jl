using Test

@time @testset "Diagnostics" begin include("diagnostic_tests.jl") end
@time @testset "Plotting" begin include("plot_test.jl") end
@time @testset "Missing values" begin include("missing_tests.jl") end
@time @testset "Sections" begin include("sections_tests.jl") end
@time @testset "Accessing parameters" begin include("get_tests.jl") end
@time @testset "Serialization" begin include("serialization_tests.jl") end
@time @testset "Sampling" begin include("sampling_tests.jl") end
@time @testset "Array constructor" begin include("arrayconstructor_tests.jl") end
@time @testset "Dataframe constructor" begin include("dfconstructor_tests.jl") end
@time @testset "Dataframe summary" begin include("summarize_tests.jl") end
@time @testset "Posterior" begin include("modelstats_test.jl") end
@time @testset "Concatenation" begin include("concatenation_tests.jl") end