using MCMCChains
using JSON
using Test

@testset "JSON Serialization (Extension)" begin
    n_iter = 100
    n_params = 5
    n_chains = 2

    val = randn(n_iter, n_params, n_chains)

    param_names = ["a", "b", "c", "d", "e"]

    # Mixed types Metadata (String, Float, Symbol, Dict) to verify JSON.lower
    info = (
        model_name = "test_model",
        start_time = 123456789.0,
        sampler = "NUTS",
        tags = [:test, :json],
        config = Dict(:step_size => 0.1, :max_depth => 10),
    )

    chn = Chains(val, param_names, Dict(:internals => ["d", "e"]); info = info)

    @testset "String Export" begin
        json_str = write_json(chn; as_string = true)
        @test json_str isa String

        parsed = JSON.parse(json_str)
        @test parsed["info"]["model_name"] == "test_model"
        @test parsed["info"]["tags"] == ["test", "json"]
        @test length(parsed["parameters"]) == n_params
        @test parsed["size"] == [n_iter, n_params, n_chains]
    end

    @testset "File Export" begin
        mktempdir() do tmpdir
            cd(tmpdir) do
                file_path_default = write_json(chn)
                @test file_path_default == "test_model.json"
                @test isfile(file_path_default)

                chn_loaded = read_json(file_path_default)

                @test size(chn_loaded) == size(chn)
                @test names(chn_loaded) == Symbol.(param_names)
                @test chn_loaded.value.data ≈ chn.value.data

                @test :d in names(chn_loaded, :internals)
                @test :a in names(chn_loaded, :parameters)

                @test chn_loaded.info.model_name == "test_model"
                @test chn_loaded.info.sampler == "NUTS"

                custom_name = "custom_chain.json"
                write_json(chn, custom_name)
                @test isfile(custom_name)
                chn_custom = read_json(custom_name)
                @test chn_custom.value.data ≈ chn.value.data

                val_missing = convert(Array{Union{Float64,Missing}}, val)
                val_missing[1, 1, 1] = missing
                chn_missing = Chains(val_missing, param_names)

                missing_file = "missing.json"
                write_json(chn_missing, missing_file)
                chn_missing_loaded = read_json(missing_file)

                @test ismissing(chn_missing_loaded[1, 1, 1])
                @test chn_missing_loaded[2, 1, 1] ≈ chn_missing[2, 1, 1]
            end
        end
    end
end
