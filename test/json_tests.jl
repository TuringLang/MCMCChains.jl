using MCMCChains
using JSON
using Test

@testset "JSON Serialization (Extension)" begin
    n_iter = 100
    n_params = 5
    n_chains = 2

    val = randn(n_iter, n_params, n_chains)

    param_names = ["a", "b", "c", "d", "e"]

    info = (
        model_name = "test_model",
        start_time = 123456789.0,
        sampler = "NUTS",
        tags = [:test, :json],
        config = Dict(:step_size => 0.1, :max_depth => 10),
    )

    chn = Chains(val, param_names, Dict(:internals => ["d", "e"]); info = info)

    @testset "String Export" begin
        lowered = JSON.lower(chn)
        @test Set(keys(lowered)) == Set([
            :size,
            :value_flat,
            :iterations,
            :parameters,
            :chains,
            :logevidence,
            :name_map,
            :info,
        ])

        json_str_std = JSON.json(chn)
        @test json_str_std isa String
        parsed_std = JSON.parse(json_str_std)
        @test parsed_std["info"]["model_name"] == "test_model"
    end

    @testset "File Export" begin
        mktempdir() do tmpdir
            cd(tmpdir) do
                json_file_std = "std_output.json"
                JSON.json(json_file_std, chn)
                @test isfile(json_file_std)

                chn_loaded_std = JSON.parsefile(json_file_std, Chains)
                @test size(chn_loaded_std) == size(chn)
                @test names(chn_loaded_std) == Symbol.(param_names)
                @test chn_loaded_std.value.data ≈ chn.value.data

                @test :d in names(chn_loaded_std, :internals)
                @test :a in names(chn_loaded_std, :parameters)

                @test chn_loaded_std.info.model_name == "test_model"
                @test chn_loaded_std.info.sampler == "NUTS"

                loaded_tags = chn_loaded_std.info.tags
                @test length(loaded_tags) == 2
                @test string(loaded_tags[1]) == "test"
                @test string(loaded_tags[2]) == "json"

                loaded_config = chn_loaded_std.info.config
                @test loaded_config["step_size"] == 0.1
                @test loaded_config["max_depth"] == 10

                str_content = read(json_file_std, String)
                chn_from_str = JSON.parse(str_content, Chains)
                @test chn_from_str.value.data ≈ chn.value.data

                val_missing = convert(Array{Union{Float64,Missing}}, val)
                val_missing[1, 1, 1] = missing
                chn_missing = Chains(val_missing, param_names)

                missing_file = "missing.json"
                JSON.json(missing_file, chn_missing)
                chn_missing_loaded = JSON.parsefile(missing_file, Chains)

                @test ismissing(chn_missing_loaded[1, 1, 1])
                @test chn_missing_loaded[2, 1, 1] ≈ chn_missing[2, 1, 1]
            end
        end
    end
end
