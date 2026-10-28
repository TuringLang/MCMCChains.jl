@testset "add_parameters" begin
    chn = Chains(
        reshape(1:48, 8, 3, 2),
        [:a, :b, :c];
        start = 101,
        thin = 2,
        evidence = 3.0,
        info = (source = :test,),
    )
    selected = chn[2:5, :, 2:2]
    scalar = selected[:a] .+ selected[:b]
    vector = map(x -> [x, 2x], scalar)

    extended = add_parameters(selected, :total => scalar, :prediction => vector)
    @test size(extended) == (4, 6, 1)
    @test names(extended) ==
          [:a, :b, :c, :total, Symbol("prediction[1]"), Symbol("prediction[2]")]
    @test range(extended) == 103:2:109
    @test chains(extended) == 2:2
    @test extended[:total] == scalar
    @test extended[Symbol("prediction[1]")] == scalar
    @test extended[Symbol("prediction[2]")] == 2 .* scalar
    @test extended.logevidence == selected.logevidence
    @test extended.info == selected.info
    @test extended.name_map == (parameters = names(extended),)
    @test names(selected) == [:a, :b, :c]

    sectioned = set_section(selected, (parameters = [:a, :b], internals = [:c]))
    sectioned_extended = add_parameters(sectioned, :total => scalar)
    @test sectioned_extended.name_map == (parameters = [:a, :b, :total], internals = [:c])
    @test sectioned.name_map == (parameters = [:a, :b], internals = [:c])

    cube = cat(reshape(scalar, 4, 1, 1), reshape(2 .* scalar, 4, 1, 1); dims = 2)
    from_cube = add_parameters(selected, "prediction" => cube)
    @test from_cube[Symbol("prediction[1]")] == scalar
    @test from_cube[Symbol("prediction[2]")] == 2 .* scalar

    generated = map(x -> (total = x, prediction = [x, 2x]), scalar)
    from_generated = add_parameters(selected, generated)
    @test from_generated.value == extended.value

    partly_missing = Matrix{Union{Missing,Float64}}(scalar)
    partly_missing[1, 1] = missing
    with_missing = add_parameters(selected, :maybe => partly_missing)
    @test ismissing(with_missing[:maybe][1, 1])
    @test with_missing[:maybe][2:end, :] == partly_missing[2:end, :]

    @test_throws DimensionMismatch add_parameters(selected, :bad => rand(3, 1))
    @test_throws DimensionMismatch add_parameters(selected, :bad => rand(4, 2, 2))
    @test_throws DimensionMismatch add_parameters(selected, generated[1:3, :])
    @test_throws ArgumentError add_parameters(selected, :bad => rand(4))
    @test_throws ArgumentError add_parameters(selected, :a => scalar)
    @test_throws ArgumentError add_parameters(selected, :total => scalar, :total => scalar)
    @test_throws ArgumentError add_parameters(
        selected,
        :bad => reshape(["a", "b", "c", "d"], 4, 1),
    )
    inconsistent = reshape([[1], [2, 3], [4], [5]], 4, 1)
    @test_throws ArgumentError add_parameters(selected, :bad => inconsistent)
    @test_throws ArgumentError add_parameters(selected)
end
