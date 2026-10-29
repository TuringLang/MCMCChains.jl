@testset "hcat with generated parameters" begin
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

    extended = hcat(selected, :total => scalar, :prediction => vector)
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
    sectioned_extended = hcat(sectioned, :total => scalar)
    @test sectioned_extended.name_map == (parameters = [:a, :b, :total], internals = [:c])
    @test sectioned.name_map == (parameters = [:a, :b], internals = [:c])

    cube = cat(reshape(scalar, 4, 1, 1), reshape(2 .* scalar, 4, 1, 1); dims = 2)
    from_cube = hcat(selected, "prediction" => cube)
    @test from_cube[Symbol("prediction[1]")] == scalar
    @test from_cube[Symbol("prediction[2]")] == 2 .* scalar

    generated = map(x -> (total = x, prediction = [x, 2x]), scalar)
    from_generated = hcat(selected, generated)
    @test from_generated.value == extended.value
    any_generated = reshape(Any[generated...], size(generated))
    @test hcat(selected, any_generated).value == extended.value
    one_draw = selected[1:1, :, :]
    @test hcat(one_draw, any_generated[1:1, :]) isa Chains
    @test hcat(one_draw, any_generated[1:1, :])[:total] == scalar[1:1, :]
    labeled_any_generated = MCMCChains.AxisArray(
        any_generated;
        iter = range(selected),
        chain = chains(selected),
    )
    @test hcat(selected, labeled_any_generated).value == extended.value
    reordered = Matrix{NamedTuple}(undef, size(scalar))
    for i in eachindex(scalar)
        x = scalar[i]
        reordered[i] =
            isodd(i) ? (total = x, prediction = [x, 2x]) : (prediction = [x, 2x], total = x)
    end
    @test hcat(selected, reordered).value == extended.value
    @test hcat(selected, reshape(Any[reordered...], size(reordered))).value ==
          extended.value
    mismatched_names = copy(reordered)
    mismatched_names[2] = (prediction = [1, 2], other = 3)
    @test_throws ArgumentError hcat(selected, mismatched_names)
    @test_throws ArgumentError hcat(
        selected,
        reshape(Any[mismatched_names...], size(mismatched_names)),
    )
    mixed_values = copy(any_generated)
    mixed_values[2] = 2
    @test_throws ArgumentError hcat(selected, mixed_values)
    @test_throws ArgumentError hcat(selected, fill(1.0, size(scalar)))

    labeled_scalar =
        MCMCChains.AxisArray(scalar; iter = range(selected), chain = chains(selected))
    @test hcat(selected, :labeled => labeled_scalar)[:labeled] == scalar
    unnamed_scalar =
        MCMCChains.AxisArray(scalar; draw = range(selected), chain = chains(selected))
    @test_throws ArgumentError hcat(selected, :bad => unnamed_scalar)
    shifted = setrange(selected, 201:2:207)
    other_chain = chn[2:5, :, 1:1]
    @test_throws ArgumentError hcat(selected, :bad => shifted[:a])
    @test_throws ArgumentError hcat(selected, :bad => other_chain[:a])
    @test_throws ArgumentError hcat(selected, :bad => shifted.value)
    matching_cube = MCMCChains.AxisArray(
        cube;
        iter = range(selected),
        component = 1:2,
        chain = chains(selected),
    )
    @test hcat(selected, :labeled_cube => matching_cube)[Symbol("labeled_cube[2]")] ==
          2 .* scalar
    labeled_generated = MCMCChains.AxisArray(
        map(x -> (generated = x,), scalar);
        iter = range(selected),
        chain = chains(selected),
    )
    @test hcat(selected, labeled_generated)[:generated] == scalar
    shifted_generated = MCMCChains.AxisArray(
        map(x -> (generated = x,), scalar);
        iter = range(shifted),
        chain = chains(selected),
    )
    @test_throws ArgumentError hcat(selected, shifted_generated)
    shifted_any_generated =
        MCMCChains.AxisArray(any_generated; iter = range(shifted), chain = chains(selected))
    @test_throws ArgumentError hcat(selected, shifted_any_generated)
    other_chain_generated = MCMCChains.AxisArray(
        map(x -> (generated = x,), scalar);
        iter = range(selected),
        chain = chains(other_chain),
    )
    @test_throws ArgumentError hcat(selected, other_chain_generated)

    partly_missing = Matrix{Union{Missing,Float64}}(scalar)
    partly_missing[1, 1] = missing
    with_missing = hcat(selected, :maybe => partly_missing)
    @test ismissing(with_missing[:maybe][1, 1])
    @test with_missing[:maybe][2:end, :] == partly_missing[2:end, :]

    @test_throws DimensionMismatch hcat(selected, :bad => rand(3, 1))
    @test_throws DimensionMismatch hcat(selected, :bad => rand(4, 2, 2))
    @test_throws DimensionMismatch hcat(selected, generated[1:3, :])
    @test_throws DimensionMismatch hcat(selected, any_generated[1:3, :])
    @test_throws ArgumentError hcat(selected, :bad => rand(4))
    @test_throws ArgumentError hcat(selected, :a => scalar)
    @test_throws ArgumentError hcat(selected, :total => scalar, :total => scalar)
    @test_throws ArgumentError hcat(selected, :bad => reshape(["a", "b", "c", "d"], 4, 1))
    inconsistent = reshape([[1], [2, 3], [4], [5]], 4, 1)
    @test_throws ArgumentError hcat(selected, :bad => inconsistent)
    @test_throws ArgumentError hcat(
        selected,
        fill((;), size(selected, 1), size(selected, 3)),
    )
end
