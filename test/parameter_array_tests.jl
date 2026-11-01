using MCMCChains, Test, Random

@testset "Array-valued parameter extraction" begin
    expected = [1000c + 100j + 10i + t for t = 1:4, i = 1:3, j = 1:2, c = 1:3]
    labels = [Symbol("x[$i,$j]") for j = 1:2 for i = 1:3]
    order = [6, 2, 4, 1, 5, 3]
    chain =
        Chains(reshape(expected, 4, 6, 3)[:, order, :], labels[order]; start = 11, thin = 2)
    @test isequal(@inferred(Array(chain, :x, (3, 2))), expected)
    @test isequal(Array(chain, "x", (3, 2)), expected)
    @test eltype(Array(chain, :x, (3, 2))) === Union{Missing,Int}
    @test size(Array(chain, :x, (3, 2))) == (4, 3, 2, 3)
    @test isequal(Array(chain, :x, (Int32(3), Int16(2))), expected)
    @test isequal(Array(chain[2:3, :, [3, 1]], :x, (3, 2)), expected[2:3, :, :, [3, 1]])
    @test size(Array(chain[1:1, :, 1:1], :x, (3, 2))) == (1, 3, 2, 1)
    @test isequal(
        Array(set_section(chain, Dict(:internals => labels)), :x, (3, 2)),
        expected,
    )
    result = Array(chain, :x, (3, 2))
    result[1, 1, 1, 1] = -1
    @test chain[1, Symbol("x[1,1]"), 1] == expected[1, 1, 1, 1]

    @testset "missing components" begin
        partial = chain[labels[[1, 3, 6]]]
        result = Array(partial, :x, (3, 2))
        for j = 1:2, i = 1:3
            target =
                Symbol("x[$i,$j]") in names(partial) ? expected[:, i, j, :] :
                fill(missing, 4, 3)
            @test isequal(result[:, i, j, :], target)
        end
        @test all(ismissing, Array(chain, :absent, (2,)))
        @test size(Array(chain[Symbol[]], :x, (3, 2))) == (4, 3, 2, 3)
        @test all(ismissing, Array(chain[Symbol[]], :x, (3, 2)))
        @test size(Array(chain, :absent, (0, 2))) == (4, 0, 2, 3)
        data = reshape(Union{Missing,Float32}[1, missing, 3, 4], 2, 1, 2)
        c = Chains(data, ["x[1]"])
        @test isequal(Array(c, :x, (1,)), data)
        @test eltype(Array(c, :x, (1,))) === Union{Missing,Float32}
        for dims in ((0, 1, 2), (3, 1, 0))
            c = Chains(zeros(dims), ["x[1]"])
            @test size(Array(c, :x, (1,))) == dims
        end
    end

    @testset "names and dimensions" begin
        c = Chains(reshape(1:12, 3, 2, 2), ["x", "x_other[1]"])
        @test Array(c, :x, ()) == c[:x]
        @test all(ismissing, Array(c, :absent, ()))
        @test_throws DimensionMismatch Array(c, :x, (1,))
        for label in ("x[1, 2]", "x[ 1 , 2 ]")
            c = Chains(reshape(1:6, 3, 1, 2), [label])
            @test Array(c, :x, (1, 2))[:, 1, 2, :] == c[Symbol(label)]
        end
        for base in ("α", "a.b", "x[1]", "x+", "a\\E")
            c = Chains(reshape(1:6, 3, 1, 2), ["$base[1]"])
            @test Array(c, base, (1,))[:, 1, :] == c[Symbol("$base[1]")]
        end
        for label in ("x[0]", "x[-1]", "x[4]", "x[1,2]")
            c = Chains(reshape(1:6, 3, 1, 2), [label])
            @test_throws DimensionMismatch Array(c, :x, (3,))
        end
        for label in ("x[a]", "x[1", "x[1][2]", "x[]")
            c = Chains(reshape(1:6, 3, 1, 2), [label])
            @test_throws ArgumentError Array(c, :x, (3,))
        end
        c = Chains(reshape(1:12, 3, 2, 2), ["x[1]", "x[01]"])
        @test_throws ArgumentError Array(c, :x, (1,))
        c = Chains(reshape(1:12, 3, 2, 2), ["x", "x[1]"])
        @test_throws DimensionMismatch Array(c, :x, ())
        @test_throws DimensionMismatch Array(c, :x, (1,))
        @test_throws ArgumentError Array(chain, :x, (-1, 2))
        @test_throws DimensionMismatch Array(chain, :x, (2, 2))
        @test_throws DimensionMismatch Array(chain, :x, ())
        @test_throws MethodError Array(chain, :x, (3.0, 2))
    end

    @testset "higher ranks" begin
        rng = MersenneTwister(31)
        for dims in ((12,), (2, 3), (2, 1, 3))
            values = reshape(1:(5*prod(dims)*2), 5, dims..., 2)
            labels =
                [Symbol("p[", join(Tuple(i), ','), "]") for i in CartesianIndices(dims)]
            order = randperm(rng, length(labels))
            c = Chains(reshape(values, 5, :, 2)[:, order, :], vec(labels)[order])
            @test isequal(Array(c, :p, dims), values)
        end
    end
end
