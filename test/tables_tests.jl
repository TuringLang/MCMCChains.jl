using Tables

@testset "Tables interface tests" begin

    @testset "Chains" begin
        @testset "Tables interface" begin
            val = rand(1000, 8, 4)
            colnames = ["a", "b", "c", "d", "e", "f", "g", "h"]
            internal_colnames = ["c", "d", "e", "f", "g", "h"]
            chn = Chains(val, colnames, Dict(:internals => internal_colnames))
            @test Tables.istable(typeof(chn))
            @test Tables.columnaccess(typeof(chn))
            @test Tables.columns(chn) === chn
            @test Tables.columnnames(chn) ==
                  [:Iteration, :Chain, :a, :b, :c, :d, :e, :f, :g, :h]
            @test Tables.getcolumn(chn, :Iteration) == [1:1000; 1:1000; 1:1000; 1:1000]
            @test Tables.getcolumn(chn, :Chain) ==
                  [fill(1, 1000); fill(2, 1000); fill(3, 1000); fill(4, 1000)]
            @test Tables.getcolumn(chn, :a) == [
                vec(chn[:, :a, 1].value)
                vec(chn[:, :a, 2].value)
                vec(chn[:, :a, 3].value)
                vec(chn[:, :a, 4].value)
            ]
            @test_throws Exception Tables.getcolumn(chn, :j)
            @test Tables.getcolumn(chn, 1) == Tables.getcolumn(chn, :Iteration)
            @test Tables.getcolumn(chn, 2) == Tables.getcolumn(chn, :Chain)
            @test Tables.getcolumn(chn, 3) == Tables.getcolumn(chn, :a)
            @test_throws Exception Tables.getcolumn(chn, :i)
            @test_throws Exception Tables.getcolumn(chn, 11)
            @test Tables.rowaccess(typeof(chn))
            @test Tables.rows(chn) === chn
            @test length(Tables.rowtable(chn)) == 4000
            nt = Tables.rowtable(chn)[1]
            @test nt ==
                  (; (k => Tables.getcolumn(chn, k)[1] for k in Tables.columnnames(chn))...)
            @test nt == collect(Iterators.take(Tables.namedtupleiterator(chn), 1))[1]
            nt = Tables.rowtable(chn)[2]
            @test nt ==
                  (; (k => Tables.getcolumn(chn, k)[2] for k in Tables.columnnames(chn))...)
            @test nt == collect(Iterators.take(Tables.namedtupleiterator(chn), 2))[2]
            @test Tables.schema(chn) isa Tables.Schema
            @test Tables.schema(chn).names ===
                  (:Iteration, :Chain, :a, :b, :c, :d, :e, :f, :g, :h)
            @test Tables.schema(chn).types === (
                Int64,
                Int64,
                Float64,
                Float64,
                Float64,
                Float64,
                Float64,
                Float64,
                Float64,
                Float64,
            )
            @test Tables.matrix(chn[:, :, 1])[:, 3:end] ≈ chn[:, :, 1].value
            @test Tables.matrix(chn[:, :, 2])[:, 3:end] ≈ chn[:, :, 2].value
        end
    end
