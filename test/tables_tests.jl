using Tables
using TableTraits
using IteratorInterfaceExtensions
using DataFrames

@testset "Tables interface tests" begin

    @testset "Chains" begin
        val = rand(1000, 8, 4)
        colnames = ["a", "b", "c", "d", "e", "f", "g", "h"]
        internal_colnames = ["c", "d", "e", "f", "g", "h"]
        chn = Chains(val, colnames, Dict(:internals => internal_colnames))

        @testset "Tables interface" begin
            @test Tables.istable(typeof(chn))

            @testset "column access" begin
                @test Tables.columnaccess(typeof(chn))
                @test Tables.columns(chn) === chn
                @test Tables.columnnames(chn) ==
                      (:iteration, :chain, :a, :b, :c, :d, :e, :f, :g, :h)
                @test Tables.getcolumn(chn, :iteration) == [1:1000; 1:1000; 1:1000; 1:1000]
                @test Tables.getcolumn(chn, :chain) ==
                      [fill(1, 1000); fill(2, 1000); fill(3, 1000); fill(4, 1000)]
                @test Tables.getcolumn(chn, :a) == [
                    vec(chn[:, :a, 1])
                    vec(chn[:, :a, 2])
                    vec(chn[:, :a, 3])
                    vec(chn[:, :a, 4])
                ]
                @test_throws Exception Tables.getcolumn(chn, :j)
                @test Tables.getcolumn(chn, 1) == Tables.getcolumn(chn, :iteration)
                @test Tables.getcolumn(chn, 2) == Tables.getcolumn(chn, :chain)
                @test Tables.getcolumn(chn, 3) == Tables.getcolumn(chn, :a)
                @test_throws Exception Tables.getcolumn(chn, :i)
                @test_throws Exception Tables.getcolumn(chn, 11)
            end

            @testset "row access" begin
                @test Tables.rowaccess(typeof(chn))
                @test Tables.rows(chn) isa Tables.RowIterator
                @test eltype(Tables.rows(chn)) <: Tables.AbstractRow
                rows = collect(Tables.rows(chn))
                @test eltype(rows) <: Tables.AbstractRow
                @test size(rows) === (4000,)
                for chainid = 1:4, iterid = 1:1000
                    row = rows[(chainid-1)*1000+iterid]
                    @test Tables.columnnames(row) ==
                          (:iteration, :chain, :a, :b, :c, :d, :e, :f, :g, :h)
                    @test Tables.getcolumn(row, 1) == iterid
                    @test Tables.getcolumn(row, 2) == chainid
                    @test Tables.getcolumn(row, 3) == chn[iterid, :a, chainid]
                    @test Tables.getcolumn(row, 10) == chn[iterid, :h, chainid]
                    @test Tables.getcolumn(row, :iteration) == iterid
                    @test Tables.getcolumn(row, :chain) == chainid
                    @test Tables.getcolumn(row, :a) == chn[iterid, :a, chainid]
                    @test Tables.getcolumn(row, :h) == chn[iterid, :h, chainid]
                end
            end

            @testset "integration tests" begin
                @test length(Tables.rowtable(chn)) == 4000
                nt = Tables.rowtable(chn)[1]
                @test nt == (;
                    (k => Tables.getcolumn(chn, k)[1] for k in Tables.columnnames(chn))...
                )
                @test nt == collect(Iterators.take(Tables.namedtupleiterator(chn), 1))[1]
                nt = Tables.rowtable(chn)[2]
                @test nt == (;
                    (k => Tables.getcolumn(chn, k)[2] for k in Tables.columnnames(chn))...
                )
                @test nt == collect(Iterators.take(Tables.namedtupleiterator(chn), 2))[2]
                @test Tables.matrix(chn[:, :, 1])[:, 3:end] ≈ chn[:, :, 1].value
                @test Tables.matrix(chn[:, :, 2])[:, 3:end] ≈ chn[:, :, 2].value
                @test Tables.matrix(Tables.rowtable(chn)) ==
                      Tables.matrix(Tables.columntable(chn))
            end

            @testset "schema" begin
                @test Tables.schema(chn) isa Tables.Schema
                @test Tables.schema(chn).names ===
                      (:iteration, :chain, :a, :b, :c, :d, :e, :f, :g, :h)
                @test Tables.schema(chn).types === (
                    Int,
                    Int,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                )
            end

            @testset "exceptions raised if reserved colname used" begin
                val2 = rand(1000, 2, 4)
                chn2 = Chains(val2, ["iteration", "a"])
                @test_throws Exception Tables.columns(chn2)
                @test_throws Exception Tables.rows(chn2)
                @test_throws Exception Tables.schema(chn2)
                chn3 = Chains(val2, ["chain", "a"])
                @test_throws Exception Tables.columns(chn3)
                @test_throws Exception Tables.rows(chn3)
                @test_throws Exception Tables.schema(chn3)
            end
        end

        @testset "TableTraits interface" begin
            @test IteratorInterfaceExtensions.isiterable(chn)
            @test TableTraits.isiterabletable(chn)
            nt = collect(Iterators.take(IteratorInterfaceExtensions.getiterator(chn), 1))[1]
            @test nt ==
                  (; (k => Tables.getcolumn(chn, k)[1] for k in Tables.columnnames(chn))...)
            nt = collect(Iterators.take(IteratorInterfaceExtensions.getiterator(chn), 2))[2]
            @test nt ==
                  (; (k => Tables.getcolumn(chn, k)[2] for k in Tables.columnnames(chn))...)

            val2 = rand(1000, 2, 4)
            chn2 = Chains(val2, ["iteration", "a"])
            @test_throws Exception IteratorInterfaceExtensions.getiterator(chn2)
            chn3 = Chains(val2, ["chain", "a"])
            @test_throws Exception IteratorInterfaceExtensions.getiterator(chn3)
        end

        @testset "DataFrames.DataFrame constructor" begin
            @inferred DataFrame(chn)
            @test DataFrame(chn) isa DataFrame
            df = DataFrame(chn)
            @test isequal(Tables.columntable(df), Tables.columntable(chn))
        end
    end
end
