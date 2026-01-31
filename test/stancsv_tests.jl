using MCMCChains
using CSV
using Test

@testset "CSV Serialization" begin

    @testset "Simple CSV via Tables.jl" begin
        chn = Chains(randn(50, 3, 2), [:a, :b, :c])

        mktempdir() do tmpdir
            f = joinpath(tmpdir, "chain.csv")
            CSV.write(f, chn)
            chn2 = Chains(CSV.File(f))
            @test :a in names(chn2)
            @test :b in names(chn2)
            @test :c in names(chn2)
        end
    end

    @testset "Parameter Name Conversion" begin
        @testset "to_stan: theta[1] → theta.1" begin
            chn = Chains(randn(10, 3, 1), [Symbol("x[1]"), Symbol("y[1,2]"), :z])
            mktempdir() do tmpdir
                f = joinpath(tmpdir, "out.csv")
                write_stancsv(f, chn; include_adaptation = false, include_timing = false)
                header = readline(f)
                @test contains(header, "x.1")
                @test contains(header, "y.1.2")
                @test contains(header, "z")
                @test !contains(header, "[")
            end
        end

        @testset "from_stan: theta.1 → theta[1]" begin
            stan_csv = "lp__,theta.1,beta.2.3\n-1.0,0.5,1.5\n-2.0,0.6,1.6"
            mktempdir() do tmpdir
                f = joinpath(tmpdir, "stan.csv")
                write(f, stan_csv)
                chn = read_stancsv(f)
                @test Symbol("theta[1]") in names(chn)
                @test Symbol("beta[2,3]") in names(chn)
            end
        end
    end

    @testset "Column Ordering" begin
        chn = Chains(
            randn(10, 4, 1),
            [:mu, :sigma, :lp__, :accept_stat__],
            Dict(:internals => [:lp__, :accept_stat__]),
        )
        mktempdir() do tmpdir
            f = joinpath(tmpdir, "out.csv")
            write_stancsv(f, chn; include_adaptation = false, include_timing = false)
            header = readline(f)
            lp_pos = findfirst("lp__", header)
            mu_pos = findfirst("mu", header)
            @test lp_pos.start < mu_pos.start
        end
    end

    @testset "Internals Classification" begin
        stan_csv = "lp__,accept_stat__,theta,sigma\n-1.0,0.9,0.5,1.0"
        mktempdir() do tmpdir
            f = joinpath(tmpdir, "stan.csv")
            write(f, stan_csv)
            chn = read_stancsv(f)
            @test :lp__ in chn.name_map.internals
            @test :accept_stat__ in chn.name_map.internals
            @test :theta in chn.name_map.parameters
            @test :sigma in chn.name_map.parameters
        end
    end

    @testset "Round Trip Data Integrity" begin
        val = randn(20, 3, 1)
        chn = Chains(val, [:a, :b, :c])
        mktempdir() do tmpdir
            f = joinpath(tmpdir, "rt.csv")
            write_stancsv(f, chn; include_adaptation = false, include_timing = false)
            chn2 = read_stancsv(f)
            @test Array(chn[:, :a, 1]) ≈ Array(chn2[:, :a, 1])
            @test Array(chn[:, :b, 1]) ≈ Array(chn2[:, :b, 1])
        end
    end

    @testset "Multi-chain Export" begin
        chn = Chains(randn(10, 2, 3), [:a, :b])
        mktempdir() do tmpdir
            base = joinpath(tmpdir, "chain.csv")
            files = write_stancsv(
                base,
                chn,
                Val(:all);
                include_adaptation = false,
                include_timing = false,
            )
            @test length(files) == 3
            for f in files
                @test isfile(f)
            end
        end
    end

    @testset "Multi-chain Read" begin
        chn = Chains(randn(10, 2, 2), [:a, :b])
        mktempdir() do tmpdir
            f1 = joinpath(tmpdir, "c1.csv")
            f2 = joinpath(tmpdir, "c2.csv")
            write_stancsv(
                f1,
                chn;
                chain_id = 1,
                include_adaptation = false,
                include_timing = false,
            )
            write_stancsv(
                f2,
                chn;
                chain_id = 2,
                include_adaptation = false,
                include_timing = false,
            )
            chn2 = read_stancsv([f1, f2])
            @test size(chn2, 3) == 2
        end
    end

    @testset "Read Real CmdStan Output" begin
        stan_csv = """# stan_version_major = 2
# stan_version_minor = 35
# model = bernoulli_model
# method = sample (Default)
#   sample
#     num_samples = 1000 (Default)
#     num_warmup = 1000 (Default)
# id = 1 (Default)
# random
#   seed = 12345 (Default)
# Adaptation terminated
# Step size = 0.932037
# Diagonal elements of inverse mass matrix:
# 0.591014
lp__,accept_stat__,stepsize__,treedepth__,n_leapfrog__,divergent__,energy__,theta
-7.81355,0.761218,0.932037,1,3,0,8.13669,0.102158
-8.71049,0.894869,0.932037,1,1,0,8.72268,0.0676546
# 
#  Elapsed Time: 0.005 seconds (Warm-up)
#                0.012 seconds (Sampling)
#                0.017 seconds (Total)
# 
"""
        mktempdir() do tmpdir
            f = joinpath(tmpdir, "stan.csv")
            write(f, stan_csv)
            chn = read_stancsv(f)
            @test size(chn, 1) == 2
            @test :lp__ in names(chn)
            @test :theta in names(chn)
            @test chn.info.model_name == "bernoulli_model"
            @test chn.info.seed == 12345
            @test chn.info.num_warmup == 1000
            @test chn.info.stepsize ≈ 0.932037
            @test Array(chn[:, :theta, 1])[:] ≈ [0.102158, 0.0676546]
        end
    end

    @testset "Adaptation Comments" begin
        chn = Chains(randn(5, 2, 1), [:a, :b])
        mktempdir() do tmpdir
            f = joinpath(tmpdir, "out.csv")
            write_stancsv(f, chn; include_timing = false)
            content = read(f, String)
            @test contains(content, "# Adaptation terminated")
            @test contains(content, "# Diagonal elements")
        end
    end

    @testset "Timing Comments" begin
        chn = Chains(randn(5, 2, 1), [:a, :b])
        mktempdir() do tmpdir
            f = joinpath(tmpdir, "out.csv")
            write_stancsv(f, chn; include_adaptation = false)
            content = read(f, String)
            @test contains(content, "Elapsed Time")
        end
    end

    @testset "No Comments Mode" begin
        chn = Chains(randn(5, 2, 1), [:a, :b])
        mktempdir() do tmpdir
            f = joinpath(tmpdir, "out.csv")
            write_stancsv(f, chn; include_adaptation = false, include_timing = false)
            content = read(f, String)
            @test !startswith(content, "#")
            @test startswith(content, "a") || startswith(content, "b")
        end
    end
end
