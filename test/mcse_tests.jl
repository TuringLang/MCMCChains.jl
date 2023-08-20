using MCMCChains

using Random
using Statistics
using Test

mymean(x) = mean(x)

@testset "mcse" begin
    x = rand(10_000, 40, 10)
    chain = Chains(x)

    for kind in (mean, std, mymean)
        if kind !== mymean
            for autocov_method in (AutocovMethod(), BDAAutocovMethod())
                # analyze chain
                mcse_df = mcse(chain; autocov_method = autocov_method, kind = kind)

                # analyze array
                mcse_array = mcse(
                    PermutedDimsArray(x, (1, 3, 2)); autocov_method = autocov_method, kind = kind,
                )
                @test mcse_df[:mcse] == mcse_array
            end
        else
            # analyze chain
            mcse_df = mcse(chain; kind = kind)

            # analyze array
            mcse_array = mcse(PermutedDimsArray(x, (1, 3, 2)); kind = kind)
            @test mcse_df[:mcse] == mcse_array
        end
    end
end
