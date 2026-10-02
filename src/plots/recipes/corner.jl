struct Corner
    c::Any
    parameters::Any
end

@recipe function f(corner::Corner)
    # Convert labels to string because `Symbol` is not supported generally supported.
    label --> permutedims(map(string, corner.parameters))
    compact --> true
    size --> (600, 600)
    # NOTE: Don't use the indices from `chains(chains)`.
    # See https://github.com/TuringLang/MCMCChains.jl/issues/413.
    ar = collect(
        Array(corner.c.value[:, corner.parameters, i]) for i = 1:length(chains(corner.c))
    )
    RecipesBase.recipetype(:cornerplot, vcat(ar...))
end
