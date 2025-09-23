"""
    summarize(
        chains[, stats_funs...];
        append_chains=true,
        [sections, var_names],
        kwargs...,
    )

Summarize `chains` in a [`PosteriorStats.SummaryStats`](@extref).

`stats_funs` is a collection of functions that reduces a matrix with shape `(draws, chains)`
to a scalar or a collection of scalars. Alternatively, an item in `stats_funs` may be a
`Pair` of the form `name => fun` specifying the name to be used for the statistic or of the
form `(name1, ...) => fun` when the function returns a collection. When the function returns
a collection, the names in this latter format must be provided.

# Keywords

- `section`: The sections of the chain to include in the summary. If not provided, defaults
  to `:parameters`.
- `append_chains`: If `true`, a single `SummaryStats` for all chains is returned. If
  `false`, a vector of `SummaryStats` (one for each chain) is returned.
- `var_names`: The names of the parameters in data. If not provided, the names are taken
  from `chains`.
- `kwargs...`: Additional keyword arguments are forwarded to
  [`PosteriorStats.summarize`](@extref).

# Examples

* `summarize(chns)` : Complete chain summary
* `summarize(chns[[:parm1, :parm2]])` : Chain summary of selected parameters
* `summarize(chns; sections=[:parameters])`  : Chain summary of :parameters section
* `summarize(chns; sections=[:parameters, :internals])` : Chain summary for multiple sections
"""
function PosteriorStats.summarize(
    chains::Chains,
    funs...;
    sections = _default_sections(chains),
    append_chains::Bool = true,
    var_names = nothing,
    name::AbstractString = "SummaryStats",
    kwargs...,
)
    # Generate a chain to work on.
    chn = Chains(chains, _clean_sections(chains, sections))

    # Obtain names of parameters.
    names_of_params = var_names === nothing ? names(chn) : var_names

    if append_chains
        # Evaluate the functions.
        data = _permutedims_diagnostics(chn.value.data)
        summarize(data, funs...; var_names = names_of_params, name, kwargs...)
    else
        # Evaluate the functions.
        data = to_vector_of_matrices(chn)
        return map(enumerate(data)) do (i, x)
            z = reshape(x, size(x, 1), 1, size(x, 2))
            name_chain = name * " (Chain $i)"
            summarize(z, funs...; var_names = names_of_params, name = name_chain, kwargs...)
        end
    end
end
