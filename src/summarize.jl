"""
    summarize(
        chains[, stats_funs...];
        append_chains=true,
        name="SummaryStats",
        [sections, var_names],
    )

Summarize `chains` in a `PosteriorStats.SummaryStats`.

`stats_funs` is a collection of functions that reduces a matrix with shape `(draws, chains)`
to a scalar or a collection of scalars. Alternatively, an item in `stats_funs` may be a
`Pair` of the form `name => fun` specifying the name to be used for the statistic or of the
form `(name1, ...) => fun` when the function returns a collection. When the function returns
a collection, the names in this latter format must be provided.

If no stats functions are provided, then those specified in [`default_summary_stats`](@ref)
are computed.

`var_names` specifies the names of the parameters in data. If not provided, the names are
inferred from data.

# Examples

* `summarize(chns)` : Complete chain summary
* `summarize(chns[[:parm1, :parm2]])` : Chain summary of selected parameters
* `summarize(chns; sections=[:parameters])`  : Chain summary of :parameters section
* `summarize(chns; sections=[:parameters, :internals])` : Chain summary for multiple sections
"""
function PosteriorStats.summarize(
    chains::Chains, funs...;
    sections = _default_sections(chains),
    append_chains::Bool = true,
    var_names=nothing,
    kwargs...
)
    # Generate a chain to work on.
    chn = Chains(chains, _clean_sections(chains, sections))

    # Obtain names of parameters.
    names_of_params = var_names === nothing ? names(chn) : var_names

    if append_chains
        # Evaluate the functions.
        data = _permutedims_diagnostics(chains.value.data)
        summarize(data, funs...; var_names=names_of_params, kwargs...)
    else
        # Evaluate the functions.
        data = to_vector_of_matrices(chn)
        return map(data) do x
            z = reshape(x, size(x, 1), 1, size(x, 2))
            summarize(z, funs...; var_names=names_of_params, kwargs...)
        end
    end
end
