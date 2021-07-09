function sort_sections(chn::Chains)
    smap = keys(chn.name_map)
    section_list = Vector{Symbol}(undef, length(smap))
    indx = 1
    if :parameters in smap
        section_list[1] = :parameters
        indx += 1
    end
    if :internals in smap
        section_list[end] = :internals
    end
    for par in smap
        if !(par in [:parameters, :internals])
            section_list[indx] = par
            indx += 1
        end
    end
    return section_list
end

"""
    Array(chains[, sections;
          append_chains = true, remove_missing_union = true])

Construct an `Array` from a chain.

Return 3 dimensionsal array or an Array of 2 dimensional Arrays. If only a single parameter
is selected for inclusion, a dimension is dropped in both cases, as is e.g. required by
cde(), etc.

# Examples

* `Array(chns)` : Array with chain values are appended
* `Array(chns[[:par]])` : Array with selected parameter chain values are appended
* `Array(chns, :parameters)`  : Array with only :parameter section
* `Array(chns, [:parameters, :internals])`  : Array includes multiple sections
* `Array(chns; append_chains=false)`  : Array of Arrays, each chain in its own array
* `Array(chns; remove_missing_union=false)` : No conversion to remove missing values
"""
function Base.Array(
    chains::Chains,
    sections = _default_sections(chains);
    append_chains = true,
    remove_missing_union = true
)
    _chains = Chains(chains, _clean_sections(chains, sections))
    chn = remove_missing_union ? concretize(_chains) : _chains

    nparams = size(chn, 2)
    if append_chains
        if nparams == 1
            return to_vector(chn)
        else
            return to_matrix(chn)
        end
    else
        if nparams == 1
            return to_vector_of_vectors(chn)
        else
            return to_vector_of_matrices(chn)
        end
    end
end

function to_matrix(chain::Chains)
    x = permutedims(chain.value.data, (1, 3, 2))
    return Matrix(reshape(x, size(x, 1) * size(x, 2), size(x, 3)))
end

function to_vector(chain::Chains)
    return Vector(vec(chain.value.data))
end

function to_vector_of_vectors(chain::Chains)
    data = chain.value.data
    return [Vector(vec(data[:, :, i])) for i in axes(data, 3)]
end

function to_vector_of_matrices(chain::Chains)
    data = chain.value.data
    return [Matrix(data[:, :, i]) for i in axes(data, 3)]
end
