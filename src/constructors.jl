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

# Array

Array constructor from a Chains object. Returns 3 dimensionsal
array or an Array of 2 dimensional Arrays. If only a single parameter is selected for
inclusion, a dimension is dropped in both cases, as is e.g. required by cde(), etc.

### Method
```julia
  Array(
    chn::Chains,
    sections::Vector{Symbol};
    append_chains::Bool,
    remove_missing_union::Bool
  )
```

### Required arguments
```julia
* `chn` : Chains object to convert to an Array
```

### Optional arguments
```julia
* `sections = Symbol[]` : Sections from the Chains object to be included
* `append_chains = true`  : Append chains into a single column
* `remove_missing_union = true`  : Convert Union{Missing, Real} to Float64
```

### Examples
```julia
* `Array(chns)` : Array with chain values are appended
* `Array(chns[:par])` : Array with selected parameter chain values are appended
* `Array(chns, [:parameters])`  : Array with only :parameter section
* `Array(chns, [:parameters, :internals])`  : Array includes multiple sections
* `Array(chns, append_chains=false)`  : Array of Arrays, each chain in its own array
* `Array(chns, remove_missing_union=false)` : No conversion to remove missing values
```

"""
function Base.Array(chain::Chains,
    sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters];
    append_chains = true,
    remove_missing_union = true,
    showall = false
)
    sections = showall ? keys(chain.name_map) : sections
    if remove_missing_union
        chn = concretize(Chains(chain, sections))
    else
        chn = Chains(chain, sections)
    end

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
    return Matrix(reshape(permutedims(chain.value.data, (1, 3, 2)), :, size(chain, 2)))
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

