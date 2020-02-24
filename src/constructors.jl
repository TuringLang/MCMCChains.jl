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
function Base.Array(chn::Chains,
        sections::Union{Symbol, Vector{Symbol}}=Symbol[:parameters];
        append_chains=true,
        showall=false,
        sorted=false
    )
    sections = showall ? keys(chn.name_map) : sections
    chn = Chains(chn, sections)

    arr = if append_chains
        mapreduce(i -> chn.value.data[:,:,i], vcat, 1:size(chn, 3))
    else
        map(i -> chn.value.data[:,:,i], 1:size(chn, 3))
    end

    arr = if append_chains
        reshape(arr, (size(arr, 1), size(arr, 2)))
    else
        arr
    end

    if size(arr, 2) == 1
        return map(identity, arr[:,1])
    else
        return map(identity, arr)
    end
end


