#################### File I/O ####################

function readcoda(output::AbstractString, index::AbstractString)
    out = readdlm(output, Any)
    ind = readdlm(index, Any)

    firstind = ind[:, 2]
    firstiter = out[firstind, 1]
    lastind = ind[:, 3]
    lastiter = out[lastind, 1]

    thin = Int((lastiter[1] - firstiter[1]) / (lastind[1] - firstind[1]))
    window = maximum(firstiter):thin:minimum(lastiter)
    startind = firstind + (first(window) - firstiter) / step(window)
    stopind = lastind - (lastiter - last(window)) / step(window)

    names = AbstractString[ind[:, 1]...]

    value = Array{Float64}(length(window), length(names))
    for i = 1:size(value, 2)
        inds = Int(startind[i]):Int(stopind[i])
        value[:, i] = out[inds, 2]
    end

    Chains(value; iterations = window, names = names)
end

"""
    write_json(c::Chains, filepath::Union{AbstractString, Nothing}=nothing; as_string::Bool=false)

Serialize the `Chains` object to JSON.

!!! note
    This function requires the `JSON package. Run `using JSON` before calling this function.

# Examples
```julia
using JSON
write_json(chains, "output.json")
```
"""
function write_json end

function write_json(c, args...; kwargs...)
    error("write_json requires the JSON package. Please run: using JSON")
end

"""
    read_json(input::AbstractString)

Deserialize a `Chains` object from a JSON file.

!!! note
    This function requires the `JSON` package. Run `using JSON` before calling this function.

# Examples
```julia
using JSON
chains = read_json("input.json")
```
"""
function read_json end

function read_json(input)
    error("read_json requires the JSON package. Please run: using JSON")
end
