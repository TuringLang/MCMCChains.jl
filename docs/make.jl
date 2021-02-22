push!(LOAD_PATH,"../src/")

using Documenter
using MCMCChains

name = "MCMCChains.jl"

makedocs(
    sitename = name,
    pages = [
        "MCMCChains" => "index.md",
        "Diagnostics" => "diagnostics.md"
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(repo = "github.com/rikhuijzer/$name.git", devbranch = "master")
