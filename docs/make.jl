using Documenter
using MCMCChains

DocMeta.setdocmeta!(MCMCChains, :DocTestSetup, :(using MCMCChains); recursive = true)

makedocs(
    sitename = "MCMCChains.jl",
    pages = [
        "MCMCChains" => "index.md",
        "Getting started" => "getting-started.md",
        "Plotting" => ["StatsPlots.jl" => "statsplots.md", "Makie.jl" => "makie.md"],
        "Serialization" => "serialization.md",
        "API" => [
            "Chains" => "chains.md",
            "Diagnostics" => "diagnostics.md",
            "Posterior statistics" => "stats.md",
            "Model selection" => "modelstats.md",
            "Summarize" => "summarize.md",
        ],
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [MCMCChains],
    checkdocs = :exports,
)
