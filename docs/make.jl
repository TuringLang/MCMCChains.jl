using Documenter
using MCMCChains

# https://github.com/JuliaDocs/Documenter.jl/issues/942
ENV["COLUMNS"] = 100

DocMeta.setdocmeta!(
    MCMCChains,
    :DocTestSetup,
    :(using MCMCChains);
    recursive=true
)

makedocs(
    sitename = "MCMCChains.jl",
    pages = [
        "MCMCChains" => "index.md",
        "Getting started" => "getting-started.md",
        "Plotting" => [
            "StatsPlots.jl" => "statsplots.md",
            "Gadfly.jl" => "gadfly.md"
        ],
        "API" => [
            "Chains" => "chains.md",
            "Diagnostics" => "diagnostics.md",
            "Posterior statistics" => "stats.md",
            "Model selection" => "modelstats.md",
            "Summarize" => "summarize.md"
        ]
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [MCMCChains],
    strict = true,
    checkdocs = :exports,
    # https://github.com/JuliaLang/julia/pull/37085#issuecomment-683356098
    doctestfilters = [
        r"{([a-zA-Z0-9]+,\s?)+[a-zA-Z0-9]+}",
        r"(Array{[a-zA-Z0-9]+,\s?1}|Vector{[a-zA-Z0-9]+})",
        r"(Array{[a-zA-Z0-9]+,\s?2}|Matrix{[a-zA-Z0-9]+})",
    ]
)

deploydocs(repo = "github.com/TuringLang/MCMCChains.jl.git", push_preview=true)
