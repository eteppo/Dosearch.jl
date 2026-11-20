# Build the documentation:
#
#     julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
#
# The HTML is written to docs/build. Doctests run as part of the build.
using Documenter, Dosearch

DocMeta.setdocmeta!(Dosearch, :DocTestSetup, :(using Dosearch); recursive = true)

makedocs(;
    sitename = "Dosearch.jl",
    modules = [Dosearch],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true", edit_link = nothing,
                             repolink = "https://github.com/eteppo/Dosearch.jl"),
    pages = [
        "Home" => "index.md",
        "Identifying causal effects" => "guide.md",
        "Tutorial: Bayesian g-formula" => "tutorial.md",
        "LDAGs" => "ldag.md",
        "Step-by-step derivations" => "derivations.md",
        "Maximal identifiable graphs" => "maximal.md",
        "Building problems in code" => "values.md",
        "Differences from R" => "differences.md",
        "API reference" => "api.md",
    ],
    checkdocs = :public,
    doctest = get(ENV, "DOCTEST_FIX", nothing) == "true" ? :fix : true,
    remotes = nothing,
)
