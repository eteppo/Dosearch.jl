# Helpers to run the R reference fixtures (test/fixtures/dosearch_fixtures.toml).

using TOML

const FIXTURE_FILE = joinpath(@__DIR__, "fixtures", "dosearch_fixtures.toml")

load_fixtures() = TOML.parsefile(FIXTURE_FILE)["case"]

# Run the Julia search on a fixture case, capturing the verbose log if needed.
function run_fixture(case::Dict)
    kwargs = (
        transportability = get(case, "transportability", nothing),
        selection_bias = get(case, "selection_bias", nothing),
        missing_data = get(case, "missing_data", nothing),
        md_sym = get(case, "md_sym", "1"),
        heuristic = case["heuristic"],
        improve = case["improve"],
        rules = get(case, "rules", nothing),
        draw_derivation = true,
        draw_all = case["draw_all"],
        warn = false,
    )
    if haskey(case, "verbose")
        io = IOBuffer()
        result = dosearch(case["data"], case["query"], case["graph"]; kwargs..., verbose = true, io)
        return result, String(take!(io))
    end
    return dosearch(case["data"], case["query"], case["graph"]; kwargs...), nothing
end

# Edges of a DOT derivation as (from label, to label, rule label) triples,
# independent of the node numbering.
function derivation_edges(dot::AbstractString)
    labels = Dict{String,String}()
    for m in eachmatch(r"(n\d+)\[shape=polygon,sides=4,label=\"([^\"]*)\"\]", dot)
        labels[m.captures[1]] = m.captures[2]
    end
    return Set((labels[m.captures[1]], labels[m.captures[2]], m.captures[3])
               for m in eachmatch(r"(n\d+) -> (n\d+)\[label=\"([^\"]*)\"\]", dot))
end
