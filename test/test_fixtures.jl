# Compare against the results of the R package on the cases in
# fixtures/dosearch_fixtures.toml (see fixtures/generate_fixtures.R): the
# reference examples and tests of the R package and seeded random problems.
# Identifiability, formulas, derivations and verbose search logs must match
# exactly. The step-by-step derivation of each identified query is checked too
# (derivation_utils.jl).

include("fixture_utils.jl")
include("derivation_utils.jl")

@testset "R fixtures" begin
    cases = load_fixtures()
    @test length(cases) > 400
    for case in cases
        if 0 in get(case, "rules", Int[])
            # Intentional difference: R accepts the non-rule code 0, which the
            # C++ code applies as a copy of the previous rule application.
            # Dosearch.jl rejects it; without it the results are the same here.
            case = merge(case, Dict("rules" => filter(!=(0), case["rules"])))
        end
        result, log = run_fixture(case)
        @testset "$(case["name"])" begin
            @test result.identifiable == case["identifiable"]
            @test something(result.formula, "") == case["formula"]
            if haskey(case, "derivation")
                if case["draw_all"]
                    # Intentional difference: with draw_all, R draws the steps
                    # in hash table order and Dosearch.jl in the order they
                    # were found, so only the sets of edges are compared.
                    @test derivation_edges(result.dot) == derivation_edges(case["derivation"])
                else
                    @test result.dot == case["derivation"]
                end
            end
            if haskey(case, "verbose")
                @test strip(log) == strip(case["verbose"])
            end
            result.identifiable && check_derivation(result)
        end
    end
end
