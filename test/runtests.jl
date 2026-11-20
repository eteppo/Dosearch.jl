using Dosearch
using Dosearch: DosearchResult, DagProblem, LdagProblem, isidentifiable, Derivation, Rule, default_rules,
                CsiRule, default_csi_rules, SearchStrategy, FoundOrder, ScoreHeuristic, Candidate,
                MaximalGraphs, CausalGraph, LabeledEdge, P, VarSet, Dist
using Test
if VERSION >= v"1.12"
    using JET
end

@testset "Dosearch.jl" begin
    @testset "API" begin
        # only the names of normal usage are exported
        exported = filter(n -> Base.isexported(Dosearch, n) && n != :Dosearch, names(Dosearch; all = true))
        @test Set(exported) == Set([:dosearch, :DosearchProblem, :latex, :typst, :maximal_identifiable_graphs,
                                    :maximal_identifiable_graph, :BidirectedPair, :LatentCause,
                                    :confounding_candidates, :with_candidates])
        if isdefined(Base, :ispublic)
            @test all(n -> Base.ispublic(Dosearch, n) && !Base.isexported(Dosearch, n), Dosearch.PUBLIC_NAMES)
        end
    end
    @testset "core" include("test_core.jl")
    @testset "dosearch" include("test_dosearch.jl")
    @testset "errors" include("test_errors.jl")
    @testset "fixtures" include("test_fixtures.jl")
    @testset "maximal graphs" include("test_maximal.jl")
    @testset "LDAGs" include("test_ldag.jl")
    @testset "derivations" include("test_derivation.jl")
    @testset "graphs and distributions as values" include("test_inputs.jl")
    # JET's analysis of the package finds no possible errors. Its reports
    # depend on the Julia version, so it runs on recent versions only.
    if VERSION >= v"1.12"
        @testset "JET" begin
            JET.test_package(Dosearch; target_modules = (Dosearch,), toplevel_logger = nothing)
        end
    end
end
