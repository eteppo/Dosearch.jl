"""
    Dosearch

Identification of causal effects from multiple, possibly incomplete,
observational and experimental data sources by a search over do-calculus and
probability manipulations. A Julia port of the R package dosearch (Tikka,
Hyttinen and Karvanen, 2021, *Journal of Statistical Software* 99(5),
doi:10.18637/jss.v099.i05).

The main entry point is [`dosearch`](@ref).
"""
module Dosearch

# `latex` extends the function of the standard library, which renders Markdown
# as LaTeX, so that both can be loaded together.
import Markdown: latex

# The names of normal usage are exported. The rest of the API is public but not
# exported: `Dosearch.CausalGraph`, or `using Dosearch: CausalGraph`.
export dosearch, DosearchProblem, latex, typst
export maximal_identifiable_graphs, maximal_identifiable_graph
export BidirectedPair, LatentCause, confounding_candidates, with_candidates

const PUBLIC_NAMES = (:DosearchResult, :DagProblem, :LdagProblem, :isidentifiable, :Derivation,
                      :Rule, :default_rules, :CsiRule, :default_csi_rules,
                      :SearchStrategy, :FoundOrder, :ScoreHeuristic,
                      :Candidate, :MaximalGraphs, :CausalGraph, :LabeledEdge, :P)
VERSION >= v"1.11.0-DEV.469" && eval(Expr(:public, PUBLIC_NAMES...))

include("sets.jl")
include("dist.jl")
include("rules.jl")
include("graph.jl")
include("parse.jl")
include("problem.jl")
include("ldag_problem.jl")
include("csi.jl")
include("frontier.jl")
include("search.jl")
include("dagsearch.jl")
include("csisearch.jl")
include("output.jl")
include("formula.jl")
include("derivation.jl")
include("interface.jl")
include("candidates.jl")
include("inputs.jl")
include("maximal.jl")

end
