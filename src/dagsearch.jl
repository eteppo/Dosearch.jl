# The search for DAGs (do-search): the model that connects a `DagProblem`
# with the search core. Port of the `dosearch` class in the R package's
# `dosearch.cpp`, apart from the rules (rules.jl) and output (output.jl).

"""
    DagModel

The search model of a [`DagProblem`](@ref): the problem, the constants the
rules need and the d-separation graph.
"""
struct DagModel <: SearchModel
    problem::DagProblem
    ctx::RuleContext
    graph::DsepGraph
    fold_separations::Bool   # some rule checks its separation with the subset mask
end

DagModel(p::DagProblem, rules::Vector{Rule}) =
    DagModel(p, rule_context(p), causal_graph(p), any(separation_by_mask, rules))

problem(m::DagModel) = m.problem
input_rule(::DagModel) = InputDistribution
no_separation(::DagModel) = NO_SEPARATION

function left_variables(m::DagModel)
    lhs = VarSet()
    for p in m.problem.data
        lhs = lhs ∪ p.a
        m.ctx.md && (lhs = lhs ∪ ((p.a ∩ m.ctx.md_p) >> 2))
    end
    return lhs
end

rule_applies(m::DagModel, r::Rule, p::Dist) = rule_applies(r, p, m.ctx)
subset_mask(m::DagModel, r::Rule, p::Dist) = subset_mask(r, p, m.ctx)
apply_rule(m::DagModel, r::Rule, p::Dist, z::VarSet) = apply_rule(r, p, z, m.ctx)
subset_count(::DagModel, ::Rule, k::Int) = (1 << k) - 1   # every subset (`rule_limit`)
separation_holds(m::DagModel, s::Separation) = dseparated(m.graph, s)

# Rules whose d-separation statement is "the left-hand side of the term is
# separated from `z` (rule 1) or from the intervention nodes of `z` (rules 2
# and 3), given the right-hand side, intervening on the do-variables". Only
# `z` varies, so the subsets that fail are exactly those that meet the nodes
# d-connected to the left-hand side, which are computed once per term.
separation_by_mask(r::Rule) = r == InsertObservation || r == ObservationToAction || r == InsertAction

# The variables that `z` must avoid for the separation of rule `r`, given the
# nodes `connected` to the left-hand side of the term.
function connected_variables(r::Rule, connected::UInt64, n::Int)
    return (r == InsertObservation ? VarSet(connected) : VarSet(connected >> n)) ∩ fullset(n)
end

term_data(m::DagModel, p::Dist) =
    m.fold_separations ? dconnected(m.graph, p.a.bits, p.b.bits, p.c.bits) : zero(UInt64)

function narrow_mask(m::DagModel, r::Rule, p::Dist, mask::VarSet, connected::UInt64)
    separation_by_mask(r) || return mask, false
    return setdiff(mask, connected_variables(r, connected, m.ctx.n)), true
end

# Terms that can serve as the required term of an enumerating rule: the
# required term itself, and the same term with any extra response indicators
# enabled that the rule adds to the left-hand side (`enumerate_candidates`).
# The C++ code pushes the candidates on a stack, so they are returned in
# reverse order of discovery, in the buffer `st.candidates`.
function candidates!(m::DagModel, st::SearchState, inf::Inference, from::Dist)
    req = inf.required
    found = empty!(st.candidates)
    exist = get(st.index, req, 0)
    exist > 0 && push!(found, exist)
    asw = setdiff(req.a, from.a) ∩ m.ctx.md_s
    if !isempty(asw)
        for ranks in 0:((1 << length(asw)) - 1)   # includes the empty set again
            # The C++ code adds the masks; this equals the union whenever the
            # added indicators are not already enabled.
            inc = deposit(VarSet(UInt64(ranks)), asw)
            rq = Dist(req.a, req.b, req.c, VarSet(req.d.bits + inc.bits))
            exist = get(st.index, rq, 0)
            exist > 0 && push!(found, exist)
        end
    end
    return reverse!(found)
end

# The derived term with the response indicators enabled in the candidate
# (`assign_candidate`).
candidate_term(::DagModel, to::Dist, req::Dist) = Dist(to.a, to.b, to.c, req.d ∩ (to.a ∪ to.b ∪ to.c))

is_primitive(::DagModel, r::Rule, primitive1::Bool, primitive2::Bool) =
    primitive1 && primitive2 && !isdocalculus(r)

# dosearch counts the query as a new index (`derive_distribution`).
target_index(::DagModel, index::Int) = index

input_score(::DagModel, p::Dist) = 0
heuristic_score(m::DagModel, p::Dist) = heuristic_score(p, m.problem.query)
