# The search for LDAGs (csi-search): rules, model and output. Port of the
# `csisearch` class in the R package's `csisearch.cpp`.
#
# A term `Dist(a, b, c, d)` of an LDAG is p(a | b) with the variables in `c`
# assigned 0 and those in `d` assigned 1; assigned variables are in `a ∪ b`.

"""
    CsiRule

The inference rules of the search for LDAGs. The integer value of each rule
is its rule code in the R package (which are different from the codes of
[`Rule`](@ref)), so `CsiRule(5)` is `CsiCaseByCaseZero`.

| Rule | Code | Derivation (current term ⇒ new term) |
|:-----|-----:|:-------------------------------------|
| `CsiMarginalize`         |  0 | ``P(y,z∣x) ⇒ P(y∣x)`` |
| `CsiCondition`           |  1 | ``P(y,z∣x) ⇒ P(y∣x,z)`` |
| `CsiProductConditional`  |  2 | ``P(y∣x,z)`` and ``P(z∣x) ⇒ P(y,z∣x)`` |
| `CsiProductMarginal`     | -2 | ``P(y∣x)`` and ``P(z∣x,y) ⇒ P(y,z∣x)`` |
| `CsiInsertObservation`   |  3 | ``P(y∣x) ⇒ P(y∣x,z)`` by a CSI |
| `CsiDeleteObservation`   | -3 | ``P(y∣x,z) ⇒ P(y∣x)`` by a CSI |
| `CsiInsertZero`          |  4 | ``P(y∣x) ⇒ P(y∣x,z=0)`` by a CSI |
| `CsiInsertOne`           | -4 | ``P(y∣x) ⇒ P(y∣x,z=1)`` by a CSI |
| `CsiCaseByCaseZero`      |  5 | ``P(y,z=0∣x)`` and ``P(y,z=1∣x) ⇒ P(y,z∣x)`` |
| `CsiCaseByCaseOne`       | -5 | ``P(y,z=1∣x)`` and ``P(y,z=0∣x) ⇒ P(y,z∣x)`` |
| `CsiGeneralByCase`       |  6 | ``P(y∣x)`` and ``P(y,z=v∣x) ⇒ P(y,z=1-v∣x)``, from the term with `z` assigned |
| `CsiGeneralByCaseZero`   |  7 | ``P(y∣x)`` and ``P(y,z=0∣x) ⇒ P(y,z=1∣x)`` |
| `CsiGeneralByCaseOne`    | -7 | ``P(y∣x)`` and ``P(y,z=1∣x) ⇒ P(y,z=0∣x)`` |
| `CsiCaseByGeneralZero`   |  8 | ``P(y,z∣x) ⇒ P(y,z=0∣x)`` |
| `CsiCaseByGeneralOne`    | -8 | ``P(y,z∣x) ⇒ P(y,z=1∣x)`` |

`CsiInputDistribution` marks the input distributions and is not a rule.
"""
@enum CsiRule::Int8 begin
    CsiInputDistribution = 99
    CsiMarginalize = 0
    CsiCondition = 1
    CsiProductConditional = 2
    CsiProductMarginal = -2
    CsiInsertObservation = 3
    CsiDeleteObservation = -3
    CsiInsertZero = 4
    CsiInsertOne = -4
    CsiCaseByCaseZero = 5
    CsiCaseByCaseOne = -5
    CsiGeneralByCase = 6
    CsiGeneralByCaseZero = 7
    CsiGeneralByCaseOne = -7
    CsiCaseByGeneralZero = 8
    CsiCaseByGeneralOne = -8
end

"""
    default_csi_rules() -> Vector{CsiRule}

The rules the search for LDAGs tries, in order, unless a rule list is given.
"""
default_csi_rules() = CsiRule.([0, 1, -3, 3, -4, 4, -5, 5, 6, -7, 7, -8, 8, -2, 2])

function rule_label(r::CsiRule)
    code = Int(r)
    code == 0 && return "M"
    code == 1 && return "C"
    abs(code) == 2 && return "P"
    code == 3 && return "I+"
    code == -3 && return "I-"
    code == 4 && return "I+0"
    code == -4 && return "I+1"
    abs(code) == 5 && return "CbC"
    (abs(code) == 6 || abs(code) == 7) && return "GbC"
    abs(code) == 8 && return "CbG"
    return ""
end

"""
    CsiSeparation(x, y, z, zero, one)

The CSI statement that a rule needs: `x` is independent of `y` given `z`, in
the context of the term's assignments `zero` and `one` (before restricting
them to context and intervention variables).
"""
struct CsiSeparation
    x::VarSet
    y::VarSet
    z::VarSet
    zero::VarSet
    one::VarSet
end

const NO_CSI = CsiSeparation(VarSet(), VarSet(), VarSet(), VarSet(), VarSet())

"""
    LdagModel

The search model of an [`LdagProblem`](@ref): the problem, its CSI graph and
the variables that appear on the left of the input distributions.
"""
struct LdagModel <: SearchModel
    problem::LdagProblem
    graph::CsiGraph
    n::Int
    con::VarSet       # context variables
    intv::VarSet      # intervention variables
    lhs::VarSet       # variables on the left of some input
end

function LdagModel(p::LdagProblem; cache::Bool = true)
    lhs = reduce(∪, (d.a for d in p.data); init = VarSet())
    return LdagModel(p, CsiGraph(p; cache), problem_nvars(p), p.context_vars, p.intervention_vars, lhs)
end

problem(m::LdagModel) = m.problem
input_rule(::LdagModel) = CsiInputDistribution
no_separation(::LdagModel) = NO_CSI
left_variables(m::LdagModel) = m.lhs

# `valid_rule`. Rule -5 is skipped when the term assigns 1 to exactly the
# first variable: the C++ code compares the set with 1, not with 0.
function rule_applies(m::LdagModel, r::CsiRule, p::Dist)
    if r == CsiMarginalize || r == CsiCondition
        return length(p.a) != 1
    elseif r == CsiProductConditional || r == CsiDeleteObservation
        return !isempty(p.b)
    elseif r == CsiProductMarginal
        return p.a ∪ p.b != m.lhs
    elseif r == CsiCaseByCaseZero
        return !isempty(p.c)
    elseif r == CsiCaseByCaseOne
        return p.d.bits != 1
    end
    return true
end

# Rules other than marginalization and conditioning only try single
# variables (`rule_limit`).
subset_count(::LdagModel, r::CsiRule, k::Int) = abs(Int(r)) <= 1 ? (1 << k) - 1 : k

function subset_mask(m::LdagModel, r::CsiRule, p::Dist)
    a, b, c, d = p.a, p.b, p.c, p.d
    others = setdiff(fullset(m.n), a ∪ b)
    if r == CsiMarginalize
        return setdiff(a, c ∪ d)
    elseif r == CsiCondition
        return a
    elseif r == CsiProductConditional
        return b ∩ m.lhs
    elseif r == CsiProductMarginal
        return others ∩ m.lhs
    elseif r == CsiInsertObservation
        return setdiff(others, m.intv)
    elseif r == CsiDeleteObservation
        return setdiff(b, a)
    elseif r == CsiInsertZero
        return others ∩ m.con
    elseif r == CsiInsertOne
        return others ∩ (m.con ∪ m.intv)
    elseif r == CsiCaseByCaseZero
        return a ∩ c ∩ m.lhs
    elseif r == CsiCaseByCaseOne
        return a ∩ d ∩ m.lhs
    elseif r == CsiGeneralByCase                 # single variables, assigned
        return a ∩ m.lhs ∩ (c ∪ d)
    elseif r == CsiGeneralByCaseZero || r == CsiGeneralByCaseOne
        return others ∩ m.con ∩ m.lhs
    elseif r == CsiCaseByGeneralZero || r == CsiCaseByGeneralOne   # single variables
        return setdiff(m.con ∩ (a ∪ b), c ∪ d)
    end
    return VarSet()
end

"""
    apply_rule(m::LdagModel, rule, p, z) -> Union{Inference{CsiSeparation},Nothing}

Apply a [`CsiRule`](@ref) to the term `p` with the subset `z` (`apply_rule`
and `get_ruleinfo` in csisearch.cpp). The C++ code adds and subtracts the
masks; on disjoint sets and subsets these are the set operations used here.
"""
function apply_rule(m::LdagModel, r::CsiRule, p::Dist, z::VarSet)
    y, x, u, v = p.a, p.b, p.c, p.d    # p(y | x) with u = 0 and v = 1
    lhs, con = m.lhs, m.con
    plain(to) = Inference(to, NO_DIST, NO_CSI, false)
    if r == CsiMarginalize
        z ⊆ y && z != y && isempty(z ∩ (u ∪ v)) || return nothing
        return plain(Dist(setdiff(y, z), x, u, v))
    elseif r == CsiCondition
        z ⊆ y && z != y && isempty(setdiff(y, z) ∩ (u ∪ v)) || return nothing
        return plain(Dist(setdiff(y, z), x ∪ z, u, v))
    elseif r == CsiProductConditional
        z ⊆ x && z ⊆ lhs || return nothing
        xmz = setdiff(x, z)
        return Inference(Dist(y ∪ z, xmz, u, v), Dist(z, xmz, u ∩ x, v ∩ x), NO_CSI, false)
    elseif r == CsiProductMarginal
        isempty(z ∩ y) && isempty(z ∩ x) && z ⊆ lhs || return nothing
        return Inference(Dist(y ∪ z, x, u, v), Dist(z, x ∪ y, u, v), NO_CSI, true)
    elseif r == CsiInsertObservation
        isempty(z ∩ y) && isempty(z ∩ x) && isempty(z ∩ m.intv) || return nothing
        return Inference(Dist(y, x ∪ z, u, v), NO_DIST, CsiSeparation(z, y, x, u ∩ x, v ∩ x), false)
    elseif r == CsiDeleteObservation
        isempty(z ∩ y) && z ⊆ x || return nothing
        xmz, umz, vmz = setdiff(x, z), setdiff(u, z), setdiff(v, z)
        return Inference(Dist(y, xmz, umz, vmz), NO_DIST, CsiSeparation(z, y, xmz, umz ∩ x, vmz ∩ x), false)
    elseif r == CsiInsertZero
        isempty(z ∩ y) && isempty(z ∩ x) && z ⊆ con || return nothing
        return Inference(Dist(y, x ∪ z, u ∪ z, v), NO_DIST, CsiSeparation(z, y, x, u ∩ x, v ∩ x), false)
    elseif r == CsiInsertOne
        isempty(z ∩ y) && isempty(z ∩ x) && z ⊆ (con ∪ m.intv) || return nothing
        return Inference(Dist(y, x ∪ z, u, v ∪ z), NO_DIST, CsiSeparation(z, y, x, u ∩ x, v ∩ x), false)
    elseif r == CsiCaseByCaseZero
        z ⊆ y && z ⊆ u && z ⊆ lhs || return nothing
        umz = setdiff(u, z)
        return Inference(Dist(y, x, umz, v), Dist(y, x, umz, v ∪ z), NO_CSI, false)
    elseif r == CsiCaseByCaseOne
        z ⊆ y && z ⊆ v && z ⊆ lhs || return nothing
        vmz = setdiff(v, z)
        return Inference(Dist(y, x, u, vmz), Dist(y, x, u ∪ z, vmz), NO_CSI, false)
    elseif r == CsiGeneralByCase
        z ⊆ y && (z ⊆ u || z ⊆ v) && z ⊆ lhs || return nothing
        umz, vmz = setdiff(u, z), setdiff(v, z)
        to = Dist(y, x, umz ∪ (v ∩ z), vmz ∪ (u ∩ z))
        return Inference(to, Dist(setdiff(y, z), x, umz, vmz), NO_CSI, false)
    elseif r == CsiGeneralByCaseZero
        isempty(z ∩ y) && isempty(z ∩ x) && z ⊆ con && z ⊆ lhs || return nothing
        return Inference(Dist(y ∪ z, x, u, v ∪ z), Dist(y ∪ z, x, u ∪ z, v), NO_CSI, false)
    elseif r == CsiGeneralByCaseOne
        isempty(z ∩ y) && isempty(z ∩ x) && z ⊆ con && z ⊆ lhs || return nothing
        return Inference(Dist(y ∪ z, x, u ∪ z, v), Dist(y ∪ z, x, u, v ∪ z), NO_CSI, false)
    elseif r == CsiCaseByGeneralZero
        (z ⊆ y || z ⊆ x) && isempty(z ∩ u) && isempty(z ∩ v) && z ⊆ con || return nothing
        return plain(Dist(y, x, u ∪ z, v))
    elseif r == CsiCaseByGeneralOne
        (z ⊆ y || z ⊆ x) && isempty(z ∩ u) && isempty(z ∩ v) && z ⊆ con || return nothing
        return plain(Dist(y, x, u, v ∪ z))
    end
    return nothing
end

# `separation_criterion`: the assignments restricted to context and
# intervention variables.
function separation_holds(m::LdagModel, s::CsiSeparation)
    u = s.zero ∩ m.con
    v = s.one ∩ m.con
    iv = s.one ∩ m.intv
    return csi_holds!(m.graph, s.x.bits, s.y.bits, s.z.bits, u.bits, v.bits, iv.bits, (u ∪ v).bits)
end

# Candidates of the reverse product rule (`enumerate_candidates`): the
# required term, and the same term with any assignments of its context
# variables, in a stack (so in reverse order).
function candidates!(m::LdagModel, st::SearchState, inf::Inference, from::Dist)
    req = inf.required
    found = empty!(st.candidates)
    exist = get(st.index, req, 0)
    exist > 0 && push!(found, exist)
    acon = req.a ∩ m.con
    if !isempty(acon)
        f = (1 << length(acon)) - 1
        for i in 0:f, j in 0:f    # assignments to zero, to one
            i & j == 0 || continue
            zinc = deposit(VarSet(UInt64(i)), acon)
            oinc = deposit(VarSet(UInt64(j)), acon)
            rq = Dist(req.a, req.b, VarSet(req.c.bits + zinc.bits), VarSet(req.d.bits + oinc.bits))
            exist = get(st.index, rq, 0)
            exist > 0 && push!(found, exist)
        end
    end
    return reverse!(found)
end

# The derived term takes the assignments of the candidate (`assign_candidate`).
candidate_term(::LdagModel, to::Dist, req::Dist) = Dist(to.a, to.b, req.c, req.d)

function is_primitive(::LdagModel, r::CsiRule, primitive1::Bool, primitive2::Bool)
    code = abs(Int(r))
    return primitive1 && primitive2 && code != 3 && code != 4
end

# csisearch does not count the query as a new index.
target_index(::LdagModel, index::Int) = index - 1

# `csisearch_heuristic::compute_score`; input distributions are scored too.
function heuristic_score(m::LdagModel, p::Dist)
    t = m.problem.query
    ca, cb, cc, cd = p.a ∩ t.a, p.b ∩ t.b, p.c ∩ t.c, p.d ∩ t.d
    return 10 * length(ca) + 5 * length(cb) + 3 * length(cc) + 3 * length(cd) -
           2 * length(setdiff(t.a, ca)) - 2 * length(setdiff(p.a, ca)) -
           2 * length(setdiff(p.b, cb)) - 2 * length(setdiff(t.b, cb)) -
           length(setdiff(t.c, cc)) - length(setdiff(p.c, cc)) -
           length(setdiff(t.d, cd)) - length(setdiff(p.d, cd))
end
input_score(m::LdagModel, p::Dist) = heuristic_score(m, p)

# Output -----------------------------------------------------------------------

# Labels of `s`, with the assignments in `zero` and `one` (`dec_to_text`).
function set_string(io::IO, p::LdagProblem, s::VarSet, zero::VarSet, one::VarSet)
    first = true
    for i in s
        first || print(io, ",")
        first = false
        print(io, p.labels[i])
        i in zero ? print(io, " = 0") : i in one && print(io, " = 1")
    end
    return nothing
end

function set_string(p::LdagProblem, s::VarSet)
    io = IOBuffer()
    set_string(io, p, s, VarSet(), VarSet())
    return String(take!(io))
end

"""
    term_string(problem::LdagProblem, p::Dist) -> String

The term `p` in the notation of the R package, e.g. `p(Y|X,I_X = 1)`.
"""
function term_string(p::LdagProblem, t::Dist)
    io = IOBuffer()
    print(io, "p(")
    set_string(io, p, t.a, t.a ∩ t.c, t.a ∩ t.d)
    if !isempty(t.b)
        print(io, "|")
        set_string(io, p, t.b, t.b ∩ t.c, t.b ∩ t.d)
    end
    print(io, ")")
    return String(take!(io))
end
