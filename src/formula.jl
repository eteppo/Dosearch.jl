# Formulas of derivations as expression trees. Each derived term of a search is
# rewritten by the rule that derived it as an expression over the terms it was
# derived from (`rewrite`). Expanding these rewrites from the query down to
# the terms that are known gives the identifying formula; the formula in the
# notation of the R package (`derive_formula` in the C++ code) is one
# rendering of it, and the step-by-step derivation (derivation.jl) renders
# the expansion one rewrite at a time.

"""
    Formula

An expression over the terms of a derivation: a [`TermFormula`](@ref), sum,
product, fraction, difference, case distinction or evaluation.
"""
abstract type Formula end

"""
    TermFormula(node)

The term of node `node` of a derivation.
"""
struct TermFormula <: Formula
    node::Int
end

# The sum of `body` over the values of `vars`.
struct SumFormula <: Formula
    vars::VarSet
    body::Formula
end

struct ProductFormula <: Formula
    left::Formula
    right::Formula
end

struct FractionFormula <: Formula
    num::Formula
    den::Formula
end

struct DifferenceFormula <: Formula
    left::Formula
    right::Formula
end

# `case0` where `var = 0` and `case1` where `var = 1` (LDAG rules
# `CsiCaseByCase*`). The R package writes the case of the first parent first.
struct CasesFormula <: Formula
    var::VarSet
    case0::Formula
    case1::Formula
    zero_first::Bool
end

# `body` evaluated at `vars = value`. The R package writes the evaluation for
# LDAGs, but not when a DAG search enables missing data mechanisms
# (`explicit = false`).
struct EvaluationFormula <: Formula
    body::Formula
    vars::VarSet
    value::String
    explicit::Bool
end

# The rewrite of a derived `node` by its rule, over the terms of its parents
# in `nodes`. A product puts first the factor whose formula in the notation of
# the R package is shorter (`rlength(i)` is that length for node `i`); the C++
# code compares byte lengths.
function ordered_product(node::Node, rlength)
    t1, t2 = TermFormula(node.parent1), TermFormula(node.parent2)
    return rlength(node.parent1) < rlength(node.parent2) ? ProductFormula(t1, t2) : ProductFormula(t2, t1)
end

function rewrite(p::DagProblem, nodes::AbstractVector{<:Node}, node::Node{Rule}, rlength)
    code = abs(Int(node.rule))
    t1 = TermFormula(node.parent1)
    if node.parent2 > 0
        t2 = TermFormula(node.parent2)
        code == 6 && return ordered_product(node, rlength)
        code == 7 && return FractionFormula(t1, t2)
        code == 8 && return FractionFormula(t2, t1)
    else
        parent = nodes[node.parent1].dist::Dist
        (code <= 3 || code == 10) && return t1
        code == 4 && return SumFormula(setdiff(parent.a, node.dist.a), t1)
        code == 5 && return FractionFormula(t1, SumFormula(node.dist.a, t1))
        code == 9 && return EvaluationFormula(t1, setdiff(node.dist.d, parent.d), p.md_sym, false)
    end
    throw(ArgumentError("rule $(node.rule) cannot derive a term"))
end

function rewrite(::LdagProblem, nodes::AbstractVector{<:Node}, node::Node{CsiRule}, rlength)
    code = abs(Int(node.rule))
    t1 = TermFormula(node.parent1)
    dist = node.dist
    parent = nodes[node.parent1].dist::Dist
    if node.parent2 > 0
        t2 = TermFormula(node.parent2)
        code == 2 && return ordered_product(node, rlength)
        if code == 5
            var = setdiff(parent.c ∪ parent.d, dist.c ∪ dist.d)
            return Int(node.rule) > 0 ? CasesFormula(var, t1, t2, true) : CasesFormula(var, t2, t1, false)
        end
        code == 6 && return DifferenceFormula(t2, t1)
        code == 7 && return DifferenceFormula(t1, t2)
    else
        (code == 3 || code == 4) && return t1
        code == 0 && return SumFormula(setdiff(parent.a, dist.a), t1)
        code == 1 && return FractionFormula(t1, SumFormula(dist.a, t1))
        if code == 8
            zero = setdiff(dist.c, parent.c)
            isempty(zero) || return EvaluationFormula(t1, zero, "0", true)
            return EvaluationFormula(t1, setdiff(dist.d, parent.d), "1", true)
        end
    end
    throw(ArgumentError("rule $(node.rule) cannot derive a term"))
end

# A formula in the notation of the R package; `term(i)` is the formula of
# node `i`.
r_formula(p, f::TermFormula, term) = term(f.node)
r_formula(p, f::SumFormula, term) = "\\sum_{" * set_string(p, f.vars) * "}" * r_formula(p, f.body, term)
r_formula(p, f::ProductFormula, term) =
    "\\left(" * r_formula(p, f.left, term) * r_formula(p, f.right, term) * "\\right)"
r_formula(p, f::FractionFormula, term) =
    "\\frac{" * r_formula(p, f.num, term) * "}{" * r_formula(p, f.den, term) * "}"
r_formula(p, f::DifferenceFormula, term) =
    "\\left(" * r_formula(p, f.left, term) * " - " * r_formula(p, f.right, term) * "\\right)"
function r_formula(p, f::CasesFormula, term)
    first, second = f.zero_first ? (f.case0, f.case1) : (f.case1, f.case0)
    return "\\left[" * r_formula(p, first, term) * " /\\ " * r_formula(p, second, term) * "\\right]"
end
function r_formula(p, f::EvaluationFormula, term)
    body = r_formula(p, f.body, term)
    f.explicit || return body
    return "\\left[" * body * "\\right]_{" * set_string(p, f.vars) * " = " * f.value * "}"
end

# The formulas of the R package of every node of `nodes`, computed on demand:
# a node whose term is known (an input distribution or a term derived from
# them by probability manipulations only) shows the term itself, and a
# derived node expands its rewrite.
struct RFormulas{P,N}
    problem::P
    nodes::N
    memo::Dict{Int,String}
end
RFormulas(p, nodes) = RFormulas(p, nodes, Dict{Int,String}())

function (rf::RFormulas)(i::Int)
    haskey(rf.memo, i) && return rf.memo[i]
    node = rf.nodes[i]
    f = if node.primitive
        term_string(rf.problem, node.dist)
    else
        r_formula(rf.problem, rewrite(rf.problem, rf.nodes, node, j -> ncodeunits(rf(j))), rf)
    end
    rf.memo[i] = f
    return f
end
