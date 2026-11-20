# Inference rules of the search: do-calculus, probability manipulations and
# the extra rules for missing data. Port of `valid_rule`, `apply_rule` and
# `get_ruleinfo` in the R package's `dosearch.cpp`.

"""
    Rule

The inference rules of the search. The integer value of each rule is the rule
code of the R package, so `Rule(-2)` is `ActionToObservation`, and a rule list
can be given either as `Rule`s or as R codes.

| Rule | Code | Derivation (current term ⇒ new term) |
|:-----|-----:|:-------------------------------------|
| `InsertObservation`      |  1 | ``P(y∣x) ⇒ P(y∣x,z)`` by do-calculus rule 1 |
| `DeleteObservation`      | -1 | ``P(y∣x,z) ⇒ P(y∣x)`` by do-calculus rule 1 |
| `ObservationToAction`    |  2 | ``P(y∣x,z) ⇒ P(y∣x,do(z))`` by do-calculus rule 2 |
| `ActionToObservation`    | -2 | ``P(y∣x,do(z)) ⇒ P(y∣x,z)`` by do-calculus rule 2 |
| `InsertAction`           |  3 | ``P(y∣x) ⇒ P(y∣x,do(z))`` by do-calculus rule 3 |
| `DeleteAction`           | -3 | ``P(y∣x,do(z)) ⇒ P(y∣x)`` by do-calculus rule 3 |
| `Marginalize`            |  4 | ``P(y,z∣x) ⇒ P(y∣x)`` |
| `Condition`              |  5 | ``P(y,z∣x) ⇒ P(y∣x,z)`` |
| `ProductConditional`     |  6 | ``P(y∣x,z)`` and ``P(z∣x) ⇒ P(y,z∣x)`` |
| `ProductMarginal`        | -6 | ``P(y∣x)`` and ``P(z∣x,y) ⇒ P(y,z∣x)`` |
| `DivideMarginal`         |  7 | ``P(y∣x) / P(z∣x) ⇒ P(y∖z∣x,z)`` |
| `DivideConditional`      | -7 | ``P(y∣x) / P(z∣x,y∖z) ⇒ P(y∖z∣x)`` |
| `DenominatorMarginal`    |  8 | ``P(y,z∣x) / P(y∣x) ⇒ P(z∣x,y)`` |
| `DenominatorConditional` | -8 | ``P(y,z∣x) / P(y∣x,z) ⇒ P(z∣x)`` |
| `EnableMechanism`        |  9 | enable (fix) a missing data mechanism |
| `ExchangeProxy`          | 10 | replace a proxy variable with its true variable |

`InputDistribution` (code 0) marks the input distributions in a derivation and
is not a rule. Rules 7 to 10 only apply to missing data problems.
"""
@enum Rule::Int8 begin
    InputDistribution = 0
    InsertObservation = 1
    DeleteObservation = -1
    ObservationToAction = 2
    ActionToObservation = -2
    InsertAction = 3
    DeleteAction = -3
    Marginalize = 4
    Condition = 5
    ProductConditional = 6
    ProductMarginal = -6
    DivideMarginal = 7
    DivideConditional = -7
    DenominatorMarginal = 8
    DenominatorConditional = -8
    EnableMechanism = 9
    ExchangeProxy = 10
end

"""
    isdocalculus(rule) -> Bool

Whether `rule` is one of the three rules of do-calculus (in either direction).
"""
isdocalculus(r::Rule) = r != InputDistribution && abs(Int(r)) <= 3

# Label of a rule in derivation graphs.
function rule_label(r::Rule)
    code = abs(Int(r))
    code == 1 && return "R1"
    code == 2 && return "R2"
    code == 3 && return "R3"
    code == 4 && return "M"
    code == 5 && return "C"
    code == 6 && return "P"
    (code == 7 || code == 8) && return "D"
    code == 9 && return "A"
    code == 10 && return "EX"
    return ""
end

"""
    default_rules(; missing_data::Bool, improve::Bool) -> Vector{Rule}

The rules the search tries, in order, unless a rule list is given. With
`improve = true` the redundant observation insertion and deletion rules are
left out.
"""
function default_rules(; missing_data::Bool, improve::Bool)
    codes = if missing_data
        improve ? [4, 5, 9, 10, -2, -3, 2, 3, 6, -6, 7, -7, 8, -8] :
                  [4, 5, 9, 10, -1, -2, -3, 1, 2, 3, 6, -6, 7, -7, 8, -8]
    else
        improve ? [4, 5, -2, -3, 2, 3, 6, -6] :
                  [4, 5, -1, -2, -3, 1, 2, 3, 6, -6]
    end
    return Rule.(codes)
end

# Constants of a model that the rules need.
struct RuleContext
    n::Int              # number of variables
    md::Bool            # missing data problem?
    md_s::VarSet        # response indicators (missing data switches)
    md_p::VarSet        # proxy variables
    md_t::VarSet        # true variables of missing data
    trsb::VarSet        # transportability and selection bias nodes
end

# Can `rule` apply to term `p` for some subset? (`valid_rule` in dosearch.cpp)
function rule_applies(r::Rule, p::Dist, ctx::RuleContext)
    a, b, c, d = p.a, p.b, p.c, p.d
    if r == DeleteObservation || r == ObservationToAction || r == ProductConditional
        return b != c                      # there must be observations
    elseif r == ActionToObservation || r == DeleteAction
        return !isempty(c)                 # there must be actions
    elseif r == Marginalize || r == Condition
        return length(a) != 1              # there must be other variables
    elseif r == DivideMarginal || r == DivideConditional
        return length(a) != 1 && !isempty(a ∩ d)
    elseif r == DenominatorMarginal || r == DenominatorConditional
        return a ∩ d == a                  # only enabled mechanisms on the left
    elseif r == EnableMechanism
        return !isempty((a ∪ b) ∩ ctx.md_s)
    elseif r == ExchangeProxy
        return !isempty(d) && !isempty((a ∪ b) ∩ ctx.md_p)
    end
    return true
end

"""
    subset_mask(rule, p, ctx) -> VarSet

A set that contains every subset `z` for which `apply_rule(rule, p, z, ctx)`
can succeed. The search only tries subsets of this mask; since the subsets of
a set in banker's order are the subsequence of the full banker's sequence
that lies inside the set, the order of the search does not change.
"""
function subset_mask(r::Rule, p::Dist, ctx::RuleContext)
    y, xw, x, d = p.a, p.b, p.c, p.d
    others = setdiff(fullset(ctx.n), y ∪ xw)
    if r == InsertObservation || r == InsertAction
        return others
    elseif r == DeleteObservation
        return setdiff(xw, x)
    elseif r == ObservationToAction
        return setdiff(setdiff(xw, x), ctx.md_p)
    elseif r == ActionToObservation || r == DeleteAction
        return setdiff(x, ctx.md_p)
    elseif r == Marginalize || r == Condition
        return setdiff(y, ctx.trsb)
    elseif r == ProductConditional
        return setdiff(setdiff(xw, x), ctx.trsb)
    elseif r == ProductMarginal || r == DenominatorMarginal
        return setdiff(others, ctx.trsb)
    elseif r == DivideMarginal || r == DivideConditional
        return setdiff(y ∩ d, ctx.trsb)
    elseif r == DenominatorConditional
        return setdiff(xw, ctx.trsb)
    elseif r == EnableMechanism
        return setdiff((y ∪ setdiff(xw, x)) ∩ ctx.md_s, d)
    elseif r == ExchangeProxy
        return (y ∪ xw) ∩ ctx.md_p
    end
    return VarSet()
end

# Missing data restriction shared by the rules that add variables: `z` may not
# add a true variable whose proxy is present or vice versa, nor both versions.
function md_addition_ok(a::VarSet, b::VarSet, z::VarSet, ctx::RuleContext)
    ab = a ∪ b
    isempty(((ab ∩ ctx.md_p) >> 2) ∩ z) || return false   # x when x* exists
    isempty(((ab ∩ ctx.md_t) << 2) ∩ z) || return false   # x* when x exists
    isempty(((z ∩ ctx.md_p) >> 2) ∩ z) || return false    # both x and x*
    isempty(((z ∩ ctx.md_t) << 2) ∩ z) || return false
    return true
end

"""
    apply_rule(rule, p, z, ctx) -> Union{Inference,Nothing}

Apply `rule` to term `p` with the subset `z`, or return `nothing` if the rule
does not apply. This is `apply_rule` followed by `get_ruleinfo` in
dosearch.cpp. The do-calculus rules use the intervention nodes `z << n`.
"""
function apply_rule(r::Rule, p::Dist, z::VarSet, ctx::RuleContext)
    y, xw, x, d = p.a, p.b, p.c, p.d   # P(y | do(x), xw \ x), enabled d
    n = ctx.n
    trsb, md_p, md_s = ctx.trsb, ctx.md_p, ctx.md_s
    E = VarSet()
    if r == InsertObservation
        isempty(z ∩ y) && isempty(z ∩ xw) || return nothing
        ctx.md && !md_addition_ok(y, xw, z, ctx) && return nothing
        to = Dist(y, xw ∪ z, x, d)
        return Inference(to, NO_DIST, Separation(y, z, xw, x), false)
    elseif r == DeleteObservation
        isempty(z ∩ y) && isempty(z ∩ x) && z ⊆ setdiff(xw, x) || return nothing
        xwmz = setdiff(xw, z)
        to = Dist(y, xwmz, x, setdiff(d, z))
        return Inference(to, NO_DIST, Separation(y, z, xwmz, x), false)
    elseif r == ActionToObservation
        isempty(z ∩ y) && isempty(z ∩ setdiff(xw, x)) && z ⊆ x && isempty(z ∩ md_p) || return nothing
        xmz = setdiff(x, z)
        to = Dist(y, xw, xmz, d)
        return Inference(to, NO_DIST, Separation(y, z << n, xw, xmz), false)
    elseif r == ObservationToAction
        isempty(z ∩ y) && isempty(z ∩ x) && z ⊆ xw && isempty(z ∩ md_p) || return nothing
        to = Dist(y, xw, x ∪ z, d ∪ (z ∩ md_s))
        return Inference(to, NO_DIST, Separation(y, z << n, xw, x), false)
    elseif r == DeleteAction
        isempty(z ∩ y) && isempty(z ∩ setdiff(xw, x)) && z ⊆ x && isempty(z ∩ md_p) || return nothing
        xwmz = setdiff(xw, z)
        xmz = setdiff(x, z)
        to = Dist(y, xwmz, xmz, setdiff(d, z))
        return Inference(to, NO_DIST, Separation(y, z << n, xwmz, xmz), false)
    elseif r == InsertAction
        isempty(z ∩ y) && isempty(z ∩ xw) || return nothing
        if ctx.md
            isempty(z ∩ md_p) || return nothing
            # dosearch.cpp tests `j` twice here instead of `j` and then `k`, so
            # the check "no x* when x exists" is skipped for this rule.
            isempty((((y ∪ xw) ∩ md_p) >> 2) ∩ z) || return nothing
            isempty(((z ∩ md_p) >> 2) ∩ z) || return nothing
            isempty(((z ∩ ctx.md_t) << 2) ∩ z) || return nothing
        end
        to = Dist(y, xw ∪ z, x ∪ z, d ∪ (z ∩ md_s))
        return Inference(to, NO_DIST, Separation(y, z << n, xw, x), false)
    elseif r == Marginalize
        z ⊆ y && z != y && isempty(z ∩ trsb) || return nothing
        # cannot marginalize over an enabled response indicator
        ctx.md && !isempty(z ∩ y ∩ md_s) && !isempty(z ∩ d) && return nothing
        return Inference(Dist(setdiff(y, z), xw, x, d), NO_DIST, NO_SEPARATION, false)
    elseif r == Condition
        z ⊆ y && z != y && isempty(z ∩ trsb) || return nothing
        # cannot condition if a response indicator is enabled
        ctx.md && !isempty(y ∩ md_s) && !isempty(y ∩ d) && return nothing
        return Inference(Dist(setdiff(y, z), xw ∪ z, x, d), NO_DIST, NO_SEPARATION, false)
    elseif r == ProductConditional
        z ⊆ xw && isempty(z ∩ x) && isempty(z ∩ trsb) || return nothing
        xwmz = setdiff(xw, z)
        to = Dist(y ∪ z, xwmz, x, d)
        required = Dist(z, xwmz, x, setdiff(d, y))
        return Inference(to, required, NO_SEPARATION, false)
    elseif r == ProductMarginal
        isempty(z ∩ y) && isempty(z ∩ xw) && isempty(z ∩ trsb) || return nothing
        ctx.md && !md_addition_ok(y, xw, z, ctx) && return nothing
        to = Dist(y ∪ z, xw, x, d)
        required = Dist(z, xw ∪ y, x, d)
        return Inference(to, required, NO_SEPARATION, true)
    elseif r == DivideMarginal
        z ⊆ y && z ⊆ d && isempty(z ∩ trsb) && z != y || return nothing
        xwz = xw ∪ z
        to = Dist(setdiff(y, z), xwz, x, d)
        required = Dist(z, xw, x, d ∩ xwz)
        return Inference(to, required, NO_SEPARATION, false)
    elseif r == DivideConditional
        z ⊆ y && z ⊆ d && isempty(z ∩ trsb) && z != y || return nothing
        ymz = setdiff(y, z)
        to = Dist(ymz, xw, x, setdiff(d, z))
        required = Dist(z, xw ∪ ymz, x, d)
        return Inference(to, required, NO_SEPARATION, false)
    elseif r == DenominatorMarginal
        isempty(z ∩ y) && isempty(z ∩ xw) && isempty(z ∩ trsb) || return nothing
        to = Dist(z, xw ∪ y, x, d)
        required = Dist(y ∪ z, xw, x, d)
        return Inference(to, required, NO_SEPARATION, true)
    elseif r == DenominatorConditional
        z ⊆ xw && isempty(z ∩ trsb) || return nothing
        xwmz = setdiff(xw, z)
        to = Dist(z, xwmz, x, setdiff(d, y))
        required = Dist(y ∪ z, xwmz, x, d)
        return Inference(to, required, NO_SEPARATION, true)
    elseif r == EnableMechanism
        z ⊆ ((y ∩ md_s) ∪ (setdiff(xw, x) ∩ md_s)) && isempty(z ∩ d) || return nothing
        return Inference(Dist(y, xw, x, d ∪ z), NO_DIST, NO_SEPARATION, false)
    elseif r == ExchangeProxy
        u = y ∩ md_p                       # proxies on the left, p(u*|.)
        v = xw ∩ md_p                      # proxies on the right, p(.|v*)
        z ⊆ (u ∪ v) || return nothing      # z must contain only proxies
        u = u ∩ z
        v = v ∩ z
        j = y ∩ d                          # enabled indicators on the left
        k = xw ∩ d                         # enabled indicators on the right
        # A proxy u* on the left can be replaced if r_u is enabled anywhere,
        # a proxy v* on the right only if r_v is enabled on the right.
        u ⊆ ((j ∪ k) << 1) && v ⊆ (k << 1) || return nothing
        yz = y ∩ z
        xwz = xw ∩ z
        to = Dist(setdiff(y, yz) ∪ (yz >> 2), setdiff(xw, xwz) ∪ (xwz >> 2), x, d)
        return Inference(to, NO_DIST, NO_SEPARATION, false)
    end
    return nothing
end
