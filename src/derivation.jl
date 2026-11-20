# The derivation of an identified query, step by step: start from the query and
# rewrite one term at a time by the rule that derived it, until only known terms
# remain and the expression is the identifying formula. Each step names its
# rule, with the d-separation or context-specific independence that justifies
# it. Rendered as LaTeX, Typst or plain text.

"""
    Derivation

How the search derived an identified query. Stored in the `derivation` field
of a [`DosearchResult`](@ref).

Written out, the derivation starts from the query and rewrites one term at a
time by the rule that derived it, until only known terms remain. The last line
is then the identifying formula. Each step names its rule. A do-calculus step
also gives the d-separation that licenses it, in the mutilated graph of the
rule. An LDAG step gives the context-specific independence instead.

`latex(d)` and `typst(d)` write the derivation for a document, and `d`
displays as plain text. The keyword `inputs = true` continues the derivation
past the terms that the formula keeps. Those terms are computed from an input
distribution by probability manipulations alone, and the rewriting stops only
at the input distributions themselves.

The fields are internal: the problem, the derivation's nodes of the search
(node 1 is the query), the rewrite of each derived node and the identifying
formula in the notation of the R package.
"""
struct Derivation{P<:DosearchProblem,N<:Node}
    problem::P
    nodes::Vector{N}
    rewrites::Vector{Formula}
    formula::String
end

"""
    Derivation(st::SearchState)

The derivation of the query found by the search `st`.
"""
function Derivation(st::SearchState)
    st.target === nothing && throw(ArgumentError("the query was not identified"))
    search_node(i) = i == 0 ? st.target : st.nodes[i]
    # The nodes that the query was derived from, by search index (the query
    # has none, 0), numbered in breadth-first order from the query.
    sources = [0]
    index = Dict(0 => 1)
    k = 1
    while k <= length(sources)
        node = search_node(sources[k])
        for j in (node.parent1, node.parent2)
            if j > 0 && !haskey(index, j)
                push!(sources, j)
                index[j] = length(sources)
            end
        end
        k += 1
    end
    renumber(i) = i == 0 ? 0 : index[i]
    nodes = map(sources) do s
        n = search_node(s)
        Node(n.dist, n.rule, renumber(n.parent1), renumber(n.parent2), n.primitive, n.sep)
    end
    p = problem(st)
    rf = RFormulas(p, nodes)
    rlength(i) = ncodeunits(rf(i))
    rewrites = Formula[n.parent1 == 0 ? TermFormula(i) : step_rewrite(p, nodes, n, rlength)
                       for (i, n) in enumerate(nodes)]
    return Derivation(p, nodes, rewrites, rf(1))
end

# The rewrite of a derived node in a derivation. A do-calculus step that
# intervenes on response indicators also enables them, as in the R package:
# the term it rewrites to is evaluated where they are enabled. The R formula
# leaves this evaluation implicit; the derivation writes it out.
step_rewrite(p, nodes::AbstractVector{<:Node}, n::Node, rlength) = rewrite(p, nodes, n, rlength)
function step_rewrite(p::DagProblem, nodes::AbstractVector{<:Node}, n::Node{Rule}, rlength)
    f = rewrite(p, nodes, n, rlength)
    isdocalculus(n.rule) || return f
    parent = nodes[n.parent1].dist
    enabled = setdiff(n.dist.d, parent.d) ∩ (parent.a ∪ parent.b)
    isempty(enabled) && return f
    return EvaluationFormula(f, enabled, p.md_sym, false)
end

"""
    rewrite_steps(d::Derivation; inputs = false) -> Vector{Tuple{Int,Formula}}

The steps of the derivation, from the query to the identifying formula. Each
step rewrites every occurrence of one term by its rule. It is returned as the
term's node and the whole expression after the step. With `inputs = true`, the
terms computed from the inputs by probability manipulations are rewritten too.
"""
function rewrite_steps(d::Derivation; inputs::Bool = false)
    known(i) = inputs ? d.nodes[i].parent1 == 0 : d.nodes[i].primitive
    expr = TermFormula(1)
    steps = Tuple{Int,Formula}[]
    while (j = first_unknown(expr, known)) > 0
        expr = substitute(expr, j, d.rewrites[j])
        push!(steps, (j, expr))
    end
    return steps
end

# The first term that is not known, in reading order (0 if none).
first_unknown(f::TermFormula, known) = known(f.node) ? 0 : f.node
first_unknown(f::SumFormula, known) = first_unknown(f.body, known)
first_unknown(f::EvaluationFormula, known) = first_unknown(f.body, known)
first_unknown(f::Union{ProductFormula,DifferenceFormula}, known) = first_unknown2(f.left, f.right, known)
first_unknown(f::FractionFormula, known) = first_unknown2(f.num, f.den, known)
first_unknown(f::CasesFormula, known) = first_unknown2(f.case0, f.case1, known)
function first_unknown2(f::Formula, g::Formula, known)
    j = first_unknown(f, known)
    return j > 0 ? j : first_unknown(g, known)
end

substitute(f::TermFormula, j::Int, g::Formula) = f.node == j ? g : f
substitute(f::SumFormula, j::Int, g::Formula) = SumFormula(f.vars, substitute(f.body, j, g))
substitute(f::ProductFormula, j::Int, g::Formula) = ProductFormula(substitute(f.left, j, g), substitute(f.right, j, g))
substitute(f::FractionFormula, j::Int, g::Formula) = FractionFormula(substitute(f.num, j, g), substitute(f.den, j, g))
substitute(f::DifferenceFormula, j::Int, g::Formula) =
    DifferenceFormula(substitute(f.left, j, g), substitute(f.right, j, g))
substitute(f::CasesFormula, j::Int, g::Formula) =
    CasesFormula(f.var, substitute(f.case0, j, g), substitute(f.case1, j, g), f.zero_first)
substitute(f::EvaluationFormula, j::Int, g::Formula) =
    EvaluationFormula(substitute(f.body, j, g), f.vars, f.value, f.explicit)

# Justifications -----------------------------------------------------------------

abstract type Justification end

# A step of probability calculus, by name.
struct NamedStep <: Justification
    name::String
end

# A do-calculus step, licensed by `(y ⫫ z | given)` in the graph whose edges
# into `overline` and out of `underline` are removed (Pearl, 1995). Rule 3 also
# removes the edges into `overline_z`, the nodes of `z` that are not
# ancestors of the observed nodes of `given` (Z(W) in Pearl's notation).
struct DoCalculusStep <: Justification
    rule::Int
    y::VarSet
    z::VarSet
    given::VarSet
    overline::VarSet
    underline::VarSet
    overline_z::VarSet
end

# A step licensed by the context-specific independence `(y ⫫ z | given)`,
# where the variables of `zero` and `one` take those values.
struct CsiStep <: Justification
    y::VarSet
    z::VarSet
    given::VarSet
    zero::VarSet
    one::VarSet
end

# The proxies `x*` that equal their true variables `x` where `r_x` is enabled.
struct ProxyStep <: Justification
    proxies::VarSet
end

# Ancestors of `w` (including `w`) in `g` without the edges into `intervened`.
function ancestors(g::DsepGraph, w::VarSet, intervened::VarSet)
    vars = fullset(g.nvars).bits
    anc = w.bits
    todo = w.bits
    while todo != 0
        v = trailing_zeros(todo) + 1
        todo &= todo - one(UInt64)
        v in intervened && continue
        new = g.parents[v] & vars & ~anc
        anc |= new
        todo |= new
    end
    return VarSet(anc)
end

# The justification of the step that rewrites `node` (the step goes from the
# derived term to the terms it was derived from).
function justification(p::DagProblem, graph, nodes, node::Node{Rule})
    code = abs(Int(node.rule))
    s = node.sep
    n = problem_nvars(p)
    if code == 1
        return DoCalculusStep(1, s.x, s.y, s.given, s.intervened, VarSet(), VarSet())
    elseif code == 2
        z = s.y >> n
        return DoCalculusStep(2, s.x, z, setdiff(s.given, z), s.intervened, z, VarSet())
    elseif code == 3
        z = s.y >> n
        w = setdiff(s.given, s.intervened)
        zw = setdiff(z, ancestors(graph, w, s.intervened))
        return DoCalculusStep(3, s.x, z, s.given, s.intervened, VarSet(), zw)
    elseif code == 4
        return NamedStep("marginalization")
    elseif code == 6
        return NamedStep("chain rule")
    elseif code == 9
        return NamedStep("evaluation")
    elseif code == 10
        parent = nodes[node.parent1].dist
        exchanged = setdiff(parent.a ∪ parent.b, node.dist.a ∪ node.dist.b) ∩ p.md_proxies
        return ProxyStep(exchanged)
    end
    return NamedStep("conditioning")              # rules 5, 7 and 8
end

function justification(::LdagProblem, graph, nodes, node::Node{CsiRule})
    code = abs(Int(node.rule))
    if code == 3 || code == 4
        s = node.sep
        return CsiStep(s.y, s.x, s.z, s.zero ∩ s.z, s.one ∩ s.z)
    end
    code == 0 && return NamedStep("marginalization")
    code == 1 && return NamedStep("conditioning")
    code == 2 && return NamedStep("chain rule")
    code == 5 && return NamedStep("cases")
    code == 8 && return NamedStep("evaluation")
    return NamedStep("law of total probability")   # rules 6 and 7
end

justification_graph(p::DagProblem) = causal_graph(p)
justification_graph(::LdagProblem) = nothing

# Notations ------------------------------------------------------------------------

abstract type Notation end

struct LatexNotation <: Notation
    labels::Dict{String,String}
end

struct TypstNotation <: Notation
    labels::Dict{String,String}
end

struct TextNotation <: Notation end

const MarkupNotation = Union{LatexNotation,TypstNotation}

# Bound variables of sums are renamed with primes where their name is already in
# use (Σ_{x'} instead of a second x). `primes` holds the renamed variables and
# `names` the names in use, as plain text.
struct Scope
    primes::Dict{Int,Int}
    names::Set{String}
end
Scope() = Scope(Dict{Int,Int}(), Set{String}())

primed(s::AbstractString, k::Int) = s * "'"^k

# The scope inside a sum over `vars`.
function bind(scope::Scope, p, vars::VarSet)
    inner = Scope(copy(scope.primes), copy(scope.names))
    for i in vars
        k = 0
        while primed(p.labels[i], k) in inner.names || (k > 0 && primed(p.labels[i], k) in p.labels)
            k += 1
        end
        inner.primes[i] = k
        push!(inner.names, primed(p.labels[i], k))
    end
    return inner
end

label(::TextNotation, s::AbstractString, k::Int = 0) = primed(s, k)
function label(n::MarkupNotation, s::AbstractString, k::Int = 0)
    haskey(n.labels, s) || return default_label(n, s, k)
    l = n.labels[s]
    k == 0 && return l
    # A prime after a superscript would be a double superscript in LaTeX.
    return n isa LatexNotation && occursin('^', l) ? "{" * l * "}" * "'"^k : l * "'"^k
end

# Variable names in math: a trailing `*` (proxy variables) becomes a superscript
# and trailing `'`s primes, trailing Unicode subscript digits (`A₁`) and the
# part after the first `_` a subscript, and names longer than one letter are
# set as a single italic word.
function default_label(n::MarkupNotation, s::AbstractString, k::Int = 0)
    if length(s) > 1 && endswith(s, '*')
        base = default_label(n, chop(s))
        return n isa LatexNotation ? base * "^{*" * "\\prime"^k * "}" : base * "^(*" * "'"^k * ")"
    elseif length(s) > 1 && endswith(s, '\'')
        return default_label(n, chop(s), k + 1)
    end
    # the subscript after `_` is a label itself, so `R_Y₁` is R_{Y_{1}}
    i = findfirst('_', s)
    if i !== nothing && i > firstindex(s) && i < lastindex(s)
        base = name_markup(n, s[firstindex(s):prevind(s, i)]) * "'"^k
        sub = default_label(n, s[nextind(s, i):end])
        return n isa LatexNotation ? base * "_{" * sub * "}" : base * "_(" * sub * ")"
    end
    if length(s) > 1 && is_subscript_digit(last(s))
        base = rstrip(is_subscript_digit, s)
        if !isempty(base)
            digits = map(c -> Char('0' + (c - '₀')), s[nextind(s, lastindex(base)):end])
            b = default_label(n, base, k)
            return n isa LatexNotation ? b * "_{" * digits * "}" : b * "_(" * digits * ")"
        end
    end
    return name_markup(n, s) * "'"^k
end

plain_name(s::AbstractString) = (length(s) == 1 && isletter(s[1])) || all(isdigit, s)
is_subscript_digit(c::Char) = '₀' <= c <= '₉'

function name_markup(::LatexNotation, s::AbstractString)
    plain_name(s) && return String(s)
    escaped = replace(s, '\\' => "\\backslash{}", '{' => "\\{", '}' => "\\}", '#' => "\\#",
                      '$' => "\\\$", '%' => "\\%", '&' => "\\&", '_' => "\\_", '^' => "\\hat{}",
                      '~' => "\\sim{}")
    return "\\mathit{" * escaped * "}"
end

function name_markup(::TypstNotation, s::AbstractString)
    plain_name(s) && return String(s)
    return "italic(\"" * replace(s, '\\' => "\\\\", '"' => "\\\"") * "\")"
end

separator(::TextNotation) = ","
separator(::MarkupNotation) = ", "

const NO_SCOPE = Scope()

# The labels of `s`, each followed by ` = v` if `value(i)` gives a value `v`.
function var_list(n::Notation, p, s::VarSet, value = i -> nothing, scope::Scope = NO_SCOPE)
    parts = String[]
    for i in s
        l = label(n, p.labels[i], get(scope.primes, i, 0))
        v = value(i)
        push!(parts, v === nothing ? l : l * " = " * v)
    end
    return join(parts, separator(n))
end

# A set in an independence statement: braces around more than one variable.
function set_markup(n::Notation, p, s::VarSet, value = i -> nothing)
    list = var_list(n, p, s, value)
    length(s) > 1 || return list
    return n isa LatexNotation ? "\\{" * list * "\\}" : "{" * list * "}"
end

# The bar of a term and of an independence statement.
bar(::LatexNotation) = " \\mid "
bar(::TypstNotation) = " | "
bar(::TextNotation) = "|"
given_bar(n::Notation) = bar(n)
given_bar(::TextNotation) = " | "

do_markup(::LatexNotation, list) = "\\mathrm{do}(" * list * ")"
do_markup(::TypstNotation, list) = "\"do\"(" * list * ")"
do_markup(::TextNotation, list) = "do(" * list * ")"

# A term; in plain text as in the R package (`term_string`).
function term_markup(n::Notation, p::DagProblem, t::Dist, scope::Scope)
    enabled(i) = i in t.d ? p.md_sym : nothing
    s = "p(" * var_list(n, p, t.a, enabled, scope)
    if !isempty(t.b)
        parts = String[]
        isempty(t.c) || push!(parts, do_markup(n, var_list(n, p, t.c, enabled, scope)))
        w = setdiff(t.b, t.c)
        isempty(w) || push!(parts, var_list(n, p, w, enabled, scope))
        s *= bar(n) * join(parts, separator(n))
    end
    return s * ")"
end

function term_markup(n::Notation, p::LdagProblem, t::Dist, scope::Scope)
    s = "p(" * var_list(n, p, t.a, assignment(t), scope)
    isempty(t.b) || (s *= bar(n) * var_list(n, p, t.b, assignment(t), scope))
    return s * ")"
end

assignment(t::Dist) = i -> i in t.c ? "0" : i in t.d ? "1" : nothing

paren(::LatexNotation, s) = "\\left(" * s * "\\right)"
paren(::Notation, s) = "(" * s * ")"

text_markup(::LatexNotation, s) = "\\text{" * s * "}"
text_markup(::TypstNotation, s) = "\"" * s * "\""
text_markup(::TextNotation, s) = s

# Formulas ------------------------------------------------------------------------

# Does the formula end with a sum that would extend over whatever follows it?
open_right(::Formula) = false
open_right(::SumFormula) = true
open_right(f::ProductFormula) = open_right(f.right)
open_right(f::DifferenceFormula) = !(f.right isa DifferenceFormula) && open_right(f.right)

factors(f::Formula) = Formula[f]
factors(f::ProductFormula) = vcat(factors(f.left), factors(f.right))

# An expression of the derivation `d`. The names of the query's variables are
# in use from the start.
function markup(n::Notation, d::Derivation, f::Formula)
    p = d.problem
    q = d.nodes[1].dist
    return markup(n, d, f, Scope(Dict{Int,Int}(), Set{String}(p.labels[i] for i in q.a ∪ q.b)))
end

markup(n::Notation, d::Derivation, f::TermFormula, scope::Scope) =
    term_markup(n, d.problem, d.nodes[f.node].dist, scope)

function markup(n::Notation, d::Derivation, f::SumFormula, scope::Scope)
    inner = bind(scope, d.problem, f.vars)
    body = markup(n, d, f.body, inner)
    f.body isa DifferenceFormula && (body = paren(n, body))
    return sum_markup(n, var_list(n, d.problem, f.vars, i -> nothing, inner), length(f.vars)) * " " * body
end

sum_markup(::LatexNotation, list, k) = "\\sum_{" * list * "}"
sum_markup(::TypstNotation, list, k) = "sum_(" * list * ")"
sum_markup(::TextNotation, list, k) = k == 1 ? "Σ_" * list : "Σ_{" * list * "}"

function markup(n::Notation, d::Derivation, f::ProductFormula, scope::Scope)
    fs = factors(f)
    parts = map(enumerate(fs)) do (k, g)
        s = markup(n, d, g, scope)
        needs = g isa DifferenceFormula || (k < length(fs) && open_right(g)) ||
                (n isa TextNotation && g isa FractionFormula)
        needs ? paren(n, s) : s
    end
    return join(parts, " ")
end

function markup(n::Notation, d::Derivation, f::FractionFormula, scope::Scope)
    num, den = markup(n, d, f.num, scope), markup(n, d, f.den, scope)
    n isa LatexNotation && return "\\frac{" * num * "}{" * den * "}"
    n isa TypstNotation && return "frac(" * num * ", " * den * ")"
    delimited(g) = g isa TermFormula || g isa CasesFormula || g isa EvaluationFormula
    delimited(f.num) || (num = paren(n, num))
    delimited(f.den) || (den = paren(n, den))
    return num * " / " * den
end

function markup(n::Notation, d::Derivation, f::DifferenceFormula, scope::Scope)
    left, right = markup(n, d, f.left, scope), markup(n, d, f.right, scope)
    open_right(f.left) && (left = paren(n, left))
    f.right isa DifferenceFormula && (right = paren(n, right))
    return left * " - " * right
end

function markup(n::Notation, d::Derivation, f::CasesFormula, scope::Scope)
    p = d.problem
    c0, c1 = markup(n, d, f.case0, scope), markup(n, d, f.case1, scope)
    v0 = var_list(n, p, f.var, i -> "0", scope)
    v1 = var_list(n, p, f.var, i -> "1", scope)
    if n isa LatexNotation
        return "\\begin{cases} " * c0 * " & \\text{if } " * v0 * " \\\\ " * c1 * " & \\text{if } " * v1 * " \\end{cases}"
    elseif n isa TypstNotation
        return "cases(" * c0 * " & quad \"if\" " * v0 * ", " * c1 * " & quad \"if\" " * v1 * ")"
    end
    return "{" * c0 * " if " * v0 * "; " * c1 * " if " * v1 * "}"
end

function markup(n::Notation, d::Derivation, f::EvaluationFormula, scope::Scope)
    # A term of a DAG evaluated where response indicators are enabled is the
    # term with them enabled.
    if d.problem isa DagProblem && f.body isa TermFormula && f.value == d.problem.md_sym
        t = d.nodes[f.body.node].dist
        return term_markup(n, d.problem, Dist(t.a, t.b, t.c, t.d ∪ f.vars), scope)
    end
    body = markup(n, d, f.body, scope)
    at = var_list(n, d.problem, f.vars, i -> f.value, scope)
    n isa LatexNotation && return "\\left[" * body * "\\right]_{" * at * "}"
    n isa TypstNotation && return "lr([" * body * "])_(" * at * ")"
    return "[" * body * "]_{" * at * "}"
end

# Justifications ------------------------------------------------------------------

justification_markup(n::Notation, p, ::Any, j::NamedStep) = text_markup(n, j.name)

independence(::LatexNotation) = " \\perp\\!\\!\\!\\perp "
independence(::TypstNotation) = " ⫫ "
independence(::TextNotation) = " ⊥⊥ "

function independence_markup(n::Notation, p, y::VarSet, z::VarSet, given::VarSet, value = i -> nothing)
    s = "(" * set_markup(n, p, y, value) * independence(n) * set_markup(n, p, z, value)
    isempty(given) || (s *= given_bar(n) * var_list(n, p, given, value))
    return s * ")"
end

# The mutilated graph of a do-calculus rule.
function graph_markup(n::LatexNotation, p, j::DoCalculusStep)
    subs = String[]
    isempty(j.overline) || push!(subs, "\\overline{" * var_list(n, p, j.overline) * "}")
    isempty(j.underline) || push!(subs, "\\underline{" * var_list(n, p, j.underline) * "}")
    isempty(j.overline_z) || push!(subs, "\\overline{" * var_list(n, p, j.overline_z) * "}")
    return isempty(subs) ? "\\mathcal{G}" : "\\mathcal{G}_{" * join(subs, "\\,") * "}"
end

# A comma in a Typst function argument would start another argument.
function graph_markup(n::TypstNotation, p, j::DoCalculusStep)
    list(s) = replace(var_list(n, p, s), ", " => "\\, ")
    subs = String[]
    isempty(j.overline) || push!(subs, "overline(" * list(j.overline) * ")")
    isempty(j.underline) || push!(subs, "underline(" * list(j.underline) * ")")
    isempty(j.overline_z) || push!(subs, "overline(" * list(j.overline_z) * ")")
    return isempty(subs) ? "cal(G)" : "cal(G)_(" * join(subs, " ") * ")"
end

# Plain text marks the removed edges with combining overlines and underlines.
function graph_markup(n::TextNotation, p, j::DoCalculusStep)
    mark(s, c) = [join(string(ch, c) for ch in p.labels[i]) for i in s]
    subs = vcat(mark(j.overline, '̅'), mark(j.underline, '̲'), mark(j.overline_z, '̅'))
    return isempty(subs) ? "G" : "G_{" * join(subs, ",") * "}"
end

function justification_markup(n::Notation, p, ::Any, j::DoCalculusStep)
    ind = independence_markup(n, p, j.y, j.z, j.given)
    g = graph_markup(n, p, j)
    rule = "rule " * string(j.rule) * ":"
    n isa LatexNotation && return "\\text{" * rule * " } " * ind * "_{" * g * "}"
    n isa TypstNotation && return "\"" * rule * "\" " * ind * "_(" * g * ")"
    return rule * " " * ind * " in " * g
end

function justification_markup(n::Notation, p, ::Any, j::CsiStep)
    value(i) = i in j.zero ? "0" : i in j.one ? "1" : nothing
    ind = independence_markup(n, p, j.y, j.z, j.given, value)
    n isa LatexNotation && return "\\text{CSI: } " * ind
    n isa TypstNotation && return "\"CSI:\" " * ind
    return "CSI: " * ind
end

function justification_markup(n::Notation, p::DagProblem, ::Any, j::ProxyStep)
    parts = map(collect(j.proxies)) do k
        # the mechanism of variable x is numbered x, r_x, x*
        eq = label(n, p.labels[k]) * " = " * label(n, p.labels[k - 2])
        r = label(n, p.labels[k - 1]) * " = " * p.md_sym
        n isa LatexNotation ? eq * " \\text{ if } " * r :
        n isa TypstNotation ? eq * " \"if\" " * r : eq * " if " * r
    end
    return join(parts, ", ")
end

# Output -----------------------------------------------------------------------------

function notation_labels(labels)
    dict = Dict{String,String}()
    for (k, v) in pairs(labels)
        dict[String(k)] = String(v)
    end
    return dict
end

# The lines of the derivation: the query, then each step's expression and
# justification.
function derivation_lines(n::Notation, d::Derivation, inputs::Bool)
    p = d.problem
    graph = justification_graph(p)
    query = markup(n, d, TermFormula(1))
    lines = Tuple{String,String}[]
    for (j, expr) in rewrite_steps(d; inputs)
        just = justification(p, graph, d.nodes, d.nodes[j])
        push!(lines, (markup(n, d, expr), justification_markup(n, p, graph, just)))
    end
    return query, lines
end

"""
    latex(d::Derivation; inputs = false, labels = Dict()) -> String
    latex(r::DosearchResult; kwargs...) -> String

The derivation `d`, or the derivation of the identified query of `r`, as a
LaTeX `align*` environment (package `amsmath`). This adds methods to
`Markdown.latex`. The first line is the query.
Each following line rewrites one term and names the rule it uses in a second
column; the last line is the identifying formula. Long derivations may need
`\\allowdisplaybreaks` to break across pages.

With `inputs = true`, terms computed from an input distribution by
probability manipulations are rewritten too, down to the input distributions.
`labels` maps variable names to LaTeX code that replaces the default: one
letter as is, longer names in `\\mathit`, the part after the first `_` as a
subscript and a trailing `*` as a superscript.

# Examples

```jldoctest
julia> r = dosearch("p(x,y,z)", "p(y|do(x))", \"\"\"
           x -> y
           z -> x
           z -> y
       \"\"\");

julia> print(latex(r))
\\begin{align*}
p(y \\mid \\mathrm{do}(x)) &= \\sum_{z} p(y, z \\mid \\mathrm{do}(x)) && \\text{marginalization} \\\\
&= \\sum_{z} p(z \\mid \\mathrm{do}(x)) p(y \\mid \\mathrm{do}(x), z) && \\text{chain rule} \\\\
&= \\sum_{z} p(z) p(y \\mid \\mathrm{do}(x), z) && \\text{rule 3: } (z \\perp\\!\\!\\!\\perp x)_{\\mathcal{G}_{\\overline{x}}} \\\\
&= \\sum_{z} p(z) p(y \\mid x, z) && \\text{rule 2: } (y \\perp\\!\\!\\!\\perp x \\mid z)_{\\mathcal{G}_{\\underline{x}}}
\\end{align*}
```
"""
function latex(d::Derivation; inputs::Bool = false, labels = Dict{String,String}())
    query, lines = derivation_lines(LatexNotation(notation_labels(labels)), d, inputs)
    io = IOBuffer()
    println(io, "\\begin{align*}")
    if isempty(lines)
        println(io, query)
    end
    for (k, (expr, just)) in enumerate(lines)
        print(io, k == 1 ? query * " &= " : "&= ", expr, " && ", just)
        println(io, k < length(lines) ? " \\\\" : "")
    end
    println(io, "\\end{align*}")
    return String(take!(io))
end

"""
    typst(d::Derivation; inputs = false, labels = Dict()) -> String
    typst(r::DosearchResult; kwargs...) -> String

The derivation `d`, or the derivation of the identified query of `r`, as a
Typst block equation. The first line is the query. Each following line
rewrites one term and names the rule it uses in a second column; the last line
is the identifying formula. Long derivations may need
`#show math.equation: set block(breakable: true)` to break across pages.

With `inputs = true`, terms computed from an input distribution by
probability manipulations are rewritten too, down to the input distributions.
`labels` maps variable names to Typst code that replaces the default: one
letter as is, longer names in `italic`, the part after the first `_` as a
subscript and a trailing `*` as a superscript.
"""
function typst(d::Derivation; inputs::Bool = false, labels = Dict{String,String}())
    query, lines = derivation_lines(TypstNotation(notation_labels(labels)), d, inputs)
    io = IOBuffer()
    print(io, "\$ ")
    isempty(lines) && print(io, query)
    for (k, (expr, just)) in enumerate(lines)
        print(io, k == 1 ? query * " &= " : "  &= ", expr, " && quad ", just)
        k < length(lines) && println(io, " \\")
    end
    println(io, " \$")
    return String(take!(io))
end

function derivation_of(r)
    r.derivation === nothing || return r.derivation
    r.identifiable && throw(ArgumentError("the derivation was not kept; search with `formula = true`"))
    throw(ArgumentError("the query $(r.problem.query_text) is not identifiable by do-search"))
end

"""
    plain_text(d::Derivation; inputs = false) -> String

The derivation `d` in plain text, as `d` is displayed.
"""
function plain_text(d::Derivation; inputs::Bool = false)
    query, lines = derivation_lines(TextNotation(), d, inputs)
    io = IOBuffer()
    print(io, query)
    isempty(lines) && return String(take!(io))
    exprs = ["  = " * e for (e, _) in lines]
    width = min(maximum(textwidth, exprs), 72)
    for ((_, just), e) in zip(lines, exprs)
        print(io, "\n", e, " "^max(width - textwidth(e), 0), "  (", just, ")")
    end
    return String(take!(io))
end

Base.show(io::IO, ::MIME"text/plain", d::Derivation) = print(io, plain_text(d))
Base.show(io::IO, ::MIME"text/latex", d::Derivation) = print(io, latex(d))
Base.show(io::IO, d::Derivation) =
    print(io, "Derivation(", term_string(d.problem, d.nodes[1].dist), " = ", d.formula, ")")
