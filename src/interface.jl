# Public entry point and result type.

"""
    DosearchResult

The result of [`dosearch`](@ref).

- `identifiable::Bool`: whether the query was identified.
- `limit_reached::Bool`: whether the search stopped at `time_limit` or
  `term_limit` before it could decide; the query is then not identified, but
  may be identifiable.
- `formula::Union{Nothing,String}`: the identifying formula in LaTeX, in the
  notation of the R package, or `nothing` if the query was not identified or
  `formula = false`.
- `derivation::Union{Nothing,Derivation}`: the [`Derivation`](@ref) of the
  formula, which [`latex`](@ref) and [`typst`](@ref) write out step by step
  (`nothing` if there is no formula).
- `dot::Union{Nothing,String}`: the derivation as a DOT graph if
  `draw_derivation = true` (`derivation` in the R package).
- `time::Float64`: seconds spent searching.
- `terms::Int`: the number of terms the search derived, input distributions
  included (what `term_limit` limits).
- `rule_times::Vector{Float64}`: seconds spent in each rule of `rules` if
  `benchmark_rules = true`.
- `rules`: the rules that the search used, in order.
- `problem::DosearchProblem`: the problem that was solved.
"""
struct DosearchResult{P<:DosearchProblem,R}
    identifiable::Bool
    limit_reached::Bool
    formula::Union{Nothing,String}
    derivation::Union{Nothing,Derivation}
    dot::Union{Nothing,String}
    time::Float64
    terms::Int
    rule_times::Vector{Float64}
    rules::Vector{R}
    problem::P
end

"""
    isidentifiable(result::DosearchResult) -> Bool

Whether the query of `result` was identified.
"""
isidentifiable(r::DosearchResult) = r.identifiable

to_rule(r::Rule) = r
to_rule(r::Integer) = rule_from_code(Rule, r, "DAGs")
to_rule(r) = throw(ArgumentError("Invalid rule $(repr(r)): give rules as integer codes or as `Rule`s."))

# The rule of an R rule code, with an error that lists the valid codes.
function rule_from_code(::Type{T}, code::Integer, models::String) where {T}
    input = T === Rule ? Int(InputDistribution) : Int(CsiInputDistribution)
    valid = sort!([Int(r) for r in instances(T) if Int(r) != input])
    (code in valid || code == input) ||
        throw(ArgumentError("Unknown rule code $code. The rule codes of $models are " *
                            "$(join(valid, ", ")); see `$T`."))
    return T(code)
end

# Keyword arguments, for errors about unknown ones.
const MODEL_KEYWORDS = (:transportability, :selection_bias, :missing_data, :md_sym, :context_vars, :warn)
const OUTPUT_KEYWORDS = (:formula, :draw_derivation, :draw_all)
const SEARCH_KEYWORDS = (:heuristic, :strategy, :improve, :rules, :time_limit, :term_limit,
                         :benchmark_rules, :verbose, :io, :cache)

function check_keywords(kwargs, allowed, all_keywords = allowed)
    for k in keys(kwargs)
        k in allowed && continue
        k in MODEL_KEYWORDS && throw(ArgumentError(
            "`$k` describes the model, so it is given when the problem is built: to " *
            "`DosearchProblem` or to the methods that take `data`, `query` and `graph`."))
        limit = min(2, (length(string(k)) - 1) ÷ 2)   # a close misspelling
        close = filter(a -> edit_distance(string(k), string(a)) <= limit, collect(all_keywords))
        hint = isempty(close) ? "" : " Did you mean `$(first(close))`?"
        throw(ArgumentError("Unknown keyword argument `$k`.$hint The keyword arguments are " *
                            join(("`$a`" for a in all_keywords), ", ") * "."))
    end
    return nothing
end

# Levenshtein distance of two strings.
function edit_distance(a::AbstractString, b::AbstractString)
    x, y = collect(a), collect(b)
    d = collect(0:length(y))
    for i in eachindex(x)
        prev, d[1] = d[1], i
        for j in eachindex(y)
            cur = d[j + 1]
            d[j + 1] = min(d[j + 1] + 1, d[j] + 1, prev + (x[i] != y[j]))
            prev = cur
        end
    end
    return d[end]
end

function rule_list(rules, problem::DagProblem, improve::Bool)
    if rules === nothing || isempty(rules)
        return default_rules(; missing_data = has_missing_data(problem), improve)
    end
    list = rules isa Union{Rule,Integer} ? [to_rule(rules)] : Rule[to_rule(r) for r in rules]
    InputDistribution in list && throw(ArgumentError("`InputDistribution` (code 0) is not a rule."))
    return list
end

"""
    dosearch(data, query, graph; kwargs...) -> DosearchResult
    dosearch(problem::DosearchProblem; kwargs...) -> DosearchResult

Identify the causal `query` from the input distributions `data` in the causal
model `graph`, by a search over do-calculus and probability manipulations
(Tikka, Hyttinen and Karvanen, 2021). The graph is either a DAG with latent
confounders or a labeled DAG (LDAG) with context-specific independences
(Tikka, Hyttinen and Karvanen, 2019). This is a port of `dosearch` of the R
package of the same name, and gives the same formulas.

# Syntax

`data` lists the input distributions, one per line, as a string or a vector of
strings. A distribution is written `p(A|do(B),C)`: `A`, `B` and `C` are comma
separated variable names, and both the conditioning variables `C` and the
do-operator are optional (`p(A|C,do(B))` is also allowed). For example

```julia
data = \"\"\"
    p(z|do(x))
    p(w,y|do(z,x))
    p(w,y,x|z)
\"\"\"
```

A distribution can also be given as a `NamedTuple` of variable roles: `0` for
the left-hand side, `1` for intervened and `2` for conditioning variables, so
`(y = 0, x = 1)` is `p(y|do(x))`.

`query` is the target distribution in the same syntax.

`graph` has one edge per line: `x -> y` for a directed edge and `x <-> y` (or
`x -- y`) for a bidirected edge, that is, an unobserved common cause of `x`
and `y`. The parts of an edge are separated by whitespace.

## LDAGs

A graph with edge labels is an LDAG. A label is a comma separated list of
value assignments of other parents of the child, and the edge vanishes in the
context that matches it. An edge that vanishes in several contexts is listed
once for each (or, as in the R package, with the labels separated by `;`):

```julia
graph = \"\"\"
    X -> Y : A = 1
    Z -> Y : A = 0
    Z -> Y : B = 1, C = 0
\"\"\"
```

LDAGs have no bidirected edges: latent variables are explicit nodes that do
not appear in the data. Interventions are explicit intervention nodes whose
names start with `I_`, such as `I_X -> X`, and distributions condition on
them instead of using the do-operator: `p(Y|X,I_X=1)`. In general,
variables in distributions of an LDAG can be assigned `0` or `1`.

## Values

For problems built in code, the inputs can be Julia values instead of text:
`data` a [`P`](@ref) or a vector of them (mixed with text if needed), `query`
a `P`, and `graph` a [`CausalGraph`](@ref).

# Keyword arguments

Model:

- `transportability`: comma separated transportability nodes, which mark the
  discrepancies between domains, or a collection of their names.
- `selection_bias`: comma separated selection bias nodes, or a collection of
  their names.
- `missing_data`: comma separated missing data mechanisms `r_x : x`, where
  `r_x` is the response indicator of `x`, or pairs `"r_x" => "x"`. The proxy
  variable `x*` of `x` is added automatically. A response indicator can be fixed to the observed value
  in distributions as `r_x = 1`.
- `md_sym = "1"`: the value that marks an enabled response indicator.
- `context_vars`: names of further context variables of an LDAG, besides those
  in the edge labels and those fixed to a value in the query or an input.
- `warn = true`: warn about inputs that are valid but probably unintended.

Transportability, selection bias and missing data are only supported for
DAGs.

Search:

- `heuristic = false`: expand the most promising terms first instead of in the
  order they were found.
- `strategy`: the [`SearchStrategy`](@ref) that orders the search,
  [`ScoreHeuristic()`](@ref ScoreHeuristic) if `heuristic = true` and
  [`FoundOrder()`](@ref FoundOrder) otherwise.
- `improve = true`: leave out redundant rules, and stop early when the query is
  trivially not identifiable.
- `rules`: the rules to use, in order, as R rule codes or as [`Rule`](@ref)s
  for DAGs and [`CsiRule`](@ref)s for LDAGs. Defaults to
  [`default_rules`](@ref) and [`default_csi_rules`](@ref).
- `time_limit = Inf`: stop the search after this many seconds (checked before
  each term is expanded; `0` stops at once). R gives the limit in hours.
- `term_limit`: stop the search once it has derived this many terms (input
  distributions included; see the field `terms` of the result), checked
  before each term is expanded. Unlike `time_limit`, the result does not
  depend on the machine.
- `cache = true`: for LDAGs, reuse context-specific independences that were
  found to hold, as the R package does by default.

Output:

- `formula = true`: compute the identifying formula and its derivation (see
  [`latex`](@ref) and [`typst`](@ref)).
- `draw_derivation = false`: return the derivation as a DOT graph, in the
  field `dot`.
- `draw_all = false`: draw every derivation step, not only the steps that lead
  to the query.
- `benchmark_rules = false`: record the time spent in each rule.
- `verbose = false`: print the steps of the search, to `io = stdout`.

Differences from the R package: times are in seconds, `draw_all` draws the
steps in the order they were found, and up to 32 variables are supported in
DAGs (R: 15) and 64 in LDAGs (R: 30). Edges of an LDAG that vanish in every
context are left out of the contexts too; R keeps them there by mistake.

# Examples

```jldoctest
julia> r = dosearch("p(x,y,z)", "p(y|do(x))", \"\"\"
           x -> y
           z -> x
           z -> y
       \"\"\");

julia> r.formula
"\\\\sum_{z}\\\\left(p(z)p(y|x,z)\\\\right)"
```
"""
function dosearch(data, query, graph::AbstractString;
                  transportability = nothing, selection_bias = nothing, missing_data = nothing,
                  md_sym::AbstractString = "1", context_vars = String[], warn::Bool = true, kwargs...)
    problem = DosearchProblem(data, query, graph; transportability, selection_bias,
                              missing_data, md_sym, context_vars, warn)
    return dosearch(problem; kwargs...)
end

function dosearch(problem::DosearchProblem; formula::Bool = true, draw_derivation::Bool = false,
                  draw_all::Bool = false, kwargs...)
    check_keywords(kwargs, SEARCH_KEYWORDS, (MODEL_KEYWORDS..., OUTPUT_KEYWORDS..., SEARCH_KEYWORDS...))
    st, elapsed = solve(problem; kwargs...)
    identifiable = st.target !== nothing
    shortcut = st.options.improve && trivially_nonidentifiable(st)
    derivation = identifiable && formula ? Derivation(st) : nothing
    f = derivation === nothing ? nothing : derivation.formula
    dot = !draw_derivation ? nothing :
          shortcut ? "" : derivation_dot(st; draw_all)
    rule_times = st.options.benchmark_rules ? st.rule_times : Float64[]
    return DosearchResult(identifiable, st.limit_reached, f, derivation, dot, elapsed, st.count,
                          rule_times, st.options.rules, problem)
end

latex(r::DosearchResult; kwargs...) = latex(derivation_of(r); kwargs...)
typst(r::DosearchResult; kwargs...) = typst(derivation_of(r); kwargs...)

search_model(p::DagProblem, rules::Vector{Rule}; cache::Bool = true) = DagModel(p, rules)
search_model(p::LdagProblem, rules::Vector{CsiRule}; cache::Bool = true) = LdagModel(p; cache)

to_csi_rule(r::CsiRule) = r
to_csi_rule(r::Integer) = rule_from_code(CsiRule, r, "LDAGs")
to_csi_rule(r) = throw(ArgumentError("Invalid rule $(repr(r)): give rules as integer codes or as `CsiRule`s."))

function rule_list(rules, problem::LdagProblem, improve::Bool)
    (rules === nothing || isempty(rules)) && return default_csi_rules()
    list = rules isa Union{CsiRule,Integer} ? [to_csi_rule(rules)] : CsiRule[to_csi_rule(r) for r in rules]
    CsiInputDistribution in list && throw(ArgumentError("`CsiInputDistribution` is not a rule."))
    return list
end

# Run the search and return its final state and the time it took.
function solve(problem::DosearchProblem; heuristic::Bool = false,
               strategy::SearchStrategy = heuristic ? ScoreHeuristic() : FoundOrder(),
               improve::Bool = true, rules = nothing, time_limit::Real = Inf,
               term_limit::Integer = typemax(Int), benchmark_rules::Bool = false,
               verbose::Bool = false, io::IO = stdout, cache::Bool = true)
    time_limit >= 0 || throw(ArgumentError("`time_limit` must be nonnegative, got $time_limit; Inf means no limit."))
    term_limit >= 0 || throw(ArgumentError("`term_limit` must be nonnegative, got $term_limit."))
    ruleset = rule_list(rules, problem, improve)
    options = SearchOptions(; rules = ruleset, time_limit = Float64(time_limit),
                            term_limit = Int(term_limit), benchmark_rules, improve, verbose, io)
    st = SearchState(search_model(problem, ruleset; cache), strategy, options)
    verbose && println(io, "Setting target: ", term_string(problem, problem.query))
    elapsed = @elapsed run_search!(st)
    return st, elapsed
end

function Base.show(io::IO, ::MIME"text/plain", r::DosearchResult)
    if r.limit_reached
        print(io, "The query ", r.problem.query_text, " was not identified before the search reached its limit.")
    elseif r.formula === nothing
        print(io, "The query ", r.problem.query_text, " is ",
              r.identifiable ? "identifiable." : "not identifiable by do-search.")
    else
        print(io, r.formula)
    end
end

function Base.show(io::IO, r::DosearchResult)
    print(io, "DosearchResult(identifiable = ", r.identifiable, ")")
end
