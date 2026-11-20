# Identification problems in labeled DAGs (LDAGs): the graph with its edge
# labels, the contexts that the labels imply, the input distributions and
# the query. Port of `get_derivation_ldag` and its helpers in the R package's
# `internal_ldag.R`.
#
# An edge label lists value assignments of other parents of the edge's child,
# separated by `;`, such as `X -> Y : A = 0; B = 1, C = 0`. The edge vanishes
# in every context that matches one of the assignments. Intervention nodes are
# explicit nodes whose names start with `I_`.

"""
    LdagContext

A context of an LDAG: the variables `zero` assigned to 0 and `one` assigned
to 1, the edges `removed` (pairs of variable indices) that vanish in it, and
the index of its `equivalence` class among the contexts of the same variables
(contexts that remove the same edges are equivalent).
"""
struct LdagContext
    zero::VarSet
    one::VarSet
    equivalence::Int
    removed::Vector{Tuple{Int,Int}}
end

"""
    ContextSet

The variables `vars` of a set of contexts, and all `contexts` over them.
"""
struct ContextSet
    vars::VarSet
    contexts::Vector{LdagContext}
end

"""
    LocalCsi

A context-specific independence that a label encodes directly: `x` is
independent of its child `y` given the other parents `z` of `y` in the
context `zero`, `one`.
"""
struct LocalCsi
    x::VarSet
    y::VarSet
    z::VarSet
    zero::VarSet
    one::VarSet
end

"""
    LdagProblem(data, query, graph; context_vars = String[], warn = true)

A validated identification problem in a labeled DAG: identify the `query`
from the input distributions `data` in the LDAG `graph`. See
[`dosearch`](@ref) for the syntax. `context_vars` names extra variables to
treat as context variables, besides those that appear in edge labels and
those fixed to a value in the query or an input.

In internal form, a term [`Dist`](@ref) over this problem has the left-hand
side `a`, the right-hand side `b`, and the variables assigned to 0 in `c` and
to 1 in `d`.
"""
struct LdagProblem <: DosearchProblem
    labels::Vector{String}
    edges::Vector{Tuple{Int,Int}}
    context_sets::Vector{ContextSet}
    local_csi::Vector{LocalCsi}
    context_vars::VarSet
    intervention_vars::VarSet
    data::Vector{Dist}
    query::Dist
    data_text::Vector{String}
    query_text::String
    graph_text::String
end

problem_nvars(p::LdagProblem) = length(p.labels)

# Assignments of `k` binary variables in the order of R's `expand.grid`: the
# first variable varies fastest, so row `r` assigns bit `i - 1` of `r - 1` to
# variable `i`.
grid_value(row::Int, i::Int) = ((row - 1) >> (i - 1)) & 1

# An edge label assignment: the parents assigned to 0 and to 1.
struct LabelRow
    zero::VarSet
    one::VarSet
    from::Int
    to::Int
    parents::VarSet
end

const LDAG_SYNTAX = "Write one edge per line, as `x -> y`, or with a label as `x -> y : z = 1`."
const LDAG_DISTRIBUTION_SYNTAX = "A distribution of an LDAG is written `p(y|x,z = 1)`: comma separated " *
                                 "variables, the conditioning variables after `|`, and values 0 or 1 assigned " *
                                 "with `=`; interventions are intervention nodes such as `I_x = 1`."

# `order` (internal) numbers the named variables first, in that order, so
# that a problem rebuilt with more edges keeps the numbering of the original.
function LdagProblem(data, query, graph::AbstractString;
                     context_vars = String[], warn::Bool = true, order::Vector{String} = String[])
    check_line_breaks(graph, "graph", "edge")
    lines = filter(!isempty, [String(strip(l)) for l in split(graph, r"\r|\n")])
    isempty(lines) && throw(ArgumentError("Invalid graph, the graph is empty. " * LDAG_SYNTAX))
    malformed = [!occursin(r"^[^:]+(:[^:]+)?$", l) for l in lines]
    if any(malformed)
        throw(ArgumentError("Invalid graph, malformed lines found: " * cs(lines[malformed]) * ". " * LDAG_SYNTAX))
    end
    parts = [split(l, ':') for l in lines]
    # One edge per line; an edge may be repeated, with or without a label.
    directed = Tuple{String,String}[]
    for p in parts
        pg = parse_graph(String(strip(first(p))))
        isempty(pg.bidirected) ||
            throw(ArgumentError("Invalid graph, bidirected edges are not supported for LDAGs."))
        append!(directed, pg.directed)
    end
    is_acyclic(directed) || throw(ArgumentError("Invalid graph, the graph contains cycles."))
    pg = ParsedGraph(directed, Tuple{String,String}[])
    dir_lhs = first.(pg.directed)
    dir_rhs = last.(pg.directed)
    vars = unique(vcat(dir_rhs, dir_lhs))
    # intervention nodes last
    interventions = Set(filter(v -> startswith(v, "I_"), dir_lhs))
    default_order = vcat(filter(!in(interventions), vars), filter(in(interventions), vars))
    nb = Numbering(vcat(order, filter(!in(Set(order)), default_order)))
    intv_names = filter(in(interventions), nb.vars)
    parents = Dict(v => String[] for v in nb.vars)
    for (from, to) in pg.directed
        from in parents[to] || push!(parents[to], from)
    end

    # Edge labels. An edge listed on several lines has the alternatives of all
    # its labels, as if they were separated by `;` on one line.
    labeled = Tuple{String,String}[]
    alternatives = Dict{Tuple{String,String},Vector{String}}()
    for (i, p) in enumerate(parts)
        length(p) == 2 || continue
        edge = pg.directed[i]
        haskey(alternatives, edge) || (push!(labeled, edge); alternatives[edge] = String[])
        append!(alternatives[edge], filter(!isempty, split(replace(p[2], r"\s+" => ""), ';')))
    end
    rows = LabelRow[]
    inferred = LabelRow[]
    vanishing = Set{Tuple{Int,Int}}()
    label_vars = String[]
    for (from, to) in labeled
        pa = filter(!=(from), parents[to])
        present = falses(1 << length(pa))
        for assignment in alternatives[(from, to)]
            tokens = split(assignment, ',')
            pairs = [split(t, '=') for t in tokens]
            names = [String(first(q)) for q in pairs]
            values = [length(q) == 2 ? String(q[2]) : "" for q in pairs]
            validate_label(from, to, pa, names, values)
            append!(label_vars, filter(v -> !startswith(v, "I_"), names))
            zero = names[values .== "0"]
            one = names[values .== "1"]
            push!(rows, LabelRow(varset(nb, zero), varset(nb, one), nb.nums[from], nb.nums[to], varset(nb, pa)))
            for r in eachindex(present)
                matches = all(v -> grid_value(r, something(findfirst(==(v), pa))) == 0, zero) &&
                          all(v -> grid_value(r, something(findfirst(==(v), pa))) == 1, one)
                matches && (present[r] = true)
            end
        end
        all(present) && push!(vanishing, (nb.nums[from], nb.nums[to]))
        infer_labels!(inferred, nb, present, from, to, pa)
    end
    input_rows = unique(vcat(rows, inferred))
    con_names = unique(label_vars)

    context_sets = ContextSet[]
    local_csi = LocalCsi[]
    for i in 2:(1 << length(con_names))
        sub = con_names[[grid_value(i, k) == 1 for k in eachindex(con_names)]]
        contexts = LdagContext[]
        classes = Vector{Vector{Tuple{Int,Int}}}()
        for j in 1:(1 << length(sub))
            zero = varset(nb, sub[[grid_value(j, k) == 0 for k in eachindex(sub)]])
            one = varset(nb, sub[[grid_value(j, k) == 1 for k in eachindex(sub)]])
            removed = Tuple{Int,Int}[]
            for row in input_rows
                (row.zero ⊆ zero && row.one ⊆ one) || continue
                (row.from, row.to) in vanishing && continue
                push!(removed, (row.from, row.to))
                if row.parents == zero ∪ one
                    push!(local_csi, LocalCsi(VarSet([row.from]), VarSet([row.to]), row.parents, zero, one))
                end
            end
            unique!(removed)
            class = findfirst(==(removed), classes)
            if class === nothing
                push!(classes, removed)
                class = length(classes)
            end
            push!(contexts, LdagContext(zero, one, class, removed))
        end
        # sets of contexts in which no edge vanishes are left out
        all(c -> isempty(c.removed), contexts) || push!(context_sets, ContextSet(varset(nb, sub), contexts))
    end
    for i in 2:(1 << length(intv_names))
        one = varset(nb, intv_names[[grid_value(i, k) == 1 for k in eachindex(intv_names)]])
        removed = unique([(row.from, row.to) for row in input_rows
                          if isempty(row.zero) && row.one ⊆ one && !((row.from, row.to) in vanishing)])
        push!(context_sets, ContextSet(one, [LdagContext(VarSet(), one, 1, removed)]))
    end
    edges = unique([(nb.nums[a], nb.nums[b]) for (a, b) in pg.directed if !((nb.nums[a], nb.nums[b]) in vanishing)])

    # distributions
    lines_data = input_lines(data)
    qtext = query_text(query)
    graph_names = Set(nb.vars)
    q = parse_ldag_distribution(qtext, "target distribution")
    add_vars!(nb, unique(vcat(q.lhs, q.rhs)))
    parsed = LdagParsedDist[]
    for line in lines_data
        d = parse_ldag_distribution(line, "input distribution")
        add_vars!(nb, unique(vcat(d.lhs, d.rhs)))
        push!(parsed, d)
    end
    to_term(d) = Dist(varset(nb, d.lhs), varset(nb, d.rhs), varset(nb, d.zero), varset(nb, d.one))
    data_dists = [to_term(d) for d in parsed]
    for (t, d) in zip(data_dists, parsed)
        isempty(t.a ∩ t.b) || throw(ArgumentError("Invalid input distribution " * d.text *
                                                  ": same variable on the left and right-hand side."))
    end
    query_dist = to_term(q)
    isempty(query_dist.a ∩ query_dist.b) ||
        throw(ArgumentError("Invalid query " * q.text * ": same variable on the left and right-hand side."))
    length(nb.vars) <= 64 ||
        throw(ArgumentError("The inputs imply a graph with more than 64 nodes."))
    extra = String.(collect(context_vars))
    unknown = filter(!in(keys(nb.nums)), extra)
    isempty(unknown) || throw(ArgumentError("Context variables $(cs(unknown)) are not in the problem."))
    # A variable fixed to a value in the query or an input is a context
    # variable too, so that the search can assign values to it. The R package
    # needs it in `con_vars` (an intentional difference): without it,
    # p(Y|C=0,I_C=1) is not identified even from p(X,C,Y) when C is in no label.
    assigned = reduce(∪, (d.c ∪ d.d for d in data_dists); init = query_dist.c ∪ query_dist.d)
    con = varset(nb, con_names) ∪ varset(nb, extra) ∪ setdiff(assigned, varset(nb, intv_names))
    warn && warn_unconnected([v for d in vcat(q, parsed) for v in vcat(d.lhs, d.rhs)], graph_names)
    warn && check_case(nb)
    return LdagProblem(copy(nb.vars), edges, context_sets, local_csi, con, varset(nb, intv_names),
                       data_dists, query_dist, lines_data, qtext, String(graph))
end

function validate_label(from, to, pa, names, values)
    prefix = "Invalid label for edge $from -> $to: "
    from in names && throw(ArgumentError(prefix * "$from cannot appear in the label."))
    to in names && throw(ArgumentError(prefix * "$to cannot appear in the label."))
    allunique(names) || throw(ArgumentError(prefix * "duplicate assignment."))
    all(in(pa), names) || throw(ArgumentError(prefix * "only other parents of $to may be assigned."))
    # R ignores values other than 0 and 1; they are rejected here.
    all(in(("0", "1")), values) || throw(ArgumentError(prefix * "values must be 0 or 1."))
    return nothing
end

# Labels implied by the explicit ones (`infer_labels`): a partial assignment
# of the other parents is a label if the edge vanishes for all of its
# completions. The partial assignments are over the non-empty proper subsets
# of the parents, in the order of R's `expand.grid`.
function infer_labels!(inferred, nb, present, from, to, pa)
    n = length(pa)
    for j in 2:((1 << n) - 1)
        sub = [k for k in 1:n if grid_value(j, k) == 1]
        for a in 1:(1 << length(sub))
            value(k) = grid_value(a, something(findfirst(==(k), sub)))
            covered = all(present[r] for r in eachindex(present)
                          if all(grid_value(r, k) == value(k) for k in sub))
            if covered
                zero = pa[[k for k in sub if value(k) == 0]]
                one = pa[[k for k in sub if value(k) == 1]]
                push!(inferred, LabelRow(varset(nb, zero), varset(nb, one), nb.nums[from],
                                         nb.nums[to], varset(nb, pa)))
            end
        end
    end
    return inferred
end

# A distribution `p(y|z)` of an LDAG, in terms of variable names.
struct LdagParsedDist
    lhs::Vector{String}
    rhs::Vector{String}
    zero::Vector{String}
    one::Vector{String}
    text::String
end

# Parse a distribution `p(y)` or `p(y|z)`, where variables can be assigned
# the values 0 or 1, as in `p(Y|X,I_X=1)`.
function parse_ldag_distribution(d::AbstractString, kind::AbstractString)
    s = replace(d, r"\s+" => "")
    parse_error() = throw(ArgumentError("Unable to parse $kind: $d. " * LDAG_DISTRIBUTION_SYNTAX))
    m = match(r"^[Pp]\((.*)\)$", s)
    m === nothing && parse_error()
    inner = something(m.captures[1])   # the group always matches
    count('|', inner) > 1 && parse_error()
    bar = findfirst('|', inner)
    lhs = parse_names(bar === nothing ? inner : inner[1:prevind(inner, bar)])
    rhs = bar === nothing ? String[] : parse_names(inner[nextind(inner, bar):end])
    (lhs === nothing || rhs === nothing) && parse_error()
    zero = String[]
    one = String[]
    for part in (lhs, rhs)
        dups = unique(part[nonunique_mask(part)])
        isempty(dups) || throw(ArgumentError("Invalid $kind: $d, duplicated variables $(cs(dups))."))
        for i in findall(contains('='), part)
            e = split(part[i], '=')
            value = length(e) == 2 ? e[2] : ""
            value in ("0", "1") || throw(ArgumentError("Invalid value assignment in $kind: $d."))
            part[i] = String(e[1])
            push!(value == "0" ? zero : one, part[i])
        end
    end
    return LdagParsedDist(lhs, rhs, zero, one, String(d))
end

# Warning for variable names that differ only by case (`check_valid_input`).
function check_case(nb::Numbering)
    lower = any(v -> !isempty(v) && all(islowercase, v), nb.vars)
    upper = any(v -> !isempty(v) && all(isuppercase, v), nb.vars)
    lower && upper && @warn "Both lower case and upper case inputs detected."
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", p::LdagProblem)
    println(io, "LdagProblem with ", problem_nvars(p), " variables")
    println(io, "  query: ", term_string(p, p.query))
    print(io, "  data:  ")
    join(io, (term_string(p, d) for d in p.data), ", ")
end
