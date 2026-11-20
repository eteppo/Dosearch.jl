# The identification problem: numbered variables, the graph, special nodes,
# input distributions and the query, validated. Port of the model building in
# `get_derivation_dag` and its helpers in `internal_dag.R`.
#
# The numbering of the variables follows the R package exactly, because it
# fixes the order of the subsets that the search tries and the order in which
# variables are printed in formulas:
#
# 1. Graph variables in order of first appearance in the list of children of
#    directed edges, then their parents, then the bidirected endpoints.
# 2. With missing data, the triples (true variable, response indicator, proxy)
#    of each mechanism come first, in the order given.
# 3. Transportability nodes and then selection bias nodes are moved last.
# 4. Variables that appear only in the query or the data are appended in order
#    of appearance, the query first.

"""
    DosearchProblem

Abstract type of validated identification problems: [`DagProblem`](@ref)
and [`LdagProblem`](@ref). `DosearchProblem(data, query, graph; kwargs...)`
builds an `LdagProblem` if the graph has edge labels and a `DagProblem`
otherwise. See [`dosearch`](@ref) for the syntax of the arguments.
"""
abstract type DosearchProblem end

function DosearchProblem(data, query, graph::AbstractString;
                         transportability = nothing, selection_bias = nothing,
                         missing_data = nothing, md_sym::AbstractString = "1",
                         context_vars = String[], warn::Bool = true)
    if occursin(':', graph)
        if transportability !== nothing || selection_bias !== nothing || missing_data !== nothing
            throw(ArgumentError("Transportability, selection bias and missing data are not supported for LDAGs."))
        end
        return LdagProblem(data, query, graph; context_vars, warn)
    end
    isempty(context_vars) || throw(ArgumentError("Context variables are only supported for LDAGs."))
    return DagProblem(data, query, graph; transportability, selection_bias, missing_data, md_sym, warn)
end

"""
    DagProblem(data, query, graph; transportability, selection_bias,
               missing_data, md_sym = "1", warn = true)

A validated causal effect identification problem in a DAG with latent
confounders: identify the `query` from the input distributions `data` in the
causal model `graph`. See [`dosearch`](@ref) for the syntax of the arguments.

The fields hold the problem in internal form: `labels` are the variable names
by index, and `data` and `query` are [`Dist`](@ref) terms over those indices.
"""
struct DagProblem <: DosearchProblem
    labels::Vector{String}
    directed::Vector{Tuple{Int,Int}}
    bidirected::Vector{Tuple{Int,Int}}
    transportability::VarSet
    selection_bias::VarSet
    md_switches::VarSet
    md_proxies::VarSet
    data::Vector{Dist}
    query::Dist
    md_sym::String
    data_text::Vector{String}
    query_text::String
end

problem_nvars(p::DagProblem) = length(p.labels)
has_missing_data(p::DagProblem) = !isempty(p.md_switches)

function rule_context(p::DagProblem)
    return RuleContext(problem_nvars(p), has_missing_data(p), p.md_switches, p.md_proxies,
                       p.md_proxies >> 2, p.transportability ∪ p.selection_bias)
end

"""
    causal_graph(problem) -> DsepGraph

The d-separation graph of `problem`, with intervention nodes.
"""
function causal_graph(p::DagProblem)
    g = DsepGraph(problem_nvars(p))
    for (from, to) in p.directed
        add_directed!(g, from, to)
    end
    for (x, y) in p.bidirected
        add_bidirected!(g, x, y)
    end
    return g
end

# Mutable variable numbering used while building a problem.
mutable struct Numbering
    vars::Vector{String}
    nums::Dict{String,Int}
end
Numbering(vars::Vector{String}) = Numbering(vars, Dict(v => i for (i, v) in enumerate(vars)))
function renumber!(nb::Numbering, vars::Vector{String})
    nb.vars = vars
    nb.nums = Dict(v => i for (i, v) in enumerate(vars))
    return nb
end
function add_vars!(nb::Numbering, names)
    for v in names
        if !haskey(nb.nums, v)
            push!(nb.vars, v)
            nb.nums[v] = length(nb.vars)
        end
    end
    return nb
end
varset(nb::Numbering, names) = VarSet([nb.nums[v] for v in names])

# Text inputs: one distribution per line, or a vector of them.
function data_lines(data::AbstractString)
    check_line_breaks(data, "data", "input distribution")
    return filter(!isempty, [replace(l, r"\s+" => "") for l in split(data, r"\r|\n")])
end
data_lines(data::NamedTuple) = data_lines(roles_to_distribution(data))
data_lines(data::AbstractVector) = reduce(vcat, [data_lines(d) for d in data]; init = String[])
function query_text(query::AbstractString)
    check_line_breaks(query, "query", "part")
    isempty(strip(query)) && throw(ArgumentError("The query is empty. " * DISTRIBUTION_SYNTAX))
    return String(strip(query))
end
query_text(query::NamedTuple) = roles_to_distribution(query)

# The input distributions, which must not be empty.
function input_lines(data)
    lines = data_lines(data)
    isempty(lines) && throw(ArgumentError("No input distributions: `data` lists the known distributions, " *
                                          "one per line. " * DISTRIBUTION_SYNTAX))
    return lines
end

# Warn about variables of the inputs that are not in the graph: they are
# unconnected, which is probably a typo.
function warn_unconnected(names, graph_names)
    extra = filter(!in(graph_names), unique(names))
    isempty(extra) && return nothing
    @warn "Variables $(cs(extra)) of the distributions are not in the graph, so they are not " *
          "connected to any other variable."
    return nothing
end

function DagProblem(data, query, graph::AbstractString; transportability = nothing,
                    selection_bias = nothing, missing_data = nothing,
                    md_sym::AbstractString = "1", warn::Bool = true)
    transportability = node_list(transportability)
    selection_bias = node_list(selection_bias)
    missing_data = mechanism_list(missing_data)
    if occursin(':', graph)
        throw(ArgumentError("Edge labels are only supported in LDAGs; use `LdagProblem` or `DosearchProblem`."))
    end
    check_line_breaks(graph, "graph", "edge")
    if isempty(strip(graph)) && missing_data === nothing
        throw(ArgumentError("Invalid graph, the graph is empty. " * GRAPH_SYNTAX))
    end
    lines = input_lines(data)
    qtext = query_text(query)
    pg = parse_graph(graph)
    dir_lhs = first.(pg.directed)
    dir_rhs = last.(pg.directed)
    bi_lhs = first.(pg.bidirected)
    bi_rhs = last.(pg.bidirected)
    nb = Numbering(unique(vcat(dir_rhs, dir_lhs, bi_rhs, bi_lhs)))

    md_s = VarSet()
    md_p = VarSet()
    if missing_data !== nothing
        mechanisms = parse_missing_data(missing_data)
        switches = first.(mechanisms)
        trues = last.(mechanisms)
        proxies = trues .* "*"
        parents_of_true = Set(dir_lhs[[r in trues for r in dir_rhs]])
        if any(in(parents_of_true), switches)
            throw(ArgumentError("A missing data mechanism cannot be a parent of a true variable."))
        end
        append!(dir_lhs, trues, switches)
        append!(dir_rhs, proxies, proxies)
        vars_md = String[]
        for (t, s, p) in zip(trues, switches, proxies)
            push!(vars_md, t, s, p)
        end
        renumber!(nb, vcat(vars_md, filter(!in(Set(vars_md)), nb.vars)))
        md_s = varset(nb, switches)
        md_p = varset(nb, proxies)
    end

    tr_vars = String[]
    if transportability !== nothing
        tr_vars = parse_node_list(transportability)
        absent = filter(!in(keys(nb.nums)), tr_vars)
        if !isempty(absent)
            throw(ArgumentError("Transportability nodes $(cs(absent)) are not present in the graph."))
        end
        if any(in(Set(vcat(dir_rhs, bi_rhs, bi_lhs))), tr_vars)
            throw(ArgumentError("Invalid graph: a transportability node cannot be a child of another node."))
        end
    end
    sb_vars = String[]
    if selection_bias !== nothing
        sb_vars = parse_node_list(selection_bias)
        absent = filter(!in(keys(nb.nums)), sb_vars)
        if !isempty(absent)
            throw(ArgumentError("Selection bias nodes $(cs(absent)) are not present in the graph."))
        end
        if any(in(Set(dir_lhs)), sb_vars)
            throw(ArgumentError("Invalid graph: a selection bias node cannot be a parent of another node."))
        end
    end
    tr = VarSet()
    sb = VarSet()
    if !isempty(tr_vars) || !isempty(sb_vars)
        special = Set(vcat(tr_vars, sb_vars))
        renumber!(nb, vcat(filter(!in(special), nb.vars), tr_vars, sb_vars))
        tr = varset(nb, tr_vars)
        sb = varset(nb, sb_vars)
    end

    has_md = missing_data !== nothing
    graph_names = Set(nb.vars)
    q = parse_distribution(qtext, "target distribution"; missing_data = has_md, md_sym)
    add_parsed_vars!(nb, q)
    parsed = ParsedDist[]
    for line in lines
        d = parse_distribution(line, "input distribution"; missing_data = has_md, md_sym)
        add_parsed_vars!(nb, d)
        push!(parsed, d)
    end
    md_t = md_p >> 2
    data_dists = [to_dist(nb, d) for d in parsed]
    for (d, pd) in zip(data_dists, parsed)
        validate_distribution(d, "Invalid input distribution ", pd, md_s, md_p, md_t, tr, sb)
    end
    query_dist = to_dist(nb, q)
    validate_distribution(query_dist, "Invalid query ", q, md_s, md_p, md_t, tr, sb)
    n = length(nb.vars)
    if n > 32
        throw(ArgumentError("The inputs imply a graph with more than 64 nodes " *
                            "(at most 32 variables and their intervention nodes are supported)."))
    end
    if warn
        warn_unconnected([v for d in vcat(q, parsed) for v in vcat(d.lhs, d.conditioning, d.intervened)],
                         graph_names)
        check_inputs(nb, parsed, md_s)
    end
    directed = [(nb.nums[a], nb.nums[b]) for (a, b) in zip(dir_lhs, dir_rhs)]
    bidirected = [(nb.nums[a], nb.nums[b]) for (a, b) in zip(bi_lhs, bi_rhs)]
    return DagProblem(copy(nb.vars), directed, bidirected, tr, sb, md_s, md_p,
                           data_dists, query_dist, String(md_sym), lines, qtext)
end

add_parsed_vars!(nb::Numbering, d::ParsedDist) =
    add_vars!(nb, unique(vcat(d.lhs, d.conditioning, d.intervened)))

function to_dist(nb::Numbering, d::ParsedDist)
    if !isempty(d.conditioning ∩ d.intervened)
        throw(ArgumentError("Invalid distribution $(d.text): variables " *
                            "$(cs(d.conditioning ∩ d.intervened)) are both conditioned on and intervened on."))
    end
    c = varset(nb, d.intervened)
    return Dist(varset(nb, d.lhs), varset(nb, d.conditioning) ∪ c, c, varset(nb, d.enabled))
end

# Semantic checks of a distribution (`validate_distribution_dag`).
function validate_distribution(p::Dist, msg::String, d::ParsedDist, md_s::VarSet, md_p::VarSet,
                               md_t::VarSet, tr::VarSet, sb::VarSet)
    fail(reason, sep = ": ") = throw(ArgumentError(msg * d.text * sep * reason))
    isempty(((p.a ∩ md_p) >> 2) ∩ (p.b ∩ md_t)) ||
        fail("proxy variable of a true variable on the left-hand side.")
    isempty(((p.a ∩ md_t) << 2) ∩ (p.b ∩ md_p)) ||
        fail("true variable of a proxy variable on the left-hand side.")
    isempty(((p.a ∩ md_p) >> 2) ∩ (p.a ∩ md_t)) ||
        fail("true and proxy versions of the same variable on the left-hand side.")
    isempty(((p.b ∩ md_p) >> 2) ∩ (p.b ∩ md_t)) ||
        fail("true and proxy versions of the same variable on the right-hand side.")
    isempty(p.a ∩ p.b) || fail("same variable on the left and right-hand side.")
    isempty(p.a ∩ tr) || fail("transportability node on the left-hand side.")
    isempty(p.a ∩ sb) || fail("selection bias node on the left-hand side.")
    isempty(p.c ∩ tr) || fail("intervention on a transportability node.", ":\n")
    isempty(p.c ∩ sb) || fail("intervention on a selection bias node.")
    p.d ⊆ md_s || fail("value assignment of a non-missing data mechanism.")
    return nothing
end

# Warnings for inputs that are valid but likely mistakes (`check_valid_input`).
function check_inputs(nb::Numbering, parsed::Vector{ParsedDist}, md_s::VarSet)
    pool = varset(nb, unique(reduce(vcat, [d.lhs for d in parsed]; init = String[])))
    absent = setdiff(md_s, pool)
    if !isempty(absent)
        @warn "There are response indicators that are not present in any input distribution: " *
              join(nb.vars[collect(absent)], ", ")
    end
    lower = any(v -> !isempty(v) && all(islowercase, v), nb.vars)
    upper = any(v -> !isempty(v) && all(isuppercase, v), nb.vars)
    if lower && upper
        @warn "Both lower case and upper case inputs detected."
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", p::DagProblem)
    println(io, "DagProblem with ", problem_nvars(p), " variables")
    println(io, "  query: ", term_string(p, p.query))
    print(io, "  data:  ")
    join(io, (term_string(p, d) for d in p.data), ", ")
end
