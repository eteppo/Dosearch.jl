# Candidates: additions to the causal model of a problem that the search for
# maximal identifiable graphs chooses from. Adding a candidate must never make
# a non-identifiable query identifiable, so that the identifiable sets of
# candidates are closed under taking subsets.

"""
    Candidate

Abstract type of the additions to a causal model that
[`maximal_identifiable_graphs`](@ref) chooses from, such as
[`BidirectedPair`](@ref).

Adding candidates to a model must only add paths to the graph, so that a
query that is identifiable with a set of candidates stays identifiable with
any subset of it. To add a kind of candidate, define
`Dosearch.confounded_pairs(problem::DagProblem, candidate)`, which returns the
bidirected edges that the candidate adds as pairs of variable indices, and
`Dosearch.latent_edges(problem::LdagProblem, candidate)`, which returns the
edges of the latent node that the candidate adds to an LDAG.
"""
abstract type Candidate end

"""
    BidirectedPair(x, y)

The candidate that adds the bidirected edge `x <-> y`, i.e. a latent common
cause of the variables named `x` and `y`. The pair is unordered. In an LDAG,
which has no bidirected edges, it adds a latent node `U_x_y` with the edges
`U_x_y -> x` and `U_x_y -> y`.
"""
struct BidirectedPair <: Candidate
    x::String
    y::String
    function BidirectedPair(x::AbstractString, y::AbstractString)
        x == y && throw(ArgumentError("a bidirected pair needs two distinct variables, got $x twice"))
        return x < y ? new(String(x), String(y)) : new(String(y), String(x))
    end
end

Base.show(io::IO, c::BidirectedPair) = print(io, "BidirectedPair(", repr(c.x), ", ", repr(c.y), ")")
Base.show(io::IO, ::MIME"text/plain", c::BidirectedPair) = print(io, c.x, " <-> ", c.y)

function variable_index(p::DagProblem, name::AbstractString)
    i = findfirst(==(name), p.labels)
    i === nothing && throw(ArgumentError("variable $name is not in the problem"))
    return i
end

proxies(p::DagProblem) = p.md_proxies

"""
    confounded_pairs(problem, candidate) -> Vector{Tuple{Int,Int}}

The bidirected edges, as pairs of variable indices of `problem`, that
`candidate` adds to the model.
"""
function confounded_pairs(p::DagProblem, c::BidirectedPair)
    x = variable_index(p, c.x)
    y = variable_index(p, c.y)
    for (i, name) in ((x, c.x), (y, c.y))
        i in p.transportability &&
            throw(ArgumentError("transportability node $name cannot have a latent common cause"))
        i in proxies(p) &&
            throw(ArgumentError("proxy variable $name cannot have a latent common cause"))
    end
    return [(x, y)]
end

"""
    LatentCause(variables...; labels = [])
    LatentCause(variables; labels = [])

The candidate that adds a latent common cause of all the named variables: a
bidirected edge between every pair of them, added together. `LatentCause("x",
"y")` adds the same edge as `BidirectedPair("x", "y")`. The order of the
variables does not matter.

In an LDAG, it adds a latent node with an edge into each of the variables.
`labels` gives some of these edges a label, as pairs `child => label` in the
syntax of the graph: `LatentCause("x", "y"; labels = ["y" => "a = 1"])` is a
latent cause of `x` and `y` whose edge into `y` vanishes where `a = 1`. The
variables of a label must be parents of the child. The search for maximal
graphs treats them as context variables throughout (see
[`maximal_identifiable_graphs`](@ref)). An edge into a variable with an
intervention node `I_x` also vanishes where `I_x = 1`, like the other edges
into an intervened variable. Labels are not supported for DAGs.
"""
struct LatentCause <: Candidate
    variables::Vector{String}
    labels::Vector{Pair{String,String}}
    function LatentCause(variables; labels = Pair{String,String}[])
        v = sort!(unique(String.(collect(variables))))
        length(v) >= 2 || throw(ArgumentError("a latent cause needs at least two distinct variables"))
        l = sort!([String(c) => canonical_label(lab) for (c, lab) in labels])
        for (c, lab) in l
            c in v || throw(ArgumentError("label for $c, which is not a child of the latent cause"))
            isempty(label_variables(lab)) && throw(ArgumentError("empty label for $c"))
        end
        allunique(first.(l)) || throw(ArgumentError("more than one label for a child of the latent cause"))
        return new(v, l)
    end
end

LatentCause(x::AbstractString, y::AbstractString, rest::AbstractString...; kwargs...) =
    LatentCause([x, y, rest...]; kwargs...)

# The variables assigned in an edge label `a = 1, b = 0; c = 1`.
label_variables(label::AbstractString) =
    unique([String(first(split(t, '='))) for t in split(replace(label, r"\s+" => ""), r"[;,]") if !isempty(t)])

# An edge label with uniform spacing: `a = 1, b = 0; c = 1`.
function canonical_label(label::AbstractString)
    rows = filter(!isempty, split(replace(label, r"\s+" => ""), ';'))
    return join([join([replace(t, "=" => " = ") for t in split(r, ',') if !isempty(t)], ", ") for r in rows], "; ")
end

Base.:(==)(a::LatentCause, b::LatentCause) = a.variables == b.variables && a.labels == b.labels
Base.hash(c::LatentCause, h::UInt) = hash(c.labels, hash(c.variables, hash(:LatentCause, h)))
function Base.show(io::IO, c::LatentCause)
    print(io, "LatentCause(", join(repr.(c.variables), ", "))
    isempty(c.labels) || print(io, "; labels = ", c.labels)
    print(io, ")")
end
function Base.show(io::IO, ::MIME"text/plain", c::LatentCause)
    print(io, "latent cause of ", join(c.variables, ", "))
    isempty(c.labels) || print(io, " (", join(("into $v : $l" for (v, l) in c.labels), "; "), ")")
end

function confounded_pairs(p::DagProblem, c::LatentCause)
    isempty(c.labels) || throw(ArgumentError("edge labels of latent causes need an LDAG"))
    idx = Int[]
    for name in c.variables
        i = variable_index(p, name)
        i in p.transportability &&
            throw(ArgumentError("transportability node $name cannot have a latent common cause"))
        i in proxies(p) &&
            throw(ArgumentError("proxy variable $name cannot have a latent common cause"))
        push!(idx, i)
    end
    return [(idx[a], idx[b]) for a in eachindex(idx) for b in (a + 1):length(idx)]
end

"""
    with_candidates(problem, candidates) -> DosearchProblem

The problem whose causal model has the `candidates` added.
"""
function with_candidates(p::DagProblem, candidates)
    extra = Tuple{Int,Int}[]
    for c in candidates
        append!(extra, confounded_pairs(p, c))
    end
    return DagProblem(p.labels, p.directed, vcat(p.bidirected, extra), p.transportability,
                           p.selection_bias, p.md_switches, p.md_proxies, p.data, p.query,
                           p.md_sym, p.data_text, p.query_text)
end

"""
    confounding_candidates(problem) -> Vector{BidirectedPair}

Every bidirected edge that can be added to the model of `problem`: all pairs
of variables that are not already bidirected, except pairs with a
transportability node or a missing data proxy. The pairs are ordered by the
internal variable numbering.
"""
function confounding_candidates(p::DagProblem)
    n = problem_nvars(p)
    excluded = p.transportability ∪ proxies(p)
    present = Set(minmax(x, y) for (x, y) in p.bidirected)
    out = BidirectedPair[]
    for x in 1:n, y in (x + 1):n
        (x in excluded || y in excluded || (x, y) in present) && continue
        push!(out, BidirectedPair(p.labels[x], p.labels[y]))
    end
    return out
end

# LDAGs -----------------------------------------------------------------------

function variable_index(p::LdagProblem, name::AbstractString)
    i = findfirst(==(name), p.labels)
    i === nothing && throw(ArgumentError("variable $name is not in the problem"))
    i in p.intervention_vars &&
        throw(ArgumentError("intervention node $name cannot have a latent common cause"))
    return i
end

"""
    latent_edges(problem::LdagProblem, candidate) -> Vector{Pair{String,String}}

The edges from the latent node that `candidate` adds to the LDAG of
`problem`, as `child => label` pairs (an empty label for an unlabeled edge).
[`with_candidates`](@ref) adds the intervention contexts of the children to
the labels.
"""
latent_edges(p::LdagProblem, c::BidirectedPair) = latent_edges(p, LatentCause(c.x, c.y))

function latent_edges(p::LdagProblem, c::LatentCause)
    labels = Dict(c.labels)
    edges = Pair{String,String}[]
    for v in c.variables
        y = variable_index(p, v)
        label = get(labels, v, "")
        for a in label_variables(label)
            i = findfirst(==(a), p.labels)
            (i !== nothing && (i, y) in p.edges) ||
                throw(ArgumentError("invalid label for the latent cause of $v: $a is not a parent of $v"))
        end
        push!(edges, v => label)
    end
    return edges
end

# The name of a latent node: `U_` and its children, then each labeled child
# with its label (`U_x_y_y_a1` for the label `a = 1` on the edge into `y`),
# primed until it is unused. Equal names stand for nodes with equal edges,
# whichever problem they are in.
function latent_name(edges::Vector{Pair{String,String}}, used)
    name = "U_" * join(first.(edges), "_")
    for (v, l) in edges
        isempty(l) || (name *= "_" * v * "_" * replace(l, " " => "", "=" => "", "," => "", ";" => "or"))
    end
    while name in used
        name *= "'"
    end
    return name
end

function with_candidates(p::LdagProblem, candidates)
    used = Set(p.labels)
    lines = String[]
    for c in candidates
        edges = latent_edges(p, c)
        u = latent_name(edges, used)
        push!(used, u)
        for (v, l) in edges
            # the edge vanishes where the child is intervened on
            y = variable_index(p, v)
            rows = vcat(isempty(l) ? String[] : [l],
                        [p.labels[i] * " = 1" for i in p.intervention_vars if (i, y) in p.edges])
            push!(lines, isempty(rows) ? u * " -> " * v : u * " -> " * v * " : " * join(rows, "; "))
        end
    end
    isempty(lines) && return p
    graph = p.graph_text * "\n" * join(lines, "\n")
    return LdagProblem(p.data_text, p.query_text, graph; context_vars = p.labels[collect(p.context_vars)],
                       warn = false, order = p.labels)
end

# The problem with the variables in the labels of `candidates` as context
# variables, so that every set of candidates is searched with the same
# context variables. (More context variables allow more case distinctions,
# which could make a larger set of candidates identifiable when a smaller one
# is not.)
with_label_context(p::DagProblem, candidates) = p

function with_label_context(p::LdagProblem, candidates)
    names = String[]
    for c in candidates
        c isa LatentCause || continue
        for (_, label) in c.labels
            append!(names, filter(v -> !startswith(v, "I_"), label_variables(label)))
        end
    end
    extra = filter(v -> (i = findfirst(==(v), p.labels); i === nothing || !(i in p.context_vars)), unique(names))
    isempty(extra) && return p
    foreach(v -> variable_index(p, v), extra)
    return LdagProblem(p.data_text, p.query_text, p.graph_text;
                       context_vars = vcat(p.labels[collect(p.context_vars)], extra), warn = false, order = p.labels)
end

"""
    confounding_candidates(problem::LdagProblem) -> Vector{BidirectedPair}

Every pair of variables of the LDAG of `problem` that can get a latent common
cause, i.e. all pairs without intervention nodes, ordered by the internal
variable numbering.
"""
function confounding_candidates(p::LdagProblem)
    vars = [i for i in 1:problem_nvars(p) if !(i in p.intervention_vars)]
    return [BidirectedPair(p.labels[x], p.labels[y]) for (k, x) in enumerate(vars) for y in vars[(k + 1):end]]
end
