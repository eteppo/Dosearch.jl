# Graphs and distributions as Julia values, for building problems in code.
# They render to the text syntax, which the problem constructors parse and
# validate, so a value and its text give the same problem, with the same
# variable numbering. The values check what the text cannot: that names are
# valid, so that the rendered text means exactly the value.

const INPUT_NAME_FORBIDDEN = ('(', ')', ',', '|', '=', ':', ';', '$')

# A variable name of a typed input: a string or symbol that the text syntax
# reads back as the same name.
function input_name(x::Union{AbstractString,Symbol})
    s = String(x)
    isempty(s) && throw(ArgumentError("A variable name cannot be empty."))
    if any(isspace, s) || any(in(INPUT_NAME_FORBIDDEN), s) || occursin("->", s) || occursin("--", s)
        throw(ArgumentError("Invalid variable name $(repr(s)): names cannot contain whitespace, " *
                            "`(`, `)`, `,`, `|`, `=`, `:`, `;`, `\$`, `->` or `--`."))
    end
    return s
end
input_name(x) = throw(ArgumentError("A variable name must be a string or a symbol, got $(repr(x))."))

# The value of an assignment: 0 or 1, or the symbol of an enabled missing data
# mechanism (`md_sym`).
function assigned_value(v::Integer)
    v == 0 || v == 1 || throw(ArgumentError("An assigned value must be 0 or 1, got $v."))
    return string(Int(v))
end
assigned_value(v::Union{AbstractString,Symbol}) = input_name(v)
assigned_value(v) = throw(ArgumentError("An assigned value must be 0, 1 or a string, got $(repr(v))."))

# A single item or a collection of them.
items(x::Union{AbstractString,Symbol,Pair}) = (x,)
items(xs) = xs

"""
    LabeledEdge(from => to, assignments...)

An edge of an LDAG that vanishes in a context: the `assignments`, pairs
`variable => value` with the values 0 or 1, of other parents of `to`.
`LabeledEdge("X" => "Y", "A" => 1)` is the line `X -> Y : A = 1` of a graph.
An edge with several labeled edges vanishes in each of their contexts.
"""
struct LabeledEdge
    from::String
    to::String
    context::Vector{Pair{String,Int}}
    function LabeledEdge(from::String, to::String, context::Vector{Pair{String,Int}})
        from == to && throw(ArgumentError("Invalid edge $from -> $to: a self-loop."))
        isempty(context) && throw(ArgumentError("The labeled edge $from -> $to needs at least one assignment; " *
                                                "an edge without a label is the pair $(repr(from)) => $(repr(to))."))
        names = first.(context)
        allunique(names) || throw(ArgumentError("The label of $from -> $to assigns a variable twice."))
        return new(from, to, context)
    end
end

function LabeledEdge(edge::Pair, assignments::Pair...)
    function value(v)
        v isa Integer && (v == 0 || v == 1) && return Int(v)
        throw(ArgumentError("A label assigns the value 0 or 1, got $(repr(v))."))
    end
    context = Pair{String,Int}[input_name(first(a)) => value(last(a)) for a in assignments]
    return LabeledEdge(input_name(first(edge)), input_name(last(edge)), context)
end

Base.:(==)(a::LabeledEdge, b::LabeledEdge) = a.from == b.from && a.to == b.to && a.context == b.context
Base.hash(e::LabeledEdge, h::UInt) = hash((e.from, e.to, e.context), hash(:LabeledEdge, h))
Base.show(io::IO, e::LabeledEdge) =
    print(io, "LabeledEdge(", repr(e.from), " => ", repr(e.to), join((", $(repr(k)) => $v" for (k, v) in e.context)), ")")
Base.show(io::IO, ::MIME"text/plain", e::LabeledEdge) = print(io, edge_text(e))

const GraphEdge = Union{Pair{String,String},BidirectedPair,LabeledEdge}

"""
    CausalGraph(edges...)

A causal graph as a value: the same graph as the text syntax of
[`dosearch`](@ref), for building and combining graphs in code. The edges are

- `from => to`, a directed edge, with names as strings or symbols;
- [`BidirectedPair`](@ref)`(x, y)`, an unobserved common cause of `x` and `y`;
- [`LabeledEdge`](@ref)`(from => to, assignments...)`, an edge of an LDAG that
  vanishes in a context;
- a graph in the text syntax, or another `CausalGraph`, whose edges are
  included in order;
- a collection or generator of edges.

Repeated edges are kept once. The order of the edges is kept: it fixes the
numbering of the variables, and so the order in which formulas list them.
Two graphs are equal if they have the same edges in the same order.
`string(g)` is the graph in the text syntax, and `g ∪ h` contains the edges
of both. Wherever a problem takes a graph, a `CausalGraph` can be given.

# Examples

```jldoctest
julia> using Dosearch: CausalGraph

julia> g = CausalGraph("x" => "y", "z" => "x", "z" => "y")
CausalGraph with 3 edges:
  x -> y
  z -> x
  z -> y

julia> g ∪ CausalGraph(BidirectedPair("x", "y"))
CausalGraph with 4 edges:
  x -> y
  z -> x
  z -> y
  x <-> y

julia> chain = CausalGraph("A\$k" => "A\$(k + 1)" for k in 1:3);

julia> print(chain)
A1 -> A2
A2 -> A3
A3 -> A4
```
"""
struct CausalGraph
    edges::Vector{GraphEdge}
    function CausalGraph(edges::Vector{GraphEdge})
        edges = unique(edges)
        directed = Tuple{String,String}[]
        bidirected = labeled = false
        for e in edges
            if e isa BidirectedPair
                bidirected = true
            else
                from, to = e isa LabeledEdge ? (e.from, e.to) : (first(e), last(e))
                from == to && throw(ArgumentError("Invalid edge $from -> $to: a self-loop."))
                push!(directed, (from, to))
                labeled |= e isa LabeledEdge
            end
        end
        bidirected && labeled && throw(ArgumentError("A graph with labeled edges is an LDAG, which has no " *
                                                     "bidirected edges: add a latent node with an edge into each " *
                                                     "of the variables instead."))
        is_acyclic(directed) || throw(ArgumentError("Invalid graph, the graph contains cycles."))
        return new(edges)
    end
end

CausalGraph(items...) = CausalGraph(add_edges!(GraphEdge[], items))

function add_edges!(edges::Vector{GraphEdge}, items)
    for item in items
        add_edge!(edges, item)
    end
    return edges
end
add_edge!(edges::Vector{GraphEdge}, e::Pair) = push!(edges, input_name(first(e)) => input_name(last(e)))
add_edge!(edges::Vector{GraphEdge}, e::Union{BidirectedPair,LabeledEdge}) = push!(edges, e)
add_edge!(edges::Vector{GraphEdge}, g::CausalGraph) = append!(edges, g.edges)
add_edge!(edges::Vector{GraphEdge}, text::AbstractString) = append!(edges, parse_edges(text))
add_edge!(edges::Vector{GraphEdge}, items) = add_edges!(edges, items)

# The edges of a graph in the text syntax.
function parse_edges(text::AbstractString)
    check_line_breaks(text, "graph", "edge")
    edges = GraphEdge[]
    for line in split(text, r"\r|\n")
        l = strip(line)
        isempty(l) && continue
        edge, label = occursin(':', l) ? split(l, ':'; limit = 2) : (l, nothing)
        pg = parse_graph(edge)
        if label === nothing
            append!(edges, [from => to for (from, to) in pg.directed])
            append!(edges, [BidirectedPair(x, y) for (x, y) in pg.bidirected])
            continue
        end
        length(pg.directed) == 1 && isempty(pg.bidirected) ||
            throw(ArgumentError("Invalid graph, malformed line: $l. " * LDAG_SYNTAX))
        from, to = only(pg.directed)
        for alternative in split(label, ';')
            assignments = Pair{String,Int}[]
            for a in split(alternative, ',')
                parts = strip.(split(a, '='))
                length(parts) == 2 && parts[2] in ("0", "1") ||
                    throw(ArgumentError("Invalid label in the line $l: a label assigns the values 0 or 1, " *
                                        "as in `z = 1`."))
                push!(assignments, input_name(parts[1]) => parse(Int, parts[2]))
            end
            push!(edges, LabeledEdge(from, to, assignments))
        end
    end
    return edges
end

edge_text(e::Pair{String,String}) = first(e) * " -> " * last(e)
edge_text(e::BidirectedPair) = e.x * " <-> " * e.y
edge_text(e::LabeledEdge) = e.from * " -> " * e.to * " : " * join(("$k = $v" for (k, v) in e.context), ", ")

Base.print(io::IO, g::CausalGraph) = foreach(e -> println(io, edge_text(e)), g.edges)
function Base.show(io::IO, ::MIME"text/plain", g::CausalGraph)
    print(io, "CausalGraph with ", length(g.edges), length(g.edges) == 1 ? " edge" : " edges", isempty(g.edges) ? "" : ":")
    foreach(e -> print(io, "\n  ", edge_text(e)), g.edges)
end
Base.show(io::IO, g::CausalGraph) = print(io, "CausalGraph(", join((repr(e) for e in g.edges), ", "), ")")
Base.:(==)(a::CausalGraph, b::CausalGraph) = a.edges == b.edges
Base.hash(g::CausalGraph, h::UInt) = hash(g.edges, hash(:CausalGraph, h))
Base.union(g::CausalGraph, hs::CausalGraph...) = CausalGraph(g, hs...)

"""
    P(variables...; given = (), intervened = ())

A distribution as a value: `P("y"; given = "z", intervened = "x")` is
`p(y|do(x),z)` in the text syntax of [`dosearch`](@ref), for building inputs
and queries in code. `given` and `intervened` take a name or a collection of
names, as strings or symbols. A variable can be assigned a value by a pair
`name => value`: 0 or 1 in an LDAG, as in `P("Y"; given = ["X", "I_X" => 1])`,
and the symbol of an enabled missing data mechanism (`md_sym`, 1 by default)
in a DAG, as in `P("x*"; given = "r_x" => 1)`.

`string(d)` is the distribution in the text syntax. Wherever a problem takes
its input distributions or its query, a `P` or a vector of them can be given.

# Examples

```jldoctest
julia> using Dosearch: CausalGraph, P

julia> P("y"; given = "z", intervened = "x")
p(y | do(x), z)

julia> dosearch([P("x", "y", "z")], P("y"; intervened = "x"), CausalGraph("x" => "y", "z" => "x", "z" => "y"))
\\sum_{z}\\left(p(z)p(y|x,z)\\right)
```
"""
struct P
    variables::Vector{String}
    given::Vector{String}
    intervened::Vector{String}
    values::Vector{Pair{String,String}}
    function P(variables::Vector{String}, given::Vector{String}, intervened::Vector{String},
               values::Vector{Pair{String,String}})
        isempty(variables) && throw(ArgumentError("A distribution needs at least one variable on the left."))
        names = vcat(variables, given, intervened)
        allunique(names) || throw(ArgumentError("A variable appears twice in the distribution " *
                                                "p($(join(variables, ", ")) | ...)."))
        return new(variables, given, intervened, values)
    end
end

function P(variables...; given = (), intervened = ())
    if length(variables) == 1 && !(only(variables) isa Union{AbstractString,Symbol,Pair})
        variables = Tuple(only(variables))   # a collection of variables
    end
    values = Pair{String,String}[]
    function name!(item)
        item isa Pair || return input_name(item)
        n = input_name(first(item))
        push!(values, n => assigned_value(last(item)))
        return n
    end
    return P(String[name!(v) for v in variables], String[name!(v) for v in items(given)],
             String[name!(v) for v in items(intervened)], values)
end

function Base.print(io::IO, d::P)
    values = Dict(d.values)
    item(v) = haskey(values, v) ? v * " = " * values[v] : v
    print(io, "p(", join(item.(d.variables), ", "))
    rhs = String[]
    isempty(d.intervened) || push!(rhs, "do(" * join(item.(d.intervened), ", ") * ")")
    append!(rhs, item.(d.given))
    isempty(rhs) || print(io, " | ", join(rhs, ", "))
    print(io, ")")
end
Base.show(io::IO, ::MIME"text/plain", d::P) = print(io, d)
function Base.show(io::IO, d::P)
    values = Dict(d.values)
    item(v) = haskey(values, v) ? repr(v) * " => " * repr(something(tryparse(Int, values[v]), values[v])) : repr(v)
    list(vs) = length(vs) == 1 ? item(only(vs)) : "[" * join(item.(vs), ", ") * "]"
    print(io, "P(", join(item.(d.variables), ", "))
    isempty(d.given) || print(io, "; given = ", list(d.given))
    isempty(d.intervened) || print(io, isempty(d.given) ? "; " : ", ", "intervened = ", list(d.intervened))
    print(io, ")")
end
Base.:(==)(a::P, b::P) = a.variables == b.variables && a.given == b.given &&
                         a.intervened == b.intervened && Set(a.values) == Set(b.values)
Base.hash(d::P, h::UInt) = hash((d.variables, d.given, d.intervened, Set(d.values)), hash(:P, h))

# The problems take the values wherever they take text.
data_lines(d::P) = data_lines(string(d))
query_text(q::P) = string(q)

# Special nodes as a comma separated list, or as a collection of names; missing
# data mechanisms as `"r_x : x, r_y : y"`, or as pairs `"r_x" => "x"`.
node_list(x::Union{Nothing,AbstractString}) = x
node_list(xs) = join((input_name(x) for x in items(xs)), ", ")
mechanism_list(x::Union{Nothing,AbstractString}) = x
mechanism_list(xs) = join((input_name(first(m)) * " : " * input_name(last(m)) for m in items(xs)), ", ")

for f in (:DosearchProblem, :DagProblem, :LdagProblem, :dosearch)
    @eval $f(data, query, graph::CausalGraph; kwargs...) = $f(data, query, string(graph); kwargs...)
end
