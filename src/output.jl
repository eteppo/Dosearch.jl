# Output of the search: terms as text and the derivation as a DOT graph. Port
# of `to_string`, `dec_to_text` and `draw` in the R package's C++ code. The
# strings must match the R package character for character. The identifying
# formula is in formula.jl.

# Comma separated labels of the variables in `s`; enabled missing data
# mechanisms are shown as `r_x = 1` (`dec_to_text`).
function set_string(io::IO, p::DagProblem, s::VarSet, enabled::VarSet)
    first = true
    for i in s
        first || print(io, ",")
        first = false
        print(io, p.labels[i])
        i in enabled && print(io, " = ", p.md_sym)
    end
    return nothing
end

"""
    term_string(problem, p::Dist) -> String

The term `p` in the notation of the R package, e.g. `p(y|do(x),z)`.
"""
function term_string(p::DagProblem, t::Dist)
    io = IOBuffer()
    print(io, "p(")
    set_string(io, p, t.a, t.d)
    isempty(t.b) || print(io, "|")
    if !isempty(t.c)
        print(io, "do(")
        set_string(io, p, t.c, t.d)
        print(io, ")")
    end
    w = setdiff(t.b, t.c)
    if !isempty(w)
        isempty(t.c) || print(io, ",")
        set_string(io, p, w, t.d)
    end
    print(io, ")")
    return String(take!(io))
end

function set_string(p::DagProblem, s::VarSet)
    io = IOBuffer()
    set_string(io, p, s, VarSet())
    return String(take!(io))
end

# DOT output of derivations (`derivation` class in derivation.cpp).
struct DotWriter
    io::IOBuffer
    ids::Dict{String,Int}
end
DotWriter() = DotWriter(IOBuffer(), Dict{String,Int}())

function node_id!(w::DotWriter, label::String)
    id = get(w.ids, label, 0)
    if id == 0
        id = length(w.ids) + 1
        w.ids[label] = id
        print(w.io, "n", id, "[shape=polygon,sides=4,label=\"", label, "\"]\n")
    end
    return id
end

function add_edge!(w::DotWriter, from::String, to::String, label::String)
    # The C++ code builds the edge line with `get_label(from) + ... +
    # get_label(to)`, whose evaluation order is unspecified; GCC evaluates
    # the `to` label first, which decides the node numbering.
    t = node_id!(w, to)
    f = node_id!(w, from)
    print(w.io, "n", f, " -> n", t, "[label=\"", label, "\"]\n")
    return w
end

# `search::draw`. Shared sub-derivations are drawn again each time they are
# reached, as in the C++ code; the graph is `strict`, so the duplicate edges
# have no effect on the drawing.
function draw!(w::DotWriter, st::SearchState, node::Node, recursive::Bool)
    node.parent1 > 0 || return w
    p = problem(st)
    to = term_string(p, node.dist)
    label = rule_label(node.rule)
    pa1 = st.nodes[node.parent1]
    add_edge!(w, term_string(p, pa1.dist), to, label)
    recursive && draw!(w, st, pa1, true)
    if node.parent2 > 0
        pa2 = st.nodes[node.parent2]
        add_edge!(w, term_string(p, pa2.dist), to, label)
        recursive && draw!(w, st, pa2, true)
    end
    return w
end

"""
    derivation_dot(st; draw_all = false) -> String

The derivation of the query as a DOT graph. With `draw_all = true`, every
derivation step of the search is drawn, not only those leading to the query.
The R package draws the steps of `draw_all` in the iteration order of a hash
table; here they are drawn in the order they were found.
"""
function derivation_dot(st::SearchState; draw_all::Bool = false)
    w = DotWriter()
    print(w.io, "strict digraph InferenceTree {\n")
    if draw_all
        for node in st.nodes
            draw!(w, st, node, false)
        end
        st.target === nothing || draw!(w, st, st.target, false)
    elseif st.target !== nothing
        draw!(w, st, st.target, true)
    end
    print(w.io, "}\n")
    return String(take!(w.io))
end
