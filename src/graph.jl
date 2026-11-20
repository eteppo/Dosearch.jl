# d-separation in a causal graph with latent confounders and intervention
# nodes, for the do-calculus rules. Replaces `dcongraph` in the R package's
# `dcongraph.cpp`.
#
# The C++ code tests one pair of nodes at a time: it intervenes on, conditions
# on and marginalizes out nodes one by one, keeping three adjacency matrices,
# and reports a dependence when the pair becomes adjacent. Here a statement
# about two sets of nodes is decided by one reachability search over active
# paths. The two agree on every statement that the rules make, where the
# intervened nodes are conditioned on and the queried nodes are not: there
# the C++ procedure is d-separation in the graph without the edges into the
# intervened nodes. (Outside that domain the C++ code keeps the bidirected
# edges of an intervened node in one direction and can differ.) The tests
# compare against a literal port of the C++ code.

"""
    DsepGraph(n)

The graph that the do-calculus rules test d-separation in: `n` variables
`1:n`, their intervention nodes `n + 1:2n` with edges `I(x) -> x`, directed
edges and bidirected edges (latent confounders). Add edges with
[`add_directed!`](@ref) and [`add_bidirected!`](@ref).
"""
struct DsepGraph
    nvars::Int
    children::Vector{UInt64}     # children[i] has bit j-1 for every edge i -> j
    parents::Vector{UInt64}      # parents[j] has bit i-1 for every edge i -> j
    confounded::Vector{UInt64}   # confounded[i] has bit j-1 for every edge i <-> j
    function DsepGraph(n::Integer)
        n <= 32 || throw(ArgumentError("at most 32 variables are supported, got $n"))
        children = zeros(UInt64, 2n)
        parents = zeros(UInt64, 2n)
        for x in 1:n
            children[n + x] |= one(UInt64) << (x - 1)     # I(x) -> x
            parents[x] |= one(UInt64) << (n + x - 1)
        end
        return new(n, children, parents, zeros(UInt64, 2n))
    end
end

"""
    add_directed!(g, from, to)

Add the edge `from -> to` between variables of `g`.
"""
function add_directed!(g::DsepGraph, from::Integer, to::Integer)
    g.children[from] |= one(UInt64) << (to - 1)
    g.parents[to] |= one(UInt64) << (from - 1)
    return g
end

"""
    add_bidirected!(g, x, y)

Add the bidirected edge `x <-> y` (a latent common cause) to `g`.
"""
function add_bidirected!(g::DsepGraph, x::Integer, y::Integer)
    g.confounded[x] |= one(UInt64) << (y - 1)
    g.confounded[y] |= one(UInt64) << (x - 1)
    return g
end

nnodes(g::DsepGraph) = 2 * g.nvars

@inline bit(i::Integer) = one(UInt64) << (i - 1)
@inline hasbit(row::UInt64, i::Integer) = (row >> (i - 1)) & one(UInt64) == one(UInt64)

"""
    dconnected(g, x, given, intervened, stop = 0) -> UInt64

The nodes that are d-connected to some node of `x` given `given` in `g` after
intervening on `intervened` (all as bitmasks), by a reachability search over
active paths (as in the Bayes-ball algorithm). Intervening on a node removes
its incoming directed edges and its bidirected edges. The search stops early,
with a partial result, once it reaches a node of `stop`.
"""
function dconnected(g::DsepGraph, x::UInt64, c::UInt64, j::UInt64, stop::UInt64 = zero(UInt64))
    # ancestors of the conditioning set, after intervening
    anc = c
    todo = c
    while todo != 0
        v = trailing_zeros(todo) + 1
        todo &= todo - one(UInt64)
        new = (hasbit(j, v) ? zero(UInt64) : g.parents[v]) & ~anc
        anc |= new
        todo |= new
    end
    # Nodes reached along an active path, by whether the last edge has an
    # arrowhead at the node (`head`) or a tail (`tail`).
    head = zero(UInt64)
    tail = zero(UInt64)
    qhead = zero(UInt64)
    qtail = zero(UInt64)
    xs = x
    while xs != 0
        v = trailing_zeros(xs) + 1
        xs &= xs - one(UInt64)
        intervened = hasbit(j, v)
        qhead |= (g.children[v] & ~j) | (intervened ? zero(UInt64) : g.confounded[v] & ~j)
        qtail |= intervened ? zero(UInt64) : g.parents[v]
    end
    while (qhead | qtail) != 0
        newh = zero(UInt64)
        newt = zero(UInt64)
        while qtail != 0     # arrived through a tail: never a collider
            v = trailing_zeros(qtail) + 1
            qtail &= qtail - one(UInt64)
            tail |= bit(v)
            hasbit(c, v) && continue
            intervened = hasbit(j, v)
            newh |= (g.children[v] & ~j) | (intervened ? zero(UInt64) : g.confounded[v] & ~j)
            newt |= intervened ? zero(UInt64) : g.parents[v]
        end
        while qhead != 0     # arrived through an arrowhead
            v = trailing_zeros(qhead) + 1
            qhead &= qhead - one(UInt64)
            head |= bit(v)
            hasbit(c, v) || (newh |= g.children[v] & ~j)        # non-collider
            if hasbit(anc, v) && !hasbit(j, v)                  # collider
                newh |= g.confounded[v] & ~j
                newt |= g.parents[v]
            end
        end
        reached = head | tail | newh | newt
        reached & stop != 0 && return reached
        qhead = newh & ~head
        qtail = newt & ~tail
    end
    return head | tail
end

"""
    dseparated(g, s::Separation) -> Bool

Whether every node of `s.x` is d-separated from every node of `s.y` given
`s.given` in `g` after intervening on `s.intervened`. See [`dconnected`](@ref).
"""
dseparated(g::DsepGraph, s::Separation) =
    dconnected(g, s.x.bits, s.given.bits, s.intervened.bits, s.y.bits) & s.y.bits == 0
