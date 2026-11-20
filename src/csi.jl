# Context-specific independence (CSI) in an LDAG. Port of `ldag` and
# `ldag_cache` in the R package's `ldag.cpp`.
#
# A CSI statement "x ⊥ y | z in the context zero = 0, one = 1" holds if it is
# encoded directly by a label, if x and y are d-separated given z in the graph
# without the edges that vanish in the context, or, recursively, if it holds
# in every context of further context variables (reasoning by cases). The
# recursion follows the C++ code exactly, including its order and its cache of
# statements that were found to hold, because the cache can change the result:
# a statement found to hold in a less restricted part of the recursion is also
# accepted in a more restricted one.
#
# Intentional difference: the contexts of a statement are looked up by the
# assignments of the variables that appear in labels (`label_vars`). The C++
# code looks them up by all assigned context variables, which also include
# context variables given with `con_vars` that appear in no label; it then
# finds no contexts, and its reasoning by cases accepts any statement, such as
# the independence of a variable and its child. The two agree whenever the
# assigned context variables all appear in labels.

struct ContextSetting
    zero::UInt64
    one::UInt64
    equiv::UInt64      # bit of the equivalence class of the context
end

"""
    CsiGraph(problem::LdagProblem; cache = true)

The graph of an LDAG problem with its contexts and local CSIs, for testing
CSI statements with [`csi_holds!`](@ref).
"""
mutable struct CsiGraph
    const n::Int
    const children::Vector{UInt64}
    const parents::Vector{UInt64}
    const con_vars::UInt64
    const label_vars::UInt64                 # context variables that appear in labels
    const intv_vars::UInt64
    const context_sets::Vector{UInt64}
    const contexts::Dict{Tuple{UInt64,UInt64},Vector{Tuple{Int,Int}}}   # (zero, one) => removed edges
    const settings::Dict{UInt64,Vector{ContextSetting}}                  # zero | one => contexts
    const local_csi::Vector{NTuple{5,UInt64}}
    const use_cache::Bool
    const cache::Set{NTuple{5,UInt64}}       # (x, y, z, zero, one | intv) that hold
    const ch::Vector{UInt64}                 # edges of the current context
    const pa::Vector{UInt64}
end

function CsiGraph(p::LdagProblem; cache::Bool = true)
    n = problem_nvars(p)
    children = zeros(UInt64, n)
    parents = zeros(UInt64, n)
    for (from, to) in p.edges
        children[from] |= bit(to)
        parents[to] |= bit(from)
    end
    contexts = Dict{Tuple{UInt64,UInt64},Vector{Tuple{Int,Int}}}()
    settings = Dict{UInt64,Vector{ContextSetting}}()
    context_sets = UInt64[]
    for cs in p.context_sets
        push!(context_sets, cs.vars.bits)
        for c in cs.contexts
            contexts[(c.zero.bits, c.one.bits)] = c.removed
            push!(get!(settings, c.zero.bits | c.one.bits, ContextSetting[]),
                  ContextSetting(c.zero.bits, c.one.bits, bit(c.equivalence)))
        end
    end
    local_csi = [(l.x.bits, l.y.bits, l.z.bits, l.zero.bits, l.one.bits) for l in p.local_csi]
    label_vars = reduce(|, context_sets; init = zero(UInt64)) & ~p.intervention_vars.bits
    return CsiGraph(n, children, parents, p.context_vars.bits, label_vars, p.intervention_vars.bits,
                    context_sets, contexts, settings, local_csi, cache, Set{NTuple{5,UInt64}}(),
                    zeros(UInt64, n), zeros(UInt64, n))
end

const NO_EDGES = Tuple{Int,Int}[]

# Ancestors of `set` (inclusive) along the edges `pa`.
function ancestors(pa::Vector{UInt64}, set::UInt64)
    anc = set
    todo = set
    while todo != 0
        v = trailing_zeros(todo) + 1
        todo &= todo - one(UInt64)
        new = pa[v] & ~anc
        anc |= new
        todo |= new
    end
    return anc
end

# d-separation of x and y given z in the current context (`ldag::d_sep`):
# Bayes-ball restricted to the ancestors of x, y and z. A node is visited
# `up` when reached from a child and `down` when reached from a parent.
function dseparated_current(g::CsiGraph, x::UInt64, y::UInt64, z::UInt64)
    ch, pa = g.ch, g.pa
    a = ancestors(pa, z)
    xyz = ancestors(pa, x | y | z)
    up = zero(UInt64)
    down = zero(UInt64)
    qup = x
    qdown = zero(UInt64)
    while (qup | qdown) != 0
        newup = zero(UInt64)
        newdown = zero(UInt64)
        while qup != 0
            v = trailing_zeros(qup) + 1
            qup &= qup - one(UInt64)
            hasbit(up, v) && continue
            up |= bit(v)
            hasbit(y, v) && return false
            if !hasbit(z, v)
                newup |= pa[v] & xyz
                newdown |= ch[v] & xyz
            end
        end
        while qdown != 0
            v = trailing_zeros(qdown) + 1
            qdown &= qdown - one(UInt64)
            hasbit(down, v) && continue
            down |= bit(v)
            hasbit(y, v) && return false
            hasbit(z, v) || (newdown |= ch[v] & xyz)
            hasbit(a, v) && (newup |= pa[v] & xyz)
        end
        qup = newup & ~up
        qdown = newdown & ~down
    end
    return true
end

# d-separation in the graph without the edges of two contexts (`csi_sep`).
function dseparated_in(g::CsiGraph, x::UInt64, y::UInt64, z::UInt64, removed1, removed2)
    copyto!(g.ch, g.children)
    copyto!(g.pa, g.parents)
    for removed in (removed1, removed2), (from, to) in removed
        g.ch[from] &= ~bit(to)
        g.pa[to] &= ~bit(from)
    end
    return dseparated_current(g, x, y, z)
end

# Is the statement encoded directly by a label (`in_label`)?
function in_label(g::CsiGraph, x, y, z, c0, c1)
    for (lx, ly, lz, lzero, lone) in g.local_csi
        if x == lx
            (y != ly || z != lz || c0 != lzero || c1 != lone) && continue
            return true
        elseif y == lx
            (x != ly || z != lz || c0 != lzero || c1 != lone) && continue
            return true
        end
    end
    return false
end

cached(g::CsiGraph, x, y, z, c0, c1) = g.use_cache && (x, y, z, c0, c1) in g.cache

function remember!(g::CsiGraph, x, y, z, c0, c1)
    if g.use_cache
        push!(g.cache, (x, y, z, c0, c1))
        push!(g.cache, (y, x, z, c0, c1))
    end
    return true
end

"""
    csi_holds!(g, x, y, z, c0, c1, intv, old_con) -> Bool

Whether `x` and `y` are independent given `z` in the context where `c0` are
0 and `c1` and the intervention nodes `intv` are 1, without reasoning by
cases over the context variables in `old_con` (`ldag::csi_criterion`, or
`ldag_cache::csi_criterion` with the cache). All sets are bitmasks.
"""
function csi_holds!(g::CsiGraph, x::UInt64, y::UInt64, z::UInt64, c0::UInt64, c1::UInt64,
                    intv::UInt64, old_con::UInt64)
    c1i = c1 | intv
    cached(g, x, y, z, c0, c1i) && return true
    in_label(g, x, y, z, c0, c1i) && return remember!(g, x, y, z, c0, c1i)
    l0 = c0 & g.label_vars
    l1 = c1 & g.label_vars
    con = get(g.contexts, (l0, l1), NO_EDGES)
    ivar = get(g.contexts, (zero(UInt64), intv), NO_EDGES)
    dseparated_in(g, x, y, z, con, ivar) && return remember!(g, x, y, z, c0, c1i)
    # Reason by cases over context variables in z.
    zc = z & ~c0 & ~c1 & g.con_vars
    if zc != 0
        for w in g.context_sets
            w & old_con != 0 && continue
            w & zc != w && continue
            valid = true
            new_con = old_con | w
            equiv_class = zero(UInt64)
            for s in get(g.settings, (l0 | l1) | w, ContextSetting[])
                (s.zero & l0 != l0 || s.one & l1 != l1) && continue
                s.equiv & equiv_class == s.equiv && continue   # class already checked
                equiv_class |= s.equiv
                s0, s1 = s.zero | c0, s.one | c1
                cached(g, x, y, z, s0, s1 | intv) && continue
                if !csi_holds!(g, x, y, z, s0, s1, intv, new_con)
                    valid = false
                    break
                end
                remember!(g, x, y, z, s0, s1 | intv)
            end
            valid && equiv_class != 0 && return remember!(g, x, y, z, c0, c1i)
        end
    end
    # Reason by cases over context variables outside the statement.
    for w in g.context_sets
        w & old_con != 0 && continue                  # cannot reuse old contexts
        w & g.intv_vars != 0 && continue              # w cannot contain interventions
        (w & x != 0 || w & y != 0 || w & z != 0) && continue
        valid = false
        new_con = old_con | w
        equiv_class = zero(UInt64)
        if csi_holds!(g, x, w, z | y, c0, c1, intv, new_con) ||
           csi_holds!(g, x, w, z, c0, c1, intv, new_con) ||
           csi_holds!(g, y, w, z | x, c0, c1, intv, new_con) ||
           csi_holds!(g, y, w, z, c0, c1, intv, new_con)
            valid = true
            for s in get(g.settings, (l0 | l1) | w, ContextSetting[])
                (s.zero & l0 != l0 || s.one & l1 != l1) && continue
                s.equiv & equiv_class == s.equiv && continue
                equiv_class |= s.equiv
                s0, s1 = s.zero | c0, s.one | c1
                cached(g, x, y, z | w, s0, s1 | intv) && continue
                if !csi_holds!(g, x, y, z | w, s0, s1, intv, new_con)
                    valid = false
                    break
                end
                remember!(g, x, y, z | w, s0, s1 | intv)
            end
        end
        valid && equiv_class != 0 && return remember!(g, x, y, z, c0, c1i)
    end
    return false
end
