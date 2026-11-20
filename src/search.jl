# The search core: derive new terms from the known ones with the inference
# rules until the query is found. Port of `search::find` and
# `search::enumerate_distribution`, which the R package shares between the
# search for DAGs (`dosearch`) and for LDAGs (`csisearch`). Here the core is
# generic over a `SearchModel`, and each model provides the rules and the
# operations that differ (see "Model interface" below).

"""
    SearchModel

Abstract type of the models that the search core works with: a problem with
its rules, its separation criterion and its output. See `DagModel` and
`LdagModel`.
"""
abstract type SearchModel end

"""
    Node{R,S}

A term found by the search and how it was derived: `rule::R` applied to the
term of node `parent1`, using the term of node `parent2` if the rule needs a
second term (`parent2 == 0` otherwise). Input distributions have no parents.
A node is `primitive` if it was derived from input distributions by standard
probability manipulations only; the formula then shows the term itself.
`sep::S` is the separation statement that the rule needed (empty if none).
"""
struct Node{R,S}
    dist::Dist
    rule::R
    parent1::Int
    parent2::Int
    primitive::Bool
    sep::S
end

"""
    SearchOptions{R}

Options of the search core. See [`dosearch`](@ref) for their meaning.
"""
Base.@kwdef struct SearchOptions{R}
    rules::Vector{R}
    time_limit::Float64 = Inf
    term_limit::Int = typemax(Int)
    benchmark_rules::Bool = false
    improve::Bool = true
    verbose::Bool = false
    io::IO = stdout       # where the verbose log goes
end

mutable struct SearchState{M<:SearchModel,R,S,F}
    const model::M
    const options::SearchOptions{R}
    const nodes::Vector{Node{R,S}}       # every term found, by index
    const index::Dict{Dist,Int}          # node index of each term
    const frontier::F
    const subsets::Vector{Vector{VarSet}}   # banker's sequences of 1:k, by k + 1
    const candidates::Vector{Int}        # buffer for `candidates!`
    const rule_times::Vector{Float64}    # seconds spent in each rule
    count::Int                           # number of stored terms
    remaining::Int                       # terms not yet expanded
    target::Union{Nothing,Node{R,S}}     # the query, once derived
    target_count::Int                    # index of the query node
    limit_reached::Bool                  # stopped at the time or term limit
end

# Model interface. Each model defines:
#
# - `problem(m)`: the problem, with `data` and `query`;
# - `input_rule(m)`, `no_separation(m)`: the rule value of input
#   distributions and the empty separation statement;
# - `left_variables(m)`: variables that can appear on the left of a term;
# - `rule_applies(m, rule, p)`, `subset_mask(m, rule, p)`,
#   `apply_rule(m, rule, p, z)`: the rules;
# - `subset_count(m, rule, k)`: how many subsets of a mask with `k` elements
#   the rule tries, in banker's order;
# - `term_data(m, p)` and `narrow_mask(m, rule, p, mask, data)`: optional
#   per-term precomputation that removes subsets that are sure to fail;
# - `separation_holds(m, sep)`: the separation criterion;
# - `candidates!(m, st, inf, from)` and `candidate_term(m, to, required)`:
#   the required terms of enumerating rules;
# - `is_primitive(m, rule, primitive1, primitive2)`;
# - `target_index(m, index)`: the index shown in the verbose log;
# - `input_score(m, p)`, `heuristic_score(m, p)`: scores for the heuristic;
# - `term_string(problem, p)` and `rule_label(rule)`: output.
#
# The output of an identified query is defined on the problem types instead:
# `rewrite` (formula.jl) and `justification` (derivation.jl).

term_data(::SearchModel, p::Dist) = nothing
narrow_mask(::SearchModel, rule, p::Dist, mask::VarSet, data) = mask, false

function SearchState(model::SearchModel, strategy::SearchStrategy, options::SearchOptions{R}) where {R}
    S = typeof(no_separation(model))
    return SearchState(model, options, Node{R,S}[], Dict{Dist,Int}(), frontier(strategy),
                       Vector{VarSet}[], Int[], zeros(length(options.rules)), 0, 0, nothing, 0, false)
end

problem(st::SearchState) = problem(st.model)
term_string(st::SearchState, p::Dist) = term_string(problem(st.model), p)

# Add an input distribution (`add_known`).
function add_input!(st::SearchState, p::Dist)
    st.count += 1
    node = Node(p, input_rule(st.model), 0, 0, true, no_separation(st.model))
    push!(st.nodes, node)
    st.index[p] = st.count
    push_frontier!(st.frontier, st, st.count, input_score(st.model, p))
    if p == problem(st).query && st.target === nothing
        st.target = node
        st.target_count = st.count
    end
    st.options.verbose && println(st.options.io, "Adding known distribution: ", term_string(st, p))
    return st
end

push_frontier!(::FoundOrderFrontier, st, i, score) = nothing
push_frontier!(f::HeapFrontier, st, i, score) = heap_push!(f, i, score)

next_node!(::FoundOrderFrontier, st, i) = i
next_node!(f::HeapFrontier, st, i) = heap_pop!(f)

new_score(::FoundOrderFrontier, st, p) = 0
new_score(::HeapFrontier, st, p) = heuristic_score(st.model, p)

"""
    trivially_nonidentifiable(st) -> Bool

The query is not identifiable if one of its left-hand side variables never
appears on the left-hand side of an input (`check_trivial`).
"""
trivially_nonidentifiable(st::SearchState) = !(problem(st).query.a ⊆ left_variables(st.model))

# Record a derived term (`derive_distribution`). Returns true if it is the query.
function derive!(st::SearchState, inf::Inference, to::Dist, from::Node, from_index::Int,
                 required::Int, rule)
    has_required = !isempty(inf.required.a)
    required_primitive = has_required ? st.nodes[required].primitive : true
    primitive = is_primitive(st.model, rule, from.primitive, required_primitive)
    node = Node(to, rule, from_index, has_required ? required : 0, primitive, inf.sep)
    if st.options.verbose
        io = st.options.io
        print(io, "Derived: ", term_string(st, to), " from ", term_string(st, from.dist))
        has_required && print(io, " and ", term_string(st, inf.required))
        println(io, " using rule: ", Int(rule))
    end
    if to == problem(st).query
        st.target = node
        st.target_count = target_index(st.model, st.count + 1)
        if st.options.verbose
            println(st.options.io, "Target found")
            println(st.options.io, "Index = ", st.target_count)
        end
        return true
    end
    st.count += 1
    st.remaining += 1
    push!(st.nodes, node)
    st.index[to] = st.count
    push_frontier!(st.frontier, st, st.count, new_score(st.frontier, st, to))
    return false
end

# Try `rule` with subset `z` on node `i` (`enumerate_distribution`). Returns
# true if the query was derived. With `prechecked = true` the caller has
# already established the separation that the rule needs.
function expand!(st::SearchState, i::Int, rule, z::VarSet, prechecked::Bool = false)
    m = st.model
    from = st.nodes[i]
    inf = apply_rule(m, rule, from.dist, z)
    inf === nothing && return false
    if inf.enumerate
        for req in candidates!(m, st, inf, from.dist)
            to = candidate_term(m, inf.to, st.nodes[req].dist)
            if !haskey(st.index, to)
                derive!(st, inf, to, from, i, req, rule) && return true
            end
        end
        return false
    end
    haskey(st.index, inf.to) && return false
    if !prechecked && !isempty(inf.sep.x) && !separation_holds(m, inf.sep)
        return false
    end
    if !isempty(inf.required.a)
        req = get(st.index, inf.required, 0)
        req > 0 || return false
        return derive!(st, inf, inf.to, from, i, req, rule)
    end
    return derive!(st, inf, inf.to, from, i, 0, rule)
end

# Banker's sequence of the non-empty subsets of 1:k, built on first use.
function subsets_of_size(st::SearchState, k::Int)
    while length(st.subsets) <= k
        push!(st.subsets, subsets_by_size(length(st.subsets)))
    end
    return st.subsets[k + 1]
end

# The main loop (`search::find`).
function find!(st::SearchState)
    m = st.model
    rules = st.options.rules
    timed = isfinite(st.options.time_limit)
    start = time()
    i = 1
    while st.remaining > 0
        if (timed && time() - start >= st.options.time_limit) || st.count >= st.options.term_limit
            st.limit_reached = true
            return false
        end
        idx = next_node!(st.frontier, st, i)
        st.remaining -= 1
        p = st.nodes[idx].dist
        data = term_data(m, p)
        for (r, rule) in enumerate(rules)
            t0 = st.options.benchmark_rules ? time() : 0.0
            found = false
            if rule_applies(m, rule, p)
                # The C++ code tries the subsets of the variables in banker's
                # order; only subsets of the mask can pass the rule. Removing
                # elements from the mask keeps the order of the other subsets.
                mask = subset_mask(m, rule, p)
                mask, prechecked = narrow_mask(m, rule, p, mask, data)
                k = length(mask)
                for ranks in view(subsets_of_size(st, k), 1:subset_count(m, rule, k))
                    if expand!(st, idx, rule, deposit(ranks, mask), prechecked)
                        found = true
                        break
                    end
                end
            end
            st.options.benchmark_rules && (st.rule_times[r] += time() - t0)
            found && return true
        end
        i += 1
    end
    return false
end

"""
    run_search!(st) -> Bool

Run the search from the input distributions of the problem. Returns whether
the query was found.
"""
function run_search!(st::SearchState)
    for p in problem(st).data
        add_input!(st, p)
    end
    st.options.verbose && println(st.options.io, "Initializing search")
    if st.options.improve && trivially_nonidentifiable(st)
        return false
    end
    st.remaining = length(st.nodes)
    st.target === nothing && find!(st)
    return st.target !== nothing
end

"""
    derivation_separations(st) -> Vector

The separation statements that the derivation of the query relies on. The
query stays derivable in every graph where all of them hold, because the
other conditions of the rules do not depend on the graph.
"""
function derivation_separations(st::SearchState{M,R,S}) where {M,R,S}
    seps = S[]
    st.target === nothing && return seps
    seen = falses(length(st.nodes))
    stack = [st.target]
    while !isempty(stack)
        node = pop!(stack)
        isempty(node.sep.x) || node.sep in seps || push!(seps, node.sep)
        for i in (node.parent1, node.parent2)
            if i > 0 && !seen[i]
                seen[i] = true
                push!(stack, st.nodes[i])
            end
        end
    end
    return seps
end
