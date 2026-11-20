# Order in which the search expands the terms it has derived.
#
# The R package has two search classes that differ only in this order: plain
# `dosearch` expands terms in the order they were found, `dosearch_heuristic`
# keeps them in a `std::priority_queue` ordered by a score. Here the order is a
# strategy value that the search core dispatches on.

"""
    SearchStrategy

Abstract type of the strategies that choose which derived term the search
expands next. See [`FoundOrder`](@ref) and [`ScoreHeuristic`](@ref).
"""
abstract type SearchStrategy end

"""
    FoundOrder()

Expand terms in the order they were derived (`heuristic = FALSE` in R).
"""
struct FoundOrder <: SearchStrategy end

"""
    ScoreHeuristic()

Expand the term that looks closest to the query first, by the score of Tikka
et al. (2021) (`heuristic = TRUE` in R).
"""
struct ScoreHeuristic <: SearchStrategy end

# Mutable frontier state for a strategy.
frontier(::FoundOrder) = FoundOrderFrontier()
frontier(::ScoreHeuristic) = HeapFrontier()

struct FoundOrderFrontier end

"""
    HeapFrontier

A binary max-heap of node indices ordered by score that reproduces
`std::priority_queue` of GCC's libstdc++ step by step, so that ties between
equal scores are broken exactly as in the R package (compiled with GCC).
"""
struct HeapFrontier
    heap::Vector{Int}       # node indices, 0-based heap layout stored 1-based
    scores::Vector{Int}     # score of each node index
end
HeapFrontier() = HeapFrontier(Int[], Int[])

# std::__push_heap: move `value` up from `hole` (0-based) towards `top`.
function heap_sift_up!(h::HeapFrontier, hole::Int, top::Int, value::Int)
    heap, scores = h.heap, h.scores
    parent = (hole - 1) ÷ 2
    while hole > top && scores[heap[parent + 1]] < scores[value]
        heap[hole + 1] = heap[parent + 1]
        hole = parent
        parent = (hole - 1) ÷ 2
    end
    heap[hole + 1] = value
    return h
end

function heap_push!(h::HeapFrontier, node::Int, score::Int)
    length(h.scores) < node && resize!(h.scores, max(node, 2 * length(h.scores)))
    h.scores[node] = score
    push!(h.heap, node)
    return heap_sift_up!(h, length(h.heap) - 1, 0, node)
end

# std::pop_heap followed by pop_back: libstdc++ moves the hole from the root
# down to a leaf, always to the larger child (the right one on ties), and then
# sifts the former last element up from there.
function heap_pop!(h::HeapFrontier)
    heap, scores = h.heap, h.scores
    top = heap[1]
    len = length(heap) - 1          # heap size after removing the top
    if len > 0
        value = heap[end]
        hole = 0
        child = 0
        while child < (len - 1) ÷ 2
            child = 2 * (child + 1)
            if scores[heap[child + 1]] < scores[heap[child]]
                child -= 1
            end
            heap[hole + 1] = heap[child + 1]
            hole = child
        end
        if len & 1 == 0 && child == (len - 2) ÷ 2
            child = 2 * (child + 1)
            heap[hole + 1] = heap[child]
            hole = child - 1
        end
        heap_sift_up!(h, hole, 0, value)
    end
    pop!(heap)
    return top
end

"""
    heuristic_score(p, target) -> Int

Score of term `p` in the search for `target`: shared variables in the same
role raise the score, unshared ones lower it (`dosearch_heuristic::compute_score`).
"""
function heuristic_score(p::Dist, target::Dist)
    common_y = p.a ∩ target.a
    common_x = p.c ∩ target.c
    p_w = setdiff(p.b, p.c)
    target_w = setdiff(target.b, target.c)
    common_z = p_w ∩ target_w
    return 10 * length(common_y) - 2 * length(setdiff(target.a, common_y)) +
           5 * length(common_x) - 2 * length(setdiff(p.c, common_x)) -
           2 * length(setdiff(target.c, common_x)) + 3 * length(common_z) -
           length(setdiff(p_w, common_z)) - length(setdiff(target_w, common_z))
end
