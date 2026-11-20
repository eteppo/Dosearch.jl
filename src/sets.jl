# Variable sets as bitmasks.
#
# The search manipulates sets of variables in its innermost loops, so sets are
# stored as a single machine word: variable `i` is bit `i - 1`. A model with
# `n` variables also uses bits `n + 1:2n` for the intervention nodes `I(x)` of
# the d-separation graph, so at most 32 variables are supported.

"""
    VarSet

An immutable set of variable indices in `1:64`, stored as a bitmask.

`VarSet` supports the usual set operations (`∪`, `∩`, `setdiff`, `⊆`, `in`,
`isempty`, `length`) and iterates over its elements in increasing order.
Construct one from a collection of indices, as in `VarSet([1, 3])` or
`VarSet(1:3)`, or from a raw bitmask with `VarSet(mask::UInt64)`.
"""
struct VarSet
    bits::UInt64
    VarSet(bits::UInt64) = new(bits)
end

VarSet() = VarSet(zero(UInt64))
function VarSet(indices::Union{AbstractVector{<:Integer},Tuple{Vararg{Integer}},AbstractSet{<:Integer}})
    bits = zero(UInt64)
    for i in indices
        1 <= i <= 64 || throw(ArgumentError("variable index $i is outside 1:64"))
        bits |= one(UInt64) << (i - 1)
    end
    return VarSet(bits)
end

"""
    fullset(n) -> VarSet

The set `{1, …, n}`.
"""
fullset(n::Integer) = n >= 64 ? VarSet(typemax(UInt64)) : VarSet((one(UInt64) << n) - one(UInt64))

Base.union(a::VarSet, b::VarSet) = VarSet(a.bits | b.bits)
Base.union(a::VarSet, b::VarSet, c::VarSet...) = union(union(a, b), c...)
Base.intersect(a::VarSet, b::VarSet) = VarSet(a.bits & b.bits)
Base.setdiff(a::VarSet, b::VarSet) = VarSet(a.bits & ~b.bits)
Base.issubset(a::VarSet, b::VarSet) = a.bits & ~b.bits == 0
Base.isempty(a::VarSet) = a.bits == 0
Base.length(a::VarSet) = count_ones(a.bits)
Base.in(i::Integer, a::VarSet) = 1 <= i <= 64 && (a.bits >> (i - 1)) & one(UInt64) == one(UInt64)
Base.hash(a::VarSet, h::UInt) = hash(a.bits, hash(:VarSet, h))
Base.isless(a::VarSet, b::VarSet) = isless(a.bits, b.bits)

# Shifts move every element by `k` positions. They map variables to their
# intervention nodes and, in missing data problems, between the true variable,
# response indicator and proxy of the same variable, which are numbered
# consecutively.
Base.:<<(a::VarSet, k::Integer) = VarSet(a.bits << k)
Base.:>>(a::VarSet, k::Integer) = VarSet(a.bits >> k)

Base.eltype(::Type{VarSet}) = Int
function Base.iterate(a::VarSet, bits::UInt64 = a.bits)
    bits == 0 && return nothing
    return trailing_zeros(bits) + 1, bits & (bits - one(UInt64))
end

function Base.show(io::IO, a::VarSet)
    print(io, "VarSet(")
    join(io, a, ", ")
    print(io, ")")
end

"""
    deposit(ranks::VarSet, mask::VarSet) -> VarSet

The subset of `mask` that contains the `i`th smallest element of `mask` for
every `i` in `ranks` (a parallel bit deposit).
"""
function deposit(ranks::VarSet, mask::VarSet)
    src, m, out = ranks.bits, mask.bits, zero(UInt64)
    while src != 0 && m != 0
        low = m & (~m + one(UInt64))      # lowest set bit of m
        src & one(UInt64) != 0 && (out |= low)
        src >>= 1
        m &= m - one(UInt64)
    end
    return VarSet(out)
end

"""
    subsets_by_size(n) -> Vector{VarSet}

All non-empty subsets of `{1, …, n}` in the order the search tries them: by
increasing cardinality, and lexicographically within one cardinality. This is
the banker's sequence of `get_subsets` in the R package's `set.cpp`.
"""
function subsets_by_size(n::Integer)
    sets = VarSet[]
    sizehint!(sets, n >= 62 ? 0 : (1 << n) - 1)
    idx = Vector{Int}(undef, n)
    for k in 1:n
        for j in 1:k
            idx[j] = j
        end
        while true
            push!(sets, VarSet(view(idx, 1:k)))
            i = k
            while i >= 1 && idx[i] == n - k + i
                i -= 1
            end
            i == 0 && break
            idx[i] += 1
            for j in (i + 1):k
                idx[j] = idx[j - 1] + 1
            end
        end
    end
    return sets
end
