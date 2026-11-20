# Search for maximal identifiable graphs: the largest sets of candidates (such
# as latent confounders) that can be added to a model while the query stays
# identifiable.
#
# Adding a candidate only adds paths to the graph, so it only removes
# d-separations, and every derivation of the search in the larger model is
# also valid in the smaller one. The identifiable sets of candidates are
# therefore closed under taking subsets, and the task is to find the maximal
# elements of this family with as few searches as possible. The oracle below
# answers from earlier searches whenever monotonicity decides the answer, and
# reuses earlier derivations: a set is identifiable without a search if every
# d-separation that some earlier derivation of the query used still holds.

# Identifiability of candidate sets, with answers inferred by monotonicity.
# Sets that need a search can be checked in batches on several threads; the
# results are recorded in a fixed order, so the answers and everything that
# depends on them are the same for any number of threads (unless a time limit
# stops searches).
#
# A search that stops at its time or term limit counts as non-identifiable,
# and the set is also recorded as `undecided`. The answers stay consistent
# with monotonicity, since a set is only searched when no earlier answer
# decides it: the oracle answers for a monotone family that lies inside the
# identifiable sets.
mutable struct MonotoneOracle{P<:DosearchProblem,C,K,T}
    const problem::P
    const candidates::Vector{C}
    const search_options::K
    const identifiable::Vector{BitSet}       # maximal known identifiable sets
    const nonidentifiable::Vector{BitSet}    # minimal known non-identifiable sets
    const certificates::Vector{T}            # of earlier derivations
    const undecided::Vector{BitSet}          # searched sets that hit a limit
    const use_certificates::Bool
    const parallel::Bool                     # search batches on several threads
    searches::Int                            # dosearch calls
    certified::Int                           # answers from certificates
    identifiable_time::Float64               # seconds in calls that identified
    nonidentifiable_time::Float64            # seconds in calls that did not
end

MonotoneOracle(problem, candidates, options; use_certificates::Bool = true,
               parallel::Bool = Threads.nthreads() > 1) =
    MonotoneOracle(problem, candidates, options, BitSet[], BitSet[], certificate_type(problem)[], BitSet[],
                   use_certificates, parallel, 0, 0, 0.0, 0.0)

# The known answer for `s`, or `nothing`.
function known(o::MonotoneOracle, s::BitSet)
    any(t -> s ⊆ t, o.identifiable) && return true
    any(t -> t ⊆ s, o.nonidentifiable) && return false
    return nothing
end

function record!(o::MonotoneOracle, s::BitSet, id::Bool)
    if id   # keep antichains: drop the sets that the new one implies
        filter!(t -> !(t ⊆ s), o.identifiable)
        push!(o.identifiable, copy(s))
    else
        filter!(t -> !(s ⊆ t), o.nonidentifiable)
        push!(o.nonidentifiable, copy(s))
    end
    return id
end

# Certificates: the separation statements that a derivation of the query
# relies on. The search uses no other property of the graph, so the
# derivation stays valid in every model where they hold. In a DAG, adding
# candidates keeps the variable numbering. In an LDAG, the latent nodes of
# the candidates depend on the set, so the statements keep the names of
# their variables.
certificate_type(::DagProblem) = Vector{Separation}
certificate(st::SearchState{DagModel}) = derivation_separations(st)

struct LdagCertificate
    labels::Vector{String}
    separations::Vector{CsiSeparation}
end
certificate_type(::LdagProblem) = LdagCertificate
certificate(st::SearchState{LdagModel}) = LdagCertificate(problem(st).labels, derivation_separations(st))

# Does some earlier derivation of the query stay valid in `problem`?
function certified(o::MonotoneOracle, problem::DagProblem)
    isempty(o.certificates) && return false
    g = causal_graph(problem)
    return any(cert -> all(sep -> dseparated(g, sep), cert), Iterators.reverse(o.certificates))
end

function certified(o::MonotoneOracle, problem::LdagProblem)
    isempty(o.certificates) && return false
    m = LdagModel(problem)
    index = Dict(l => i for (i, l) in enumerate(problem.labels))
    return any(Iterators.reverse(o.certificates)) do cert
        all(cert.separations) do sep
            s = translate(sep, cert.labels, index)
            s !== nothing && separation_holds(m, s)
        end
    end
end

# A variable set of a problem with variables `labels` in the numbering
# `index` of another problem, or `nothing` if a variable is missing there.
function translate(s::VarSet, labels::Vector{String}, index::Dict{String,Int})
    bits = zero(UInt64)
    for i in s
        j = get(index, labels[i], 0)
        j == 0 && return nothing
        bits |= one(UInt64) << (j - 1)
    end
    return VarSet(bits)
end

function translate(s::CsiSeparation, labels::Vector{String}, index::Dict{String,Int})
    parts = map(v -> translate(v, labels, index), (s.x, s.y, s.z, s.zero, s.one))
    any(isnothing, parts) && return nothing
    return CsiSeparation(parts...)
end

# Search `problem`: (identifiable, seconds, certificate or nothing, whether the
# search stopped at a limit).
function search(o::MonotoneOracle, problem::DosearchProblem)
    st, elapsed = solve(problem; o.search_options...)
    id = st.target !== nothing
    cert = id && o.use_certificates ? certificate(st) : nothing
    return (id, elapsed, cert, st.limit_reached)
end

"""
    check!(o, sets)

Decide the identifiability of every set in `sets`: from known answers and
certificates where possible, and otherwise by searches, in parallel if
`o.parallel`.
"""
function check!(o::MonotoneOracle{P,C,K,T}, sets::AbstractVector{BitSet}) where {P,C,K,T}
    pending = BitSet[]
    problems = P[]
    for s in sets
        (known(o, s) === nothing && !(s in pending)) || continue
        problem = with_candidates(o.problem, o.candidates[collect(s)])
        if certified(o, problem)
            o.certified += 1
            record!(o, s, true)
        else
            push!(pending, s)
            push!(problems, problem)
        end
    end
    results = Vector{Tuple{Bool,Float64,Union{Nothing,T},Bool}}(undef, length(pending))
    if o.parallel && length(pending) > 1
        Threads.@threads for i in eachindex(problems)
            results[i] = search(o, problems[i])
        end
    else
        for i in eachindex(problems)
            results[i] = search(o, problems[i])
        end
    end
    for (s, (id, elapsed, cert, limited)) in zip(pending, results)
        o.searches += 1
        if id
            o.identifiable_time += elapsed
            cert === nothing || push!(o.certificates, cert)
        else
            o.nonidentifiable_time += elapsed
            limited && push!(o.undecided, copy(s))
        end
        record!(o, s, id)
    end
    return o
end

function isidentifiable!(o::MonotoneOracle, s::BitSet)
    check!(o, [s])
    return known(o, s)::Bool
end

with_element(s::BitSet, c::Int) = push!(copy(s), c)

# Extend the identifiable set `s` to a maximal one by adding the candidates in
# order whenever the query stays identifiable. (Trying candidates in halving
# blocks gives the same result but needed more searches on examples.) With
# several threads, first check each candidate added to `s` alone: a candidate
# that cannot be added now cannot be added later either.
function extend_greedily!(o::MonotoneOracle, s::BitSet)
    all_candidates = BitSet(eachindex(o.candidates))
    isidentifiable!(o, all_candidates) && return all_candidates
    rest = [c for c in eachindex(o.candidates) if !(c in s)]
    o.parallel && check!(o, [with_element(s, c) for c in rest])
    for c in rest
        t = with_element(s, c)
        isidentifiable!(o, t) && (s = t)
    end
    return s
end

# The index of the first identifiable set in `sets`, or `nothing`. With
# several threads, the undecided sets are checked a batch at a time.
function first_identifiable!(o::MonotoneOracle, sets::Vector{BitSet})
    k = o.parallel ? Threads.nthreads() : 1
    n = length(sets)
    i = 1
    while i <= n
        batch = Int[]   # undecided sets from i on, up to a known identifiable one
        j = i
        while j <= n && length(batch) < k
            answer = known(o, sets[j])
            answer === true && break
            answer === nothing && push!(batch, j)
            j += 1
        end
        hit = j <= n && known(o, sets[j]) === true
        check!(o, sets[batch])
        for b in i:(j - 1)
            known(o, sets[b]) === true && return b
        end
        hit && return j
        i = j
    end
    return nothing
end

# Update the minimal transversals (hitting sets) of a hypergraph after adding
# the edge `e` (one step of Berge's algorithm).
function add_edge_transversals(transversals::Vector{BitSet}, e::BitSet)
    new = BitSet[]
    for t in transversals
        if isempty(intersect(t, e))
            for v in e
                push!(new, with_element(t, v))
            end
        else
            push!(new, t)
        end
    end
    unique!(new)
    return [t for t in new if !any(u -> u != t && u ⊆ t, new)]
end

# All maximal identifiable sets by dualize and advance (Gunopulos et al.,
# 2003): a set that hits the complement of every maximal set found so far is
# either identifiable, and then extends to a new maximal set, or a minimal
# non-identifiable set. When every minimal hitting set is non-identifiable,
# all maximal sets have been found, and the hitting sets are exactly the
# minimal non-identifiable sets. Each maximal and each minimal
# non-identifiable set is searched once, plus the searches of the greedy
# extensions, which monotonicity often answers. Returns the maximal and the
# minimal non-identifiable sets.
function dualize_and_advance!(o::MonotoneOracle)
    universe = BitSet(eachindex(o.candidates))
    maximal = [extend_greedily!(o, BitSet())]
    transversals = add_edge_transversals([BitSet()], setdiff(universe, maximal[1]))
    while true
        i = first_identifiable!(o, transversals)
        i === nothing && return maximal, transversals
        m = extend_greedily!(o, transversals[i])
        push!(maximal, m)
        transversals = add_edge_transversals(transversals, setdiff(universe, m))
    end
end

const MAXIMAL_SEARCH_KEYWORDS = (:heuristic, :improve, :rules, :time_limit, :term_limit)

function search_options(; kwargs...)
    check_keywords(kwargs, MAXIMAL_SEARCH_KEYWORDS,
                   (MODEL_KEYWORDS..., :candidates, :formulas, MAXIMAL_SEARCH_KEYWORDS...))
    return search_options_(; kwargs...)
end

function search_options_(; heuristic::Bool = false, improve::Bool = true, rules = nothing,
                         time_limit::Real = Inf, term_limit::Integer = typemax(Int))
    time_limit >= 0 || throw(ArgumentError("`time_limit` must be nonnegative, got $time_limit; Inf means no limit."))
    term_limit >= 0 || throw(ArgumentError("`term_limit` must be nonnegative, got $term_limit."))
    return (; heuristic, improve, rules, time_limit, term_limit)
end

"""
    MaximalGraphs{C}

The result of [`maximal_identifiable_graphs`](@ref): the maximal identifiable
sets of candidates of type `C`. It is a vector of these sets, and also has the
fields

- `maximal::Vector{Vector{C}}`: the maximal sets of candidates that can be
  added while the query stays identifiable;
- `minimal_nonidentifiable::Vector{Vector{C}}`: the minimal sets of candidates
  that make the query non-identifiable. Every set that contains none of them
  (and none of the `undecided` sets) is identifiable, and no other set is;
- `undecided::Vector{Vector{C}}`: the minimal sets whose search stopped at
  `time_limit` or `term_limit`. They count as non-identifiable, so if there
  are any, the maximal sets are identifiable but may not be maximal, and some
  minimal non-identifiable sets may be missing;
- `formulas::Vector{String}`: the identifying formula with each maximal set
  added, if `formulas = true` was given, and empty otherwise;
- `candidates::Vector{C}` and `problem::DosearchProblem`: the inputs;
- `searches::Int`: the number of searches that were needed.

Each set lists its candidates in the order of `candidates`.
"""
struct MaximalGraphs{C} <: AbstractVector{Vector{C}}
    maximal::Vector{Vector{C}}
    minimal_nonidentifiable::Vector{Vector{C}}
    undecided::Vector{Vector{C}}
    formulas::Vector{String}
    candidates::Vector{C}
    problem::DosearchProblem
    searches::Int
end

Base.size(r::MaximalGraphs) = size(r.maximal)
Base.IndexStyle(::Type{<:MaximalGraphs}) = IndexLinear()
Base.getindex(r::MaximalGraphs, i::Int) = r.maximal[i]

function show_sets(io::IO, sets, limit::Int)
    for s in first(sets, limit)
        print(io, "\n  ")
        isempty(s) ? print(io, "(none)") :
            join(io, (sprint(show, MIME"text/plain"(), c) for c in s), "; ")
    end
    length(sets) > limit && print(io, "\n  ⋮ (", length(sets) - limit, " more)")
end

plural(n, word) = string(n, " ", word, n == 1 ? "" : "s")

function Base.show(io::IO, ::MIME"text/plain", r::MaximalGraphs)
    print(io, plural(length(r.maximal), "maximal identifiable set"), " of ",
          plural(length(r.candidates), "candidate"), ":")
    show_sets(io, r.maximal, 10)
    print(io, "\n", plural(length(r.minimal_nonidentifiable), "minimal non-identifiable set"), ":")
    show_sets(io, r.minimal_nonidentifiable, 10)
    if !isempty(r.undecided)
        print(io, "\n", plural(length(r.undecided), "minimal undecided set"), " (search limit reached):")
        show_sets(io, r.undecided, 10)
    end
end

"""
    maximal_identifiable_graphs(problem; candidates, formulas = false, heuristic,
                                improve, rules, time_limit, term_limit) -> MaximalGraphs
    maximal_identifiable_graphs(data, query, graph; kwargs...)

All maximal sets of `candidates` that can be added to the model of `problem`
while the query stays identifiable by [`dosearch`](@ref), and all minimal sets
that make it non-identifiable. The result is a vector of the maximal sets; see
[`MaximalGraphs`](@ref) for its other fields. Build the problem with a set
added with [`with_candidates`](@ref).

`candidates` defaults to [`confounding_candidates`](@ref): every latent
confounder between a pair of variables that can be added. Candidates of other
kinds, such as [`LatentCause`](@ref), can be given instead. If the query is
not identifiable in the original model there are no maximal sets, and the
empty set is the only minimal non-identifiable set. With `formulas = true`,
the identifying formula of each maximal graph is computed too. `heuristic`,
`improve` and `rules` are passed to `dosearch`. The first maximal set is the
one that [`maximal_identifiable_graph`](@ref) finds.

For an LDAG, a candidate adds a latent node with edges into its variables,
which may be labeled (see [`LatentCause`](@ref)). The variables in the labels
of the candidates are context variables in every search, including the one of
the original model, so that larger sets of candidates cannot gain case
distinctions; `problem` in the result is the problem with them.

`time_limit` (seconds) and `term_limit` (derived terms) limit each search, as
in `dosearch`. A set whose search stops at a limit counts as
non-identifiable and is reported as `undecided`. Every reported maximal set is
still identifiable, but it may not be maximal. With `term_limit` the result is
reproducible; with `time_limit` it depends on the machine, its load and the
number of threads. The formulas are computed without the limits.

Limits pay off when some searches cost far more than the others. A limit that
also stops searches that would identify the query makes the result coarse,
since every superset of an undecided set counts as non-identifiable too.
Adding confounding shrinks the space of derivable terms, so the search in the
original model is a guide: set the limit above its `terms` or `time` (from
`dosearch(problem)`).

The sets are found by dualize and advance (Gunopulos et al., 2003), which
checks each maximal set and each minimal non-identifiable set of candidates
once, plus the checks needed to extend sets to maximal ones. It relies on
adding candidates never making a query identifiable, which holds for latent
confounders: they only remove d-separations. A set is known to be identifiable
without a search when all d-separations used by an earlier derivation of the
query still hold, so most searches are of non-identifiable sets. Searches run
on several threads when Julia has them (e.g. `julia --threads=auto`); the
result does not depend on the number of threads.

# Examples

In the chain `x -> z -> y`, the effect of `x` on `y` stays identifiable if
`x` and `y` are confounded (by the front-door formula) or if `z` and `y` are
confounded, but not both:

```jldoctest
julia> graph = \"\"\"
           x -> z
           z -> y
       \"\"\";

julia> maximal_identifiable_graphs("p(x,y,z)", "p(y|do(x))", graph)
2 maximal identifiable sets of 3 candidates:
  y <-> z
  x <-> y
2 minimal non-identifiable sets:
  x <-> z
  y <-> z; x <-> y
```

Confounding `x` and `z` alone, or `y` and `z` together with `x` and `y`,
makes the effect non-identifiable.
"""
function maximal_identifiable_graphs(problem::DosearchProblem;
                                     candidates::AbstractVector = confounding_candidates(problem),
                                     formulas::Bool = false, kwargs...)
    options = search_options(; kwargs...)
    problem = with_label_context(problem, candidates)
    o = MonotoneOracle(problem, collect(candidates), options)
    C = eltype(o.candidates)
    pick(s) = o.candidates[collect(s)]
    if !isidentifiable!(o, BitSet())
        empty = [C[]]
        limited = !isempty(o.undecided)
        return MaximalGraphs(Vector{C}[], limited ? Vector{C}[] : empty, limited ? empty : Vector{C}[],
                             String[], o.candidates, problem, o.searches)
    end
    maximal, minimal_nonid = dualize_and_advance!(o)
    # Every minimal non-identifiable set was searched; the undecided ones hit a limit.
    undecided = filter(s -> s in o.undecided, minimal_nonid)
    proven = filter(s -> !(s in o.undecided), minimal_nonid)
    sets = [pick(s) for s in maximal]
    unlimited = merge(options, (; time_limit = Inf, term_limit = typemax(Int)))
    f = formulas ? [something(dosearch(with_candidates(problem, s); unlimited...).formula) for s in sets] :
                   String[]
    return MaximalGraphs(sets, [pick(s) for s in proven], [pick(s) for s in undecided], f,
                         o.candidates, problem, o.searches)
end

"""
    maximal_identifiable_graph(problem; candidates, heuristic, improve, rules,
                               time_limit, term_limit) -> Union{Nothing,Vector{<:Candidate}}
    maximal_identifiable_graph(data, query, graph; kwargs...)

One maximal set of `candidates` that can be added to the model of `problem`
while the query stays identifiable, found greedily: each candidate is added in
turn if the query stays identifiable. Earlier candidates are preferred, so the
order of `candidates` decides which maximal set is found. Needs at most one
search per candidate. Returns `nothing` if the query is not identified in
the original model. See [`maximal_identifiable_graphs`](@ref) for the
arguments; with a limit, a candidate whose search stops at it is left out,
so the set may not be maximal.
"""
function maximal_identifiable_graph(problem::DosearchProblem;
                                    candidates::AbstractVector = confounding_candidates(problem),
                                    kwargs...)
    problem = with_label_context(problem, candidates)
    o = MonotoneOracle(problem, collect(candidates), search_options(; kwargs...))
    isidentifiable!(o, BitSet()) || return nothing
    return o.candidates[collect(extend_greedily!(o, BitSet()))]
end

for f in (:maximal_identifiable_graphs, :maximal_identifiable_graph)
    @eval function $f(data, query, graph::Union{AbstractString,CausalGraph};
                      transportability = nothing, selection_bias = nothing, missing_data = nothing,
                      md_sym::AbstractString = "1", context_vars = String[], warn::Bool = true, kwargs...)
        problem = DosearchProblem(data, query, graph; transportability, selection_bias,
                                  missing_data, md_sym, context_vars, warn)
        return $f(problem; kwargs...)
    end
end
