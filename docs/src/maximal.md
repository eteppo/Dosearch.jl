# Maximal identifiable graphs

```@meta
CurrentModule = Dosearch
DocTestSetup = :(using Dosearch)
```

How much unmeasured confounding can a causal model contain before the query
stops being identifiable? [`maximal_identifiable_graphs`](@ref) answers this
for a set of candidate additions to the graph, such as latent confounders. It
finds every maximal set of candidates that can be added while the query stays
identifiable, and every minimal set that makes it non-identifiable. Every
set that contains none of the minimal sets is identifiable, and no other set
is.

In the chain `x -> z -> y`, the effect of `x` on `y` stays identifiable if `x`
and `y` are confounded (by the front-door formula) or if `z` and `y` are
confounded, but not both:

```jldoctest
julia> graph = """
           x -> z
           z -> y
       """;

julia> r = maximal_identifiable_graphs("p(x,y,z)", "p(y|do(x))", graph; formulas = true)
2 maximal identifiable sets of 3 candidates:
  y <-> z
  x <-> y
2 minimal non-identifiable sets:
  x <-> z
  y <-> z; x <-> y

julia> r.formulas[2]
"\\sum_{z}\\left(p(z|x)\\sum_{x}\\left(p(x)p(y|z,x)\\right)\\right)"
```

The result is a vector of the maximal sets, with further fields; see
[`MaximalGraphs`](@ref). [`maximal_identifiable_graph`](@ref) finds one
maximal set greedily, with at most one search per candidate.

## Candidates

Candidates are [`Candidate`](@ref)s. By default they are
[`confounding_candidates`](@ref): a [`BidirectedPair`](@ref) for every pair
of variables, a latent common cause of the two. A [`LatentCause`](@ref) of
several variables adds one latent common cause of all of them.
[`with_candidates`](@ref) builds the problem with a set of candidates added.

## How it works

Adding confounding only adds paths to the graph, so it only removes
independences, and every derivation in the larger model is valid in the
smaller one. The identifiable sets of candidates are therefore closed under
taking subsets. The search uses this:

- It finds all maximal sets by *dualize and advance* (Gunopulos et al., 2003):
  each maximal and each minimal non-identifiable set is searched once, plus
  the searches that extend sets to maximal ones.
- It answers from earlier searches whenever monotonicity decides the answer.
- A set is identifiable without a search if every independence that an
  earlier derivation of the query used still holds.
- Searches run on all threads (`julia --threads=auto`). The result does not
  depend on the number of threads.

Most of the time goes to the searches that end non-identifiable, which have
to derive every term they can.

## Limits

`time_limit` (seconds) and `term_limit` (derived terms) stop each search
that takes too long. A set whose search stops counts as non-identifiable and
is listed in `undecided`. The reported maximal sets are then still
identifiable but may not be maximal. With `term_limit` the result is the
same on any machine.

A limit pays off when some searches cost far more than the others. A limit
that also stops searches that would identify the query makes the result
coarse, since every superset of an undecided set counts as non-identifiable
too. The search of the original model is a guide: adding confounding shrinks
the space of derivable terms, so set the limit above the `terms` (or `time`)
of `dosearch(problem)`.

## LDAGs

In an LDAG, which has no bidirected edges, a candidate adds a latent node
with edges into its variables. The edges can be labeled, which makes the
confounding context-specific. Edges into an intervened variable vanish under
the intervention, as the other edges into it do.

```jldoctest
julia> p = DosearchProblem("p(A,X,Y)", "p(Y|X,I_X=1)", """
           I_X -> X
           X -> Y : A = 1
           A -> Y
           A -> X : I_X = 1
       """);

julia> specific = LatentCause("X", "Y"; labels = ["X" => "A = 0"]);

julia> maximal_identifiable_graphs(p; candidates = [specific, BidirectedPair("X", "Y")])
1 maximal identifiable set of 2 candidates:
  latent cause of X, Y (into X : A = 0)
1 minimal non-identifiable set:
  X <-> Y
```

Confounding of `X` and `Y` keeps the effect identifiable if it vanishes
where `A = 0`, since `X` has no effect on `Y` where `A = 1`. The variables in
the labels of candidates are context variables in every search, so that
larger sets of candidates cannot gain case distinctions.
