# API reference

```@meta
CurrentModule = Dosearch
```

```@docs
Dosearch
```

`using Dosearch` brings the names of normal usage: [`dosearch`](@ref),
[`DosearchProblem`](@ref), [`latex`](@ref), [`typst`](@ref),
[`maximal_identifiable_graphs`](@ref), [`maximal_identifiable_graph`](@ref),
[`BidirectedPair`](@ref), [`LatentCause`](@ref),
[`confounding_candidates`](@ref) and [`with_candidates`](@ref). The other names
below are public but not exported: write `Dosearch.CausalGraph`, or import
them with `using Dosearch: CausalGraph`.

## Identification

```@docs
dosearch
DosearchResult
isidentifiable
DosearchProblem
DagProblem
LdagProblem
```

## Graphs and distributions as values

```@docs
CausalGraph
LabeledEdge
P
```

## Derivations

```@docs
Derivation
latex
typst
```

## Rules and search strategies

```@docs
Rule
default_rules
CsiRule
default_csi_rules
SearchStrategy
FoundOrder
ScoreHeuristic
```

## Maximal identifiable graphs

```@docs
maximal_identifiable_graphs
maximal_identifiable_graph
MaximalGraphs
Candidate
BidirectedPair
LatentCause
confounding_candidates
with_candidates
```

## Internals

```@docs
VarSet
Dist
```
