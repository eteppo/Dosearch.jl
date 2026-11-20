# LDAGs

```@meta
CurrentModule = Dosearch
DocTestSetup = :(using Dosearch)
```

A labeled DAG (LDAG) encodes context-specific independences: an edge can
vanish in some contexts, i.e. for some values of other parents of its child.
An edge label is a comma separated list of value assignments, and the edge
vanishes in the context that matches it. An edge that vanishes in several
contexts is listed once for each:

```julia
graph = """
    X -> Y : A = 1
    Z -> Y : A = 0
    Z -> Y : B = 1, C = 0
"""
```

Here `X -> Y` vanishes where `A = 1`, and `Z -> Y` where `A = 0` or where
`B = 1` and `C = 0`. (The R package writes the alternatives of a label on one
line, separated by `;`, as in `Z -> Y : A = 0; B = 1, C = 0`. That works
too.) The labeled variables (`A`, `B` and `C`) are *context
variables*, and they are binary. So is every variable fixed to a value in the
query or an input distribution. More context variables can be named with
`context_vars`.

LDAGs have no bidirected edges: latent variables are explicit nodes that
appear in no distribution. Interventions are explicit intervention nodes
whose names start with `I_`, and the edges into an intervened variable vanish
where its intervention node is 1. A distribution conditions on the
intervention node instead of using the do-operator: `p(Y|X,I_X=1)` is the
effect of `X` on `Y`. In general, variables in distributions of an LDAG can be
assigned the values `0` and `1`.

```jldoctest ldag
julia> graph = """
           I_X -> X
           X -> Y : A = 1
           L -> X : I_X = 1
           L -> X : A = 0
           L -> Y
           A -> X : I_X = 1
           A -> Y
       """;

julia> r = dosearch("p(A,X,Y)", "p(Y|X,I_X=1)", graph)
\sum_{A}\left[\left(p(Y|X,A = 0)\left[p(A)\right]_{A = 0}\right) /\ \left(p(Y|A = 1)\left[p(A)\right]_{A = 1}\right)\right]
```

The latent `L` confounds `X` and `Y` except where `A = 0`, and `X` has no
effect on `Y` where `A = 1`, so the effect is identified case by case: the
formula `[f /\ g]` is `f` where `A = 0` and `g` where `A = 1`. The
step-by-step derivation spells this out:

```jldoctest ldag
julia> r.derivation
p(Y|X,I_X = 1)
  = Σ_A p(Y,A|X,I_X = 1)                                                  (marginalization)
  = Σ_A {p(Y,A = 0|X,I_X = 1) if A = 0; p(Y,A = 1|X,I_X = 1) if A = 1}    (cases)
  = Σ_A {p(Y|X,A = 0,I_X = 1) p(A = 0|X,I_X = 1) if A = 0; p(Y,A = 1|X,I_X = 1) if A = 1}  (chain rule)
  = Σ_A {p(Y|X,A = 0) p(A = 0|X,I_X = 1) if A = 0; p(Y,A = 1|X,I_X = 1) if A = 1}  (CSI: (Y ⊥⊥ I_X | X,A = 0))
  = Σ_A {p(Y|X,A = 0) [p(A|X,I_X = 1)]_{A = 0} if A = 0; p(Y,A = 1|X,I_X = 1) if A = 1}  (evaluation)
  = Σ_A {p(Y|X,A = 0) [p(A|I_X = 1)]_{A = 0} if A = 0; p(Y,A = 1|X,I_X = 1) if A = 1}  (CSI: (A ⊥⊥ X | I_X = 1))
  = Σ_A {p(Y|X,A = 0) [p(A)]_{A = 0} if A = 0; p(Y,A = 1|X,I_X = 1) if A = 1}  (CSI: (A ⊥⊥ I_X))
  = Σ_A {p(Y|X,A = 0) [p(A)]_{A = 0} if A = 0; p(Y|X,A = 1,I_X = 1) p(A = 1|X,I_X = 1) if A = 1}  (chain rule)
  = Σ_A {p(Y|X,A = 0) [p(A)]_{A = 0} if A = 0; p(Y|A = 1,I_X = 1) p(A = 1|X,I_X = 1) if A = 1}  (CSI: (Y ⊥⊥ X | A = 1,I_X = 1))
  = Σ_A {p(Y|X,A = 0) [p(A)]_{A = 0} if A = 0; p(Y|A = 1) p(A = 1|X,I_X = 1) if A = 1}  (CSI: (Y ⊥⊥ I_X | A = 1))
  = Σ_A {p(Y|X,A = 0) [p(A)]_{A = 0} if A = 0; p(Y|A = 1) [p(A|X,I_X = 1)]_{A = 1} if A = 1}  (evaluation)
  = Σ_A {p(Y|X,A = 0) [p(A)]_{A = 0} if A = 0; p(Y|A = 1) [p(A|I_X = 1)]_{A = 1} if A = 1}  (CSI: (A ⊥⊥ X | I_X = 1))
  = Σ_A {p(Y|X,A = 0) [p(A)]_{A = 0} if A = 0; p(Y|A = 1) [p(A)]_{A = 1} if A = 1}  (CSI: (A ⊥⊥ I_X))
```

The search for LDAGs uses its own rules ([`CsiRule`](@ref)), which add case
distinctions over context variables to the rules of probability, and the
independences of a context instead of do-calculus. Transportability,
selection bias and missing data are only supported for DAGs.
