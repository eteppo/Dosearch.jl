# Dosearch.jl

Dosearch.jl identifies causal effects from multiple, possibly incomplete data
sources. Given input distributions (observational or experimental), a causal
graph and a target distribution (the query), it searches for a way to
express the query in terms of the inputs with the rules of do-calculus and
probability theory. If the search succeeds, it returns the identifying
formula.

It is a Julia port of the R package
[dosearch](https://github.com/santikka/dosearch) (Tikka, Hyttinen and
Karvanen, 2021), and gives the same formulas as the R package. Beyond the port, it

- writes the derivation of a formula step by step in LaTeX or Typst, ready
  for a supplement or appendix ([Step-by-step derivations](derivations.md));
- searches for the maximal amount of latent confounding that can be added to
  a graph while the query stays identifiable
  ([Maximal identifiable graphs](maximal.md));
- takes graphs and distributions as Julia values, for problems built in code
  ([Building problems in code](values.md)).

## Installation

The package is not registered yet:

```julia
using Pkg
Pkg.add(url = "https://github.com/eteppo/Dosearch.jl")
```

## A first example

The causal effect of `x` on `y` in a graph where `z` confounds them is
identified by adjusting for `z`:

```jldoctest
julia> using Dosearch

julia> r = dosearch("p(x,y,z)", "p(y|do(x))", """
           x -> y
           z -> x
           z -> y
       """);

julia> r.identifiable
true

julia> r.formula
"\\sum_{z}\\left(p(z)p(y|x,z)\\right)"

julia> r.derivation
p(y|do(x))
  = Σ_z p(y,z|do(x))             (marginalization)
  = Σ_z p(z|do(x)) p(y|do(x),z)  (chain rule)
  = Σ_z p(z) p(y|do(x),z)        (rule 3: (z ⊥⊥ x) in G_{x̅})
  = Σ_z p(z) p(y|x,z)            (rule 2: (y ⊥⊥ x | z) in G_{x̲})
```

## Contents

```@contents
Pages = ["guide.md", "tutorial.md", "ldag.md", "derivations.md", "maximal.md", "values.md", "differences.md", "api.md"]
Depth = 2
```

## References

- Tikka, S., Hyttinen, A. and Karvanen, J. (2021). Causal effect
  identification from multiple incomplete data sources: a general
  search-based approach. *Journal of Statistical Software* 99(5).
  [doi:10.18637/jss.v099.i05](https://doi.org/10.18637/jss.v099.i05)
- Tikka, S., Hyttinen, A. and Karvanen, J. (2019). Identifying causal effects
  via context-specific independence relations. *Advances in Neural
  Information Processing Systems* 32.
- Pearl, J. (1995). Causal diagrams for empirical research. *Biometrika*
  82(4), 669–688.
