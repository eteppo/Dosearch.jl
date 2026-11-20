# Dosearch.jl

Experimental. A Julia port of the R package
[dosearch](https://github.com/santikka/dosearch) (Tikka, Hyttinen & Karvanen,
2021, [doi:10.18637/jss.v099.i05](https://doi.org/10.18637/jss.v099.i05)):
identification of causal effects from multiple, possibly incomplete data
sources by a search over do-calculus. It gives the same formulas as the R
package, and adds step-by-step derivations in LaTeX or Typst and a search for
maximal identifiable graphs.

**Documentation:** https://eteppo.github.io/Dosearch.jl/

```julia
using Pkg
Pkg.add(url = "https://github.com/eteppo/Dosearch.jl")
```

```julia
using Dosearch

r = dosearch("p(x,y,z)", "p(y|do(x))", """
    x -> z
    z -> y
    x <-> y
""")
r.formula   # the front-door formula
```

Licensed under GPL (>= 3), as dosearch.
