# Step-by-step derivations

```@meta
DocTestSetup = :(using Dosearch)
```

The identifying formula alone does not show why it is right. The
`derivation` field of a result holds how the search derived it, and
[`latex`](@ref) and [`typst`](@ref) write it out for a document. The
derivation starts from the query. Each line rewrites one term and names the
rule it uses. The last line is the identifying formula.

For the front-door model, `latex(r)` gives

```jldoctest frontdoor
julia> r = dosearch("p(x,y,z)", "p(y|do(x))", """
           x -> z
           z -> y
           x <-> y
       """);

julia> print(latex(r))
\begin{align*}
p(y \mid \mathrm{do}(x)) &= \sum_{z} p(z, y \mid \mathrm{do}(x)) && \text{marginalization} \\
&= \sum_{z} p(z \mid \mathrm{do}(x)) p(y \mid \mathrm{do}(x), z) && \text{chain rule} \\
&= \sum_{z} p(z \mid x) p(y \mid \mathrm{do}(x), z) && \text{rule 2: } (z \perp\!\!\!\perp x)_{\mathcal{G}_{\underline{x}}} \\
&= \sum_{z} p(z \mid x) p(y \mid \mathrm{do}(z, x)) && \text{rule 2: } (y \perp\!\!\!\perp z \mid x)_{\mathcal{G}_{\overline{x}\,\underline{z}}} \\
&= \sum_{z} p(z \mid x) p(y \mid \mathrm{do}(z)) && \text{rule 3: } (y \perp\!\!\!\perp x \mid z)_{\mathcal{G}_{\overline{z}\,\overline{x}}} \\
&= \sum_{z} p(z \mid x) \sum_{x'} p(y, x' \mid \mathrm{do}(z)) && \text{marginalization} \\
&= \sum_{z} p(z \mid x) \sum_{x'} p(x' \mid \mathrm{do}(z)) p(y \mid \mathrm{do}(z), x') && \text{chain rule} \\
&= \sum_{z} p(z \mid x) \sum_{x'} p(x') p(y \mid \mathrm{do}(z), x') && \text{rule 3: } (x \perp\!\!\!\perp z)_{\mathcal{G}_{\overline{z}}} \\
&= \sum_{z} p(z \mid x) \sum_{x'} p(x') p(y \mid z, x') && \text{rule 2: } (y \perp\!\!\!\perp z \mid x)_{\mathcal{G}_{\underline{z}}}
\end{align*}
```

which renders as

```@eval
using Dosearch, Markdown
graph = """
    x -> z
    z -> y
    x <-> y
"""
r = dosearch("p(x,y,z)", "p(y|do(x))", graph)
Markdown.MD(Markdown.LaTeX(replace(latex(r), "align*" => "aligned")))
```

and `typst(r)` gives the same as a Typst equation:

```jldoctest frontdoor
julia> print(typst(r))
$ p(y | "do"(x)) &= sum_(z) p(z, y | "do"(x)) && quad "marginalization" \
  &= sum_(z) p(z | "do"(x)) p(y | "do"(x), z) && quad "chain rule" \
  &= sum_(z) p(z | x) p(y | "do"(x), z) && quad "rule 2:" (z ⫫ x)_(cal(G)_(underline(x))) \
  &= sum_(z) p(z | x) p(y | "do"(z, x)) && quad "rule 2:" (y ⫫ z | x)_(cal(G)_(overline(x) underline(z))) \
  &= sum_(z) p(z | x) p(y | "do"(z)) && quad "rule 3:" (y ⫫ x | z)_(cal(G)_(overline(z) overline(x))) \
  &= sum_(z) p(z | x) sum_(x') p(y, x' | "do"(z)) && quad "marginalization" \
  &= sum_(z) p(z | x) sum_(x') p(x' | "do"(z)) p(y | "do"(z), x') && quad "chain rule" \
  &= sum_(z) p(z | x) sum_(x') p(x') p(y | "do"(z), x') && quad "rule 3:" (x ⫫ z)_(cal(G)_(overline(z))) \
  &= sum_(z) p(z | x) sum_(x') p(x') p(y | z, x') && quad "rule 2:" (y ⫫ z | x)_(cal(G)_(underline(z))) $
```

In the REPL, the derivation displays as plain text:

```jldoctest frontdoor
julia> r.derivation
p(y|do(x))
  = Σ_z p(z,y|do(x))                           (marginalization)
  = Σ_z p(z|do(x)) p(y|do(x),z)                (chain rule)
  = Σ_z p(z|x) p(y|do(x),z)                    (rule 2: (z ⊥⊥ x) in G_{x̲})
  = Σ_z p(z|x) p(y|do(z,x))                    (rule 2: (y ⊥⊥ z | x) in G_{x̅,z̲})
  = Σ_z p(z|x) p(y|do(z))                      (rule 3: (y ⊥⊥ x | z) in G_{z̅,x̅})
  = Σ_z p(z|x) Σ_x' p(y,x'|do(z))              (marginalization)
  = Σ_z p(z|x) Σ_x' p(x'|do(z)) p(y|do(z),x')  (chain rule)
  = Σ_z p(z|x) Σ_x' p(x') p(y|do(z),x')        (rule 3: (x ⊥⊥ z) in G_{z̅})
  = Σ_z p(z|x) Σ_x' p(x') p(y|z,x')            (rule 2: (y ⊥⊥ z | x) in G_{z̲})
```

## Notation

- Do-calculus steps give the d-separation that licenses them in Pearl's form:
  ``(Y \perp\!\!\!\perp Z \mid X, W)_{\mathcal{G}_{\overline{X}\underline{Z}}}``
  is the d-separation in the graph without the edges into ``X`` and out of
  ``Z``. For rule 3, the second overlined set is Pearl's ``Z(W)``: the nodes
  of ``Z`` that are not ancestors of ``W``.
- LDAG steps give the context-specific independence instead, such as
  ``(Y \perp\!\!\!\perp X \mid A = 1)``.
- Steps of probability calculus are named: marginalization, conditioning,
  chain rule, cases, law of total probability and evaluation (of a function
  of ``A`` at ``A = 0``, written ``[f]_{A = 0}``). Missing data steps replace
  a proxy by its true variable where it is observed. A do-calculus step that
  intervenes on a response indicator sets it to 1, as in the R package, and
  the derivation writes out that the term it rewrites to is evaluated there,
  such as ``p(y^* \mid r_y = 1, x)``. The formula of the R package leaves
  this implicit: ``p(y^* \mid r_y, x)``.
- A summation variable that would repeat a name in use is primed, as the
  inner ``x'`` above, where the formula of the R package sums over ``x``
  again.

## Options

- `inputs = true`: the formula keeps terms that are computed from one input
  distribution by probability manipulations alone, such as ``p(z \mid x)``
  above. With `inputs = true` the derivation rewrites those too, down to the
  input distributions.
- `labels`: a dictionary from variable names to code that replaces their
  default rendering. By default, names longer than one letter are set in
  italics, the part after the first `_` becomes a subscript (rendered the same
  way, so `R_Y₁` is ``R_{Y_1}``), trailing subscript digits such as `₁` become
  a subscript and a trailing `*` (proxy variables) a superscript.

```jldoctest
julia> r = dosearch("p(x,y)", "p(y|x)", "x -> y");

julia> print(typst(r; inputs = true, labels = Dict("x" => "X", "y" => "Y")))
$ p(Y | X) &= frac(p(Y, X), sum_(Y') p(Y', X)) && quad "conditioning" $
```

## Documents

The LaTeX output is an `align*` environment of the `amsmath` package. Long
derivations need `\allowdisplaybreaks` in the preamble to break across pages.
The Typst output is a block equation; long derivations need
`#show math.equation: set block(breakable: true)`.
