# Building problems in code

```@meta
CurrentModule = Dosearch
DocTestSetup = :(using Dosearch)
```

The text syntax is the main way to write a problem: it reads like the
mathematics and matches the R package. When a program builds the problem
instead, for example a graph over many time points or parts combined in
several ways, text assembled by interpolation and concatenation is easy to get
wrong. The inputs can then be Julia values:

- a [`CausalGraph`](@ref) of edges: `"x" => "y"`,
  [`BidirectedPair`](@ref)`("x", "y")` and, in an LDAG,
  [`LabeledEdge`](@ref)`("x" => "y", "z" => 1)`;
- a distribution [`P`](@ref): `P("y"; given = "z", intervened = "x")` is
  `p(y|do(x),z)`.

They are not exported, so a program imports them:
`using Dosearch: CausalGraph, LabeledEdge, P`. The values render to the text
syntax and go through the same parser, so they give the same problems as the
text. Their constructors check names, self-loops and cycles where the values
are built, so that a name cannot change the meaning of the text. Every
function that takes a graph, input distributions or a query takes them, and
the input distributions can mix values and text.

## Generated graphs

A treatment at each visit, with a covariate that the earlier treatment
affects and that shares an unmeasured cause with the outcome:

```jldoctest values
julia> using Dosearch: CausalGraph, LabeledEdge, P

julia> visits = 0:1;

julia> g = CausalGraph(("L$k" => "A$k" for k in visits), ("A$k" => "Y" for k in visits),
                       "L0" => "Y", "A0" => "L1", "L0" => "L1", BidirectedPair("L1", "Y"))
CausalGraph with 8 edges:
  L0 -> A0
  L1 -> A1
  A0 -> Y
  A1 -> Y
  L0 -> Y
  A0 -> L1
  L0 -> L1
  L1 <-> Y

julia> data = P(["L$k" for k in visits]..., ["A$k" for k in visits]..., "Y")
p(L0, L1, A0, A1, Y)

julia> query = P("Y"; intervened = ["A$k" for k in visits])
p(Y | do(A0, A1))

julia> dosearch(data, query, g)
\sum_{L1,L0}\left(p(Y|A0,A1,L1,L0)\left(p(L0)p(L1|A0,L0)\right)\right)
```

## Combining graphs

Graphs combine with `∪`. A graph in the text syntax converts to a
`CausalGraph`, and `string` gives the text back, so the two forms mix. The
order of the edges is kept, because it fixes the order of the variables in
formulas.

```jldoctest values
julia> print(CausalGraph("x -> y") ∪ CausalGraph("z" => "y"))
x -> y
z -> y
```

## LDAGs and values

An edge of an LDAG with a label is a [`LabeledEdge`](@ref), and a variable of
a distribution is assigned a value by a pair:

```jldoctest values
julia> lg = CausalGraph("I_X" => "X", LabeledEdge("L" => "X", "I_X" => 1), "L" => "Y", "X" => "Y");

julia> dosearch(P("L", "X", "Y"), P("Y"; given = ["X", "I_X" => 1]), lg)
\sum_{L}\left(p(L)p(Y|X,L)\right)
```

## Special nodes

Transportability and selection bias nodes are collections of names, as in
`selection_bias = ["s"]`, and missing data mechanisms are pairs of a response
indicator and its variable:

```jldoctest values
julia> dosearch([P("x*", "y*", "r_x", "r_y"), P("y")], P("y"; intervened = "x"),
                CausalGraph("x" => "y", "y" => "r_y", "r_y" => "r_x");
                missing_data = ["r_x" => "x", "r_y" => "y"])
\frac{\left(p(y)p(x|r_x = 1,y,r_y = 1)\right)}{\sum_{y}\left(p(y)p(x|r_x = 1,y,r_y = 1)\right)}
```
