# Identifying causal effects

```@meta
CurrentModule = Dosearch
DocTestSetup = :(using Dosearch)
```

[`dosearch`](@ref) takes three inputs, written as text in the syntax of the R
package: the input distributions (`data`), the target distribution (`query`)
and the causal graph (`graph`). Problems built in code can give them as Julia
values instead; see [Building problems in code](values.md).

## Distributions

A distribution is written `p(A|do(B),C)`, where `A`, `B` and `C` are comma
separated variable names. The conditioning variables `C` and the
do-operator are optional, and `p(A|C,do(B))` means the same. `data` lists one
distribution per line, as a string or a vector of strings:

```julia
data = """
    p(x,y,z)
    p(y|do(x),z)
"""
```

A distribution can also be given as a `NamedTuple` of variable roles: `0`
for the left-hand side, `1` for intervened and `2` for conditioning
variables, so `(y = 0, x = 1)` is `p(y|do(x))`.

## Graphs

The graph has one edge per line: `x -> y` for a directed edge and `x <-> y`
(or `x -- y`) for a bidirected edge, i.e. a latent common cause of `x` and
`y`. The three parts of an edge are separated by whitespace. A graph whose
edges have labels is a labeled DAG (LDAG); see [LDAGs](ldag.md).

Inputs with several lines are easiest to write as multi-line strings, as in
the examples here. The package reports an error with a hint for malformed
input, such as `x->y` without spaces or several edges on one line.

The front-door model, where a latent variable confounds `x` and `y`, but the
effect of `x` on `y` goes through `z`:

```jldoctest
julia> r = dosearch("p(x,y,z)", "p(y|do(x))", """
           x -> z
           z -> y
           x <-> y
       """)
\sum_{z}\left(p(z|x)\sum_{x}\left(p(x)p(y|z,x)\right)\right)

julia> dosearch("p(x,y)", "p(y|do(x))", """
           x -> y
           x <-> y
       """)
The query p(y|do(x)) is not identifiable by do-search.
```

## Results

[`dosearch`](@ref) returns a [`DosearchResult`](@ref). Its fields are:

- `identifiable` (also [`isidentifiable`](@ref)): whether the query was
  identified;
- `formula`: the identifying formula in LaTeX, in the notation of the R
  package;
- `derivation`: the [`Derivation`](@ref) of the formula, which
  [`latex`](@ref) and [`typst`](@ref) write out step by step (see
  [Step-by-step derivations](derivations.md));
- `limit_reached`: whether the search stopped at `time_limit` or
  `term_limit` before it could decide;
- `time` and `terms`: the seconds spent searching and the number of terms
  the search derived;
- `dot`: the derivation as a DOT graph, with `draw_derivation = true`;
- `rule_times`, `rules` and `problem`.

A query is *not identifiable by do-search* if the search derives every term
it can without reaching the query. Do-search is complete for the problems
that do-calculus is known to solve completely, such as the identification of
causal effects from a single observational distribution.

## Several data sources

The inputs can mix observational and experimental distributions, also over
different sets of variables. Here an experiment on `z` and an observational
distribution of `x` and `z` identify the effect of `x` on `y`:

```jldoctest
julia> data = """
           p(x,z)
           p(y|do(z))
       """;

julia> graph = """
           x -> z
           z -> y
           x <-> y
       """;

julia> dosearch(data, "p(y|do(x))", graph)
\sum_{z}\left(p(z|x)p(y|do(z))\right)
```

## Transportability, selection bias and missing data

Special nodes mark discrepancies between data sources:

- `transportability`: comma separated transportability nodes, which mark the
  differences between a source domain and the target domain. A distribution
  that conditions on a transportability node comes from the source domain.
- `selection_bias`: comma separated selection bias nodes. A distribution that
  conditions on a selection node is observed under selection.
- `missing_data`: comma separated mechanisms `r_x : x`, where `r_x` is the
  response indicator of `x`. The proxy variable `x*` of `x` (`x` where
  observed, missing otherwise) is added automatically. `r_x = 1` in a
  distribution means that `x` is observed.

A source domain where `s` changes the distribution of `z`, with an
experiment on `x` there:

```jldoctest
julia> data = """
           p(x,y,z)
           p(y,z|do(x),s)
       """;

julia> graph = """
           x -> z
           z -> y
           s -> z
           x <-> y
       """;

julia> dosearch(data, "p(y|do(x))", graph; transportability = "s")
\sum_{z}\left(p(z|x)p(y|do(x),z,s)\right)
```

Missing values of `x` and `y`, whose response indicators depend on `y`, with
the marginal distribution of `y` known from another source:

```jldoctest
julia> data = """
           p(x*,y*,r_x,r_y)
           p(y)
       """;

julia> graph = """
           x -> y
           y -> r_y
           r_y -> r_x
       """;

julia> dosearch(data, "p(y|do(x))", graph; missing_data = "r_x : x, r_y : y")
\frac{\left(p(y)p(x|r_x = 1,y,r_y = 1)\right)}{\sum_{y}\left(p(y)p(x|r_x = 1,y,r_y = 1)\right)}
```

## Search options

- `heuristic = true` expands the most promising terms first. It often finds
  identifiable queries much faster; a non-identifiable query costs the same
  either way.
- `rules` restricts the search to some of the rules, given as R rule codes or
  as [`Rule`](@ref)s ([`CsiRule`](@ref)s for LDAGs).
- `time_limit` (seconds) and `term_limit` (derived terms) stop a search that
  takes too long. The result then has `limit_reached = true`, and the query
  is not identified although it may be identifiable. `term_limit` gives the
  same result on any machine.
- `improve = false` searches with all rules, as the original algorithm does.

## Problems

[`DosearchProblem`](@ref) parses and validates the inputs once, and
`dosearch(problem; kwargs...)` searches. It returns a [`DagProblem`](@ref) or
an [`LdagProblem`](@ref), depending on the graph.

```jldoctest
julia> p = DosearchProblem("p(x,y,z)", "p(y|do(x))", """
           x -> y
           z -> x
           z -> y
       """)
DagProblem with 3 variables
  query: p(y|do(x))
  data:  p(y,x,z)

julia> dosearch(p; heuristic = true).identifiable
true
```
