# Differences from the R package

Dosearch.jl gives the same formulas as R dosearch 1.0.12, character for
character: the test suite compares against the R package on its examples
and tests and on several hundred random problems. The rule order, the order in
which subsets of variables are tried and the search order all follow the R
package. The differences are these.

## Interface

- Options are keyword arguments instead of a `control` list.
- Times are in seconds (`time_limit` in R is in hours). `time_limit = 0`
  stops at once; no limit is `Inf`. `term_limit` is a deterministic
  alternative, and the result tells whether a search stopped at its limit
  (`limit_reached`).
- The DOT graph of `draw_derivation = true` is in the field `dot`. The field
  `derivation` holds the step-by-step derivation.
- `draw_all = true` draws the derivation steps in the order they were found,
  not in hash table order.
- Up to 32 variables are supported in DAGs (R: 15) and 64 in LDAGs (R: 30).

## Inputs

Dosearch.jl rejects some inputs that R accepts, with an error that explains
the syntax:

- an empty graph or no input distributions, where R reports the query as
  not identifiable;
- the characters `\n` (a backslash and `n`, as in a raw string), which are
  not a line break;
- variable names in the graph with `(`, `)`, `,`, `|` or `=`, which no
  distribution could refer to;
- unknown keyword arguments and rule codes, and negative limits.

It warns about variables of the distributions that are not in the graph,
which are unconnected and probably a typo (`warn = false` turns warnings
off). Other differences:

- The non-rule code `0` in `rules` is rejected; R applies it as a copy of
  the previous rule.
- Value assignments such as `x = 1` without missing data mechanisms are
  rejected.
- Variable names may contain `do`, as in `dose`.
- LDAG label values other than 0 and 1 are rejected; R ignores them.

## LDAGs

- An edge can be listed on several lines, one label on each, which means the
  same as the labels separated by `;` on one line.
- Edges that vanish in every context are left out of the contexts. R means to
  do this too, but keeps them because of a comparison bug. All the R fixtures
  give the same results anyway.
- The CSI criterion looks up the contexts of a statement by the variables that
  appear in labels. R also uses context variables given with `con_vars` that
  appear in no label. It then finds no contexts, and its reasoning by cases
  accepts false independences, such as that of a variable and its parent, so
  it can identify a query that is not identifiable. The two agree when the
  assigned context variables all appear in labels, as in all the R fixtures.
- A variable fixed to a value in the query or an input distribution is a
  context variable, as if it were in `con_vars`. R assigns values only to the
  variables in labels and in `con_vars`, so without `con_vars` it cannot
  identify, for example, `p(Y|C=0,I_C=1)` in the LDAG `I_C -> C`,
  `X -> C : I_C = 1`, `X -> Y`, even from `p(X,C,Y)`, although it equals
  `p(Y)`. With `con_vars`, R gives the same result as Dosearch.jl.
