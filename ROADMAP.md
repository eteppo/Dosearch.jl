# Roadmap

Planned work, roughly in order of priority.

## Features

- **Closest identifiable graph** (requested). When a query is not
  identifiable, or its search reaches a limit, search for the identifiable
  graphs closest to the user's: those that remove as few of its unobserved
  variables (bidirected edges, latent causes, latent nodes of an LDAG) as
  possible. This is the maximal identifiable graph search run from the other
  end: the graph without its confounding is the base, and its own confounding
  is the candidates, so the machinery exists. Each result says which
  confounding had to go, which tells the user what assumption identification
  needs.

- **Unmeasured causes of loss to follow-up.** With selection nodes, the
  candidates of the maximal identifiable graph search can include pairs with
  a selection node, which asks which unmeasured causes of selection the
  identification tolerates (the tutorial's problem has 14 maximal sets then).
  Check the formulas of such sets before documenting it.

