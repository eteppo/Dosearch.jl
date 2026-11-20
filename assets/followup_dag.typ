// The causal graph of the follow-up study in the tutorial (two intervals),
// with loss to follow-up C, the latent access to care V and frailty U.
// Compile with `typst compile followup_dag.typ followup_dag.svg`.
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge
#set page(width: auto, height: auto, margin: 8pt, fill: white)
#set text(size: 13pt)

#let treat = (fill: rgb("#dbe9f6"), stroke: rgb("#3b6ea8"))
#let dead = (fill: rgb("#f8dcd6"), stroke: rgb("#b5432e"))
#let latent = (stroke: (paint: rgb("#777777"), dash: "dashed"))
#let hidden = (stroke: (paint: rgb("#777777"), dash: "dashed", thickness: 0.7pt))

#diagram(
  spacing: (16mm, 13mm),
  node-stroke: 0.6pt,
  node-shape: circle,
  edge-stroke: 0.7pt,
  mark-scale: 80%,

  // columns: visit 0, interval 1, visit 1, interval 2
  // rows: severity, treatment, death, loss to follow-up
  node((0, 0), $L_0$, name: <L0>),
  node((0, 1), $A_0$, name: <A0>, ..treat),
  node((1.5, 2), $Y_1$, name: <Y1>, ..dead),
  node((3, 0), $L_1$, name: <L1>),
  node((3, 1), $A_1$, name: <A1>, ..treat),
  node((4.5, 2), $Y_2$, name: <Y2>, ..dead),
  node((1.5, 3), $C_1$, name: <C1>, ..treat),
  node((4.5, 3), $C_2$, name: <C2>, ..treat),
  node((1.5, -1), $V$, name: <V>, ..latent),
  node((6, 0.6), $U$, name: <U>, ..latent),

  edge(<L0>, <A0>, "-|>"),
  edge(<L0>, <Y1>, "-|>"),
  edge(<A0>, <Y1>, "-|>"),
  edge(<L0>, <L1>, "-|>"),
  edge(<A0>, <L1>, "-|>"),
  edge(<Y1>, <L1>, "-|>"),
  edge(<A0>, <A1>, "-|>"),
  edge(<L1>, <A1>, "-|>"),
  edge(<Y1>, <A1>, "-|>"),
  edge(<L0>, <Y2>, "-|>", bend: 42deg),
  edge(<L1>, <Y2>, "-|>"),
  edge(<A1>, <Y2>, "-|>"),
  edge(<Y1>, <Y2>, "-|>"),
  edge(<L0>, <C1>, "-|>", bend: -55deg),
  edge(<A0>, <C1>, "-|>"),
  edge(<Y1>, <C2>, "-|>"),
  edge(<L1>, <C2>, "-|>", bend: -42deg),
  edge(<A1>, <C2>, "-|>"),

  edge(<V>, <L0>, "-|>", ..hidden),
  edge(<V>, <A0>, "-|>", ..hidden),
  edge(<V>, <A1>, "-|>", ..hidden),
  edge(<U>, <Y1>, "-|>", ..hidden, bend: 12deg),
  edge(<U>, <Y2>, "-|>", ..hidden),
  edge(<U>, <L1>, "-|>", ..hidden),
)
