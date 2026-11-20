# Generate the reference results of the R package dosearch that the Julia tests
# compare against. The output, dosearch_fixtures.toml, is committed, so the
# tests never need R.
#
# Usage, from the repository root, with dosearch 1.0.12 installed (build it
# from the reference/dosearch submodule):
#
#   LC_ALL=C Rscript test/fixtures/generate_fixtures.R
#
# The C locale matters: R sorts the endpoints of bidirected edges with the
# collation of the locale, and Dosearch.jl sorts by code point.
#
# Each case records the inputs, the control options and the R results:
# identifiability, the formula, the derivation (DOT) and, for a few cases, the
# verbose search log and the `draw_all` derivation.

suppressPackageStartupMessages(library(dosearch))
stopifnot(as.character(packageVersion("dosearch")) == "1.0.12")
invisible(Sys.setlocale("LC_COLLATE", "C"))

out_file <- file.path("test", "fixtures", "dosearch_fixtures.toml")

toml_str <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  x <- gsub("\r", "\\r", x, fixed = TRUE)
  x <- gsub("\t", "\\t", x, fixed = TRUE)
  paste0("\"", x, "\"")
}

# Format p(a|do(b),c), with the do-operator first or last.
fmt <- function(a, b = character(0), c = character(0), do_first = TRUE) {
  rhs <- character(0)
  do_part <- if (length(b) > 0) paste0("do(", paste(b, collapse = ","), ")") else character(0)
  c_part <- if (length(c) > 0) paste(c, collapse = ",") else character(0)
  rhs <- if (do_first) c(do_part, c_part) else c(c_part, do_part)
  if (length(rhs) == 0) {
    paste0("p(", paste(a, collapse = ","), ")")
  } else {
    paste0("p(", paste(a, collapse = ","), "|", paste(rhs, collapse = ","), ")")
  }
}

cases <- list()
add_case <- function(name, data, query, graph, transportability = NULL,
                     selection_bias = NULL, missing_data = NULL,
                     heuristic = FALSE, improve = TRUE, rules = integer(0),
                     md_sym = "1", verbose = FALSE, draw_all = FALSE) {
  cases[[length(cases) + 1L]] <<- list(
    name = name, data = data, query = query, graph = graph,
    transportability = transportability, selection_bias = selection_bias,
    missing_data = missing_data, heuristic = heuristic, improve = improve,
    rules = as.integer(rules), md_sym = md_sym, verbose = verbose,
    draw_all = draw_all
  )
}

# Reference cases: examples and tests of the R package -------------------------

g_backdoor <- "x -> y\nz -> x\nz -> y"
g_frontdoor <- "x -> z\nz -> y\nx <-> y"
for (h in c(FALSE, TRUE)) {
  for (im in c(TRUE, FALSE)) {
    sfx <- paste0(" (heuristic=", h, ", improve=", im, ")")
    add_case(paste0("backdoor", sfx), "p(x,y,z)", "p(y|do(x))", g_backdoor,
             heuristic = h, improve = im)
    add_case(paste0("frontdoor", sfx), "p(x,y,z)", "p(y|do(x))", g_frontdoor,
             heuristic = h, improve = im)
    add_case(paste0("bow-arc", sfx), "p(x,y)", "p(y|do(x))", "x -> y\nx <-> y",
             heuristic = h, improve = im)
  }
}
add_case("backdoor verbose", "p(x,y,z)", "p(y|do(x))", g_backdoor, verbose = TRUE)
add_case("frontdoor verbose heuristic", "p(x,y,z)", "p(y|do(x))", g_frontdoor,
         heuristic = TRUE, verbose = TRUE)
add_case("draw_all", "p(x,y,z, w)", "p(y|do(x))",
         "x -> y \n z -> x \n w -> z \n x <-> w \n w <-> y", draw_all = TRUE)

data_all_rules <- "
  p(w|do(x_2),y,x_1)
  p(y|do(x_2),z_1,z_2,x_1)
  p(x_1|do(x_2),w)
  p(z_2,x_2|do(x_1))
  p(z_1|do(x_1,y),x_2)
"
graph_all_rules <- "
  x_1 -> z_2
  x_1 -> z_1
  x_1 -> w
  z_1 -> w
  z_2 -> w
  x_2 -> w
  x_2 -> z_1
  x_2 -> z_2
  z_2 -> y
  z_1 -> y
"
add_case("all rules are needed (heuristic)", data_all_rules, "p(y,x_1|do(x_2),w)",
         graph_all_rules, heuristic = TRUE)
add_case("all rules are needed", data_all_rules, "p(y,x_1|do(x_2),w)",
         graph_all_rules)
all_rules <- c(-2, 2, -3, 3, 4, 5, -6, 6)
for (r in all_rules) {
  add_case(paste0("all rules are needed, without rule ", r), data_all_rules,
           "p(y,x_1|do(x_2),w)", graph_all_rules, rules = setdiff(all_rules, r))
}

data_tr_sb <- "\n  p(x,z,y|s)\n  p(y,z|t,do(x))\n"
graph_tr_sb <- "\n  x -> z\n  z -> y\n  x -> s\n  t -> z\n  x <-> y\n"
for (h in c(FALSE, TRUE)) {
  for (im in c(TRUE, FALSE)) {
    add_case(paste0("transportability and selection bias (heuristic=", h,
                    ", improve=", im, ")"),
             data_tr_sb, "p(y|do(x))", graph_tr_sb, transportability = "t",
             selection_bias = "s", heuristic = h, improve = im)
  }
}

g_md <- "x -> y\ny -> r_y\nr_y -> r_x"
md <- "r_x : x, r_y : y"
for (h in c(FALSE, TRUE)) {
  add_case(paste0("case-control without p(y) (heuristic=", h, ")"),
           "p(x*,y*,r_x,r_y)", "p(y|do(x))", g_md, missing_data = md, heuristic = h)
  add_case(paste0("case-control (heuristic=", h, ")"),
           "p(x*,y*,r_x,r_y)\np(y)", "p(y|do(x))", g_md, missing_data = md, heuristic = h)
  add_case(paste0("case-control, all rules (heuristic=", h, ")"),
           "p(x*,y*,r_x,r_y)\np(y)", "p(y|do(x))", "x -> z \n y -> z",
           missing_data = md, heuristic = h,
           rules = c(seq(-3, 3), 4, 5, -6, 6, -7, 7, -8, 8, 9, 10))
}
add_case("case-control verbose heuristic", "p(x*,y*,r_x,r_y)\np(y)", "p(y|do(x))",
         g_md, missing_data = md, heuristic = TRUE, verbose = TRUE)
add_case("case-control md_sym", "p(x*,y*,r_x,r_y)\np(y)", "p(y|do(x))", g_md,
         missing_data = md, md_sym = "0")
g_div <- "x -> y\ny -> r_x\nr_x -> r_y"
add_case("division rule 7", "p(r_x=1,r_y=1) \n p(r_x=1)", "p(r_y=1|r_x=1)", g_div,
         missing_data = md, rules = 7)
add_case("division rule -7", "p(r_x=1,r_y=1) \n p(r_x=1|r_y=1)", "p(r_y=1)", g_div,
         missing_data = md, rules = -7)
add_case("division rule 8", "p(r_x=1) \n p(r_x=1,r_y=1)", "p(r_y = 1|r_x = 1)", g_div,
         missing_data = md, rules = 8)
add_case("division rule -8", "p(r_x=1|r_y=1) \n p(r_x=1,r_y=1)", "p(r_y = 1)", g_div,
         missing_data = md, rules = -8)
add_case("division rules -1 and 8", "p(r_x=1|z) \n p(r_x=1,r_y=1)",
         "p(r_y = 1|r_x = 1)", g_div, missing_data = md, rules = c(-1, 8))

add_case("trivially non-identifiable", "p(x)", "p(y)", "x -> y")
add_case("trivially non-identifiable, improve=FALSE", "p(x)", "p(y)", "x -> y",
         improve = FALSE)
add_case("trivially identifiable", "p(y)", "p(y)", "x -> y")
add_case("trivially identifiable (heuristic)", "p(y)", "p(y)", "x -> y", heuristic = TRUE)
add_case("rule 1 insertion", "p(y|x)", "p(y)", "x -> z \n y -> z", rules = c(-1, 1))
add_case("rule 1 deletion", "p(y)", "p(y|x)", "x -> z \n y -> z", rules = c(-1, 1))
add_case("extra variables", "p(x,y,z,w)", "p(y|do(x))", "x -> z \n z -> y \n x <-> y")
add_case("alternative data format", "p(x)\np(y|x)\np(z|x,y)", "p(y|do(x))", g_backdoor)
add_case("upper case", "P(X,Y,Z)", "P(Y|do(X),Z)", "X -> Y\nZ -> X\nZ -> Y\nX <-> Y")

data3 <- "
  p(z_2,x_2|do(x_1))
  p(z_1|x_2,do(x_1,y))
  p(x_1|w_1,do(x_2))
  p(y|z_1,z_2,x_1,do(x_2))
  p(w|y,x_1,do(x_2))
"
graph3 <- "
  x_1 -> z_2
  x_1 -> z_1
  x_2 -> z_1
  x_2 -> z_2
  z_1 -> y
  z_2 -> y
  x_1 -> w
  x_2 -> w
  z_1 -> w
  z_2 -> w
"
add_case("example 3: multiple experiments", data3, "p(y,x_1|w,do(x_2))", graph3)
add_case("example 4: selection bias", "p(x,y,z_1,z_2|s)\np(z_1,z_2)", "p(y|do(x))",
         "x -> z_1\nz_1 -> z_2\nx -> y\ny -- z_2\nz_2 -> s", selection_bias = "s")
data5 <- "
  p(x,y,z_1,z_2)
  p(x,y,z_1|t_1,t_2,do(z_2))
  p(x,y,z_2|t_3,do(z_1))
"
graph5 <- "
  z_1 -> x
  x   -> z_2
  z_2 -> y
  z_1 <-> x
  z_1 <-> z_2
  z_1 <-> y
  t_1 -> z_1
  t_2 -> z_2
  t_3 -> y
"
add_case("example 5: transportability", data5, "p(y|do(x))", graph5,
         transportability = "t_1, t_2, t_3")
add_case("example 5: transportability (heuristic)", data5, "p(y|do(x))", graph5,
         transportability = "t_1, t_2, t_3", heuristic = TRUE)
add_case("example 6: missing data", "p(x*,y*,z*,m_x,m_y,m_z)", "p(x,y,z)",
         "z -> x\nx -> y\nx -> m_z\ny -> m_z\ny -> m_x\nz <-> y",
         missing_data = "m_x : x, m_y : y, m_z : z")
add_case("gets can be got", "p(x,y,z, w)", "p(y|do(x))",
         "x -> y \n z -> x \n w -> z \n x <-> w \n w <-> y")

# Random cases -----------------------------------------------------------------

set.seed(20261118)
pool <- c("a", "b", "c", "d", "e", "f", "g", "h", "u", "v", "w", "x", "y", "z")

random_graph <- function(vars, p_dir, p_bi) {
  ord <- sample(vars)
  edges <- character(0)
  n <- length(vars)
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      if (runif(1) < p_dir) edges <- c(edges, paste(ord[i], "->", ord[j]))
      if (runif(1) < p_bi) {
        arrow <- if (runif(1) < 0.2) "--" else "<->"
        pair <- sample(c(ord[i], ord[j]))
        edges <- c(edges, paste(pair[1], arrow, pair[2]))
      }
    }
  }
  if (length(edges) == 0L) edges <- paste(ord[1], "->", ord[2])
  sample(edges)
}

random_query <- function(vars) {
  y <- sample(vars, sample(1:min(2, length(vars) - 1L), 1))
  rest <- setdiff(vars, y)
  x <- sample(rest, sample(1:min(2, length(rest)), 1))
  rest <- setdiff(rest, x)
  w <- if (length(rest) > 0L && runif(1) < 0.25) sample(rest, 1) else character(0)
  fmt(y, x, w, do_first = runif(1) < 0.5)
}

random_term <- function(vars) {
  a <- sample(vars, sample(seq_along(vars), 1))
  rest <- setdiff(vars, a)
  b <- character(0)
  c <- character(0)
  if (length(rest) > 0L) {
    b <- if (runif(1) < 0.6) sample(rest, sample(0:min(2, length(rest)), 1)) else character(0)
    rest <- setdiff(rest, b)
    if (length(rest) > 0L && runif(1) < 0.4) c <- sample(rest, 1)
  }
  fmt(a, b, c, do_first = runif(1) < 0.5)
}

random_data <- function(vars) {
  style <- sample(1:3, 1)
  if (style == 1L) {
    fmt(sample(vars))
  } else if (style == 2L) {
    obs <- sample(vars, length(vars) - 1L)
    x <- sample(vars, 1)
    paste(fmt(obs), fmt(sample(setdiff(vars, x)), x), sep = "\n")
  } else {
    paste(replicate(sample(2:3, 1), random_term(vars)), collapse = "\n")
  }
}

for (i in seq_len(260)) {
  vars <- sample(pool, sample(3:6, 1))
  graph <- paste(random_graph(vars, 0.45, 0.2), collapse = "\n")
  add_case(paste("random dag", i), random_data(vars), random_query(vars), graph,
           heuristic = runif(1) < 0.5, improve = runif(1) < 0.75)
}

for (i in seq_len(80)) {
  vars <- sample(setdiff(pool, c("s", "t")), sample(3:5, 1))
  edges <- random_graph(vars, 0.45, 0.2)
  use_t <- runif(1) < 0.7
  use_s <- !use_t || runif(1) < 0.5
  data <- character(0)
  if (use_t) {
    edges <- c(edges, paste("t ->", sample(vars, 1)))
    x <- sample(vars, 1)
    data <- c(data, fmt(sample(setdiff(vars, x), sample(1:(length(vars) - 1L), 1)),
                        x, "t", do_first = runif(1) < 0.5))
  }
  if (use_s) {
    edges <- c(edges, paste(sample(vars, 1), "-> s"))
    data <- c(data, fmt(sample(vars), character(0), "s"))
  }
  if (runif(1) < 0.4) data <- c(data, fmt(sample(vars, sample(1:length(vars), 1))))
  add_case(paste("random transportability/selection bias", i),
           paste(sample(data), collapse = "\n"), random_query(vars),
           paste(sample(edges), collapse = "\n"),
           transportability = if (use_t) "t" else NULL,
           selection_bias = if (use_s) "s" else NULL,
           heuristic = runif(1) < 0.5, improve = runif(1) < 0.75)
}

for (i in seq_len(80)) {
  vars <- sample(c("x", "y", "z", "w"), sample(2:3, 1))
  md_vars <- sample(vars, sample(1:min(2, length(vars)), 1))
  sw <- paste0("r_", md_vars)
  edges <- random_graph(vars, 0.5, 0.15)
  for (s in sw) {
    if (runif(1) < 0.7) edges <- c(edges, paste(sample(vars, 1), "->", s))
  }
  if (length(sw) == 2L && runif(1) < 0.4) edges <- c(edges, paste(sw[1], "->", sw[2]))
  obs <- c(ifelse(vars %in% md_vars, paste0(vars, "*"), vars), sw)
  data <- fmt(sample(obs))
  if (runif(1) < 0.5) data <- paste(data, fmt(sample(vars, 1)), sep = "\n")
  query <- if (runif(1) < 0.5) random_query(vars) else fmt(sample(vars))
  add_case(paste("random missing data", i), data, query,
           paste(sample(edges), collapse = "\n"),
           missing_data = paste(paste(sw, ":", md_vars), collapse = ", "),
           heuristic = runif(1) < 0.5, improve = runif(1) < 0.75)
}

# Run and write -----------------------------------------------------------------

run_case <- function(cs) {
  control <- list(
    heuristic = cs$heuristic, improve = cs$improve, rules = cs$rules,
    md_sym = cs$md_sym, draw_derivation = TRUE, draw_all = cs$draw_all,
    benchmark = TRUE, warn = FALSE
  )
  log <- NULL
  if (cs$verbose) {
    control$verbose <- TRUE
    log <- capture.output(
      res <- dosearch(cs$data, cs$query, cs$graph, cs$transportability,
                      cs$selection_bias, cs$missing_data, control)
    )
  } else {
    res <- dosearch(cs$data, cs$query, cs$graph, cs$transportability,
                    cs$selection_bias, cs$missing_data, control)
  }
  list(res = res, log = log)
}

con <- file(out_file, open = "w")
writeLines(c(
  "# Generated by test/fixtures/generate_fixtures.R from R dosearch 1.0.12.",
  "# Do not edit by hand.",
  ""
), con)
n_written <- 0L
n_skipped <- 0L
for (cs in cases) {
  out <- tryCatch(run_case(cs), error = function(e) e)
  if (inherits(out, "error")) {
    n_skipped <- n_skipped + 1L
    next
  }
  res <- out$res
  if (res$time > 2000) { # milliseconds; keep the test suite fast
    n_skipped <- n_skipped + 1L
    next
  }
  lines <- c(
    "[[case]]",
    paste0("name = ", toml_str(cs$name)),
    paste0("data = ", toml_str(cs$data)),
    paste0("query = ", toml_str(cs$query)),
    paste0("graph = ", toml_str(cs$graph))
  )
  for (key in c("transportability", "selection_bias", "missing_data")) {
    if (!is.null(cs[[key]])) lines <- c(lines, paste0(key, " = ", toml_str(cs[[key]])))
  }
  lines <- c(
    lines,
    paste0("heuristic = ", tolower(cs$heuristic)),
    paste0("improve = ", tolower(cs$improve)),
    paste0("draw_all = ", tolower(cs$draw_all))
  )
  if (length(cs$rules) > 0L) {
    lines <- c(lines, paste0("rules = [", paste(cs$rules, collapse = ", "), "]"))
  }
  if (cs$md_sym != "1") lines <- c(lines, paste0("md_sym = ", toml_str(cs$md_sym)))
  lines <- c(
    lines,
    paste0("identifiable = ", tolower(res$identifiable)),
    paste0("formula = ", toml_str(res$formula))
  )
  if (nchar(res$derivation) <= 20000) {
    lines <- c(lines, paste0("derivation = ", toml_str(res$derivation)))
  }
  if (!is.null(out$log)) {
    lines <- c(lines, paste0("verbose = ", toml_str(paste(out$log, collapse = "\n"))))
  }
  writeLines(c(lines, ""), con)
  n_written <- n_written + 1L
}
close(con)
cat("Wrote", n_written, "cases to", out_file, "(skipped", n_skipped, "invalid or slow cases)\n")
