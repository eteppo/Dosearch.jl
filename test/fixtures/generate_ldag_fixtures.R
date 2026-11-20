# Generate the reference results of the R package dosearch for labeled DAGs
# (LDAGs). The output, ldag_fixtures.toml, is committed, so the tests never
# need R.
#
# Usage, from the repository root, with dosearch 1.0.12 installed:
#
#   LC_ALL=C Rscript test/fixtures/generate_ldag_fixtures.R
#
# Each case records the inputs, the control options, the R results
# (identifiability, formula, derivation and for some cases the verbose log)
# and the internal form of the problem that the R code passes to C++: the
# variables, edges, contexts, local CSIs and distributions.

suppressPackageStartupMessages(library(dosearch))
stopifnot(as.character(packageVersion("dosearch")) == "1.0.12")
invisible(Sys.setlocale("LC_COLLATE", "C"))

out_file <- file.path("test", "fixtures", "ldag_fixtures.toml")

toml_str <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  paste0("\"", x, "\"")
}
toml_ints <- function(x) paste0("[", paste(as.integer(x), collapse = ", "), "]")
toml_strs <- function(x) paste0("[", paste(vapply(x, toml_str, ""), collapse = ", "), "]")

cases <- list()
add_case <- function(name, data, query, graph, heuristic = FALSE, cache = TRUE,
                     rules = integer(0), con_vars = character(0), verbose = FALSE) {
  cases[[length(cases) + 1L]] <<- list(
    name = name, data = data, query = query, graph = graph, heuristic = heuristic,
    cache = cache, rules = as.integer(rules), con_vars = con_vars, verbose = verbose
  )
}

# The internal problem, as built by get_derivation_ldag.
parse_internal <- function(cs) {
  args <- list(
    dir_lhs = integer(0L), dir_rhs = integer(0L), vars = character(0L),
    nums = integer(0L), n = 0L, con_vars = character(0L),
    intv_vars = character(0L), parents = list(), contexts = c(),
    label_map = NULL, local_csi = list()
  )
  args <- dosearch:::transform_graph_ldag(args, cs$graph)
  args <- dosearch:::parse_labels(args)
  args <- dosearch:::parse_query_ldag(args, cs$query)
  args <- dosearch:::parse_data_ldag(args, cs$data)
  args <- dosearch:::validate_data_ldag(args)
  args <- dosearch:::validate_query_ldag(args)
  args$con_vars <- union(cs$con_vars, args$con_vars)
  args
}

# Reference cases: tests and examples of the R package -------------------------

add_case("fig 1(e)", "p(A,X,Y)", "p(Y|X,I_X=1)",
         "I_X -> X\nX -> Y : A = 1\nL -> X : I_X = 1; A = 0\nL -> Y\nA -> X : I_X = 1\nA -> Y",
         heuristic = TRUE)
add_case("fig 1(e) found order", "p(A,X,Y)", "p(Y|X,I_X=1)",
         "I_X -> X\nX -> Y : A = 1\nL -> X : I_X = 1; A = 0\nL -> Y\nA -> X : I_X = 1\nA -> Y")
add_case("fig 6(a)", "p(X,Y,W)", "p(Y|X,I_X=1)",
         "I_X -> X\nW -> X : I_X = 1\nZ -> X : I_X = 1; W = 1\nX -> Y\nZ -> Y\nZ -> X")
add_case("fig 6(b)", "p(X,Y,Z)", "p(Y|X,I_X=1)",
         "I_X -> X\nI_Z -> Z\nA -> Z : I_Z = 1\nA -> Y\nH -> X : I_X = 1\nH -> Y\nX -> Z : A = 0; I_Z = 1\nZ -> Y : A = 1")
add_case("fig 6(c)", "p(X,Y,Z)", "p(Y|X,I_X=1)",
         "X -> Y : Z = 1\nZ -> Y\nZ -> X : I_X = 1\nI_X -> X\nH -> X : I_X = 1\nH -> Z\nQ -> Z\nQ -> Y : Z = 0")
add_case("fig 6(d)", "p(X,Y,Z,A,W)", "p(Y|X,I_X=1)",
         "I_X -> X\nI_Z -> Z\nA -> W\nZ -> Y\nA -> Z\nX -> Z : I_Z = 1; A = 1\nX -> Y : A = 0\nW -> X : I_X = 1\nW -> Y : A = 0\nA -> Y\nU -> X : I_X = 1\nU -> Y : A = 1")
add_case("fig 6(e)", "p(X,Y,Z,A)", "p(Y|X,I_X=1)",
         "I_X -> X\nI_W -> W\nI_Z -> Z\nA -> W : I_W = 1\nA -> Z : I_Z = 1\nW -> Z : I_Z = 1\nZ -> X : I_X = 1\nX -> Y\nL -> W : A = 0; I_W = 1\nL -> X : I_X = 1\nM -> W : I_W = 1\nM -> Y\nN -> Z : I_Z = 1\nN -> Y")
g_nested <- "X -> Z : A = 0\nA -> Z\nA -> Y\nX -> W : B = 0\nW -> A : B = 1\nB -> A\nB -> W\nZ -> Y : A = 1"
add_case("nested csi", "p(Y)", "p(Y|X)", g_nested)
add_case("nested csi without cache", "p(Y)", "p(Y|X)", g_nested, cache = FALSE)
add_case("nested csi heuristic", "p(Y)", "p(Y|X)", g_nested, heuristic = TRUE)
add_case("trivial non-identifiability", "p(x)", "p(y)", "x -> y\nz -> y : x = 1")
add_case("trivial identifiability", "p(y)", "p(y)", "x -> y\nz -> y : x = 1")
add_case("trivial identifiability heuristic", "p(y)", "p(y)", "x -> y\nz -> y : x = 1", heuristic = TRUE)
add_case("verbose", "p(X,Y,Z)", "p(Y|X,I_X = 1)", "X -> Y\nI_X -> X\nZ -> X : I_X = 1\nZ -> Y",
         verbose = TRUE)
add_case("verbose heuristic", "p(X,Y,Z)", "p(Y|X,I_X = 1)", "X -> Y\nI_X -> X\nZ -> X : I_X = 1\nZ -> Y",
         verbose = TRUE, heuristic = TRUE)
add_case("edge vanishes", "p(X,A,Y)", "p(Y)", "X -> Y\nA -> X\nL -> X : A = 0; A = 1\nL -> Y",
         cache = FALSE)
add_case("edge vanishes with two parents", "p(X,A,Y)", "p(Y)",
         "x -> Y\nA -> X\nB -> X\nL -> X : A = 0, B = 0; A = 1, B = 0; A = 0, B = 1; A = 1, B = 1\nL -> Y",
         cache = FALSE)
add_case("derivation", "p(X,Y,Z,A,W)", "p(Y|X,V=0)",
         "V -> X\nI_Z -> Z\nA -> W\nZ -> Y\nA -> Z\nX -> Z : I_Z = 1; A = 1\nX -> Y : A = 0\nW -> X : V = 0\nW -> Y : A = 0\nA -> Y\nU -> X : V = 0\nU -> Y : A = 1")
add_case("non-primitive conditioning", "p(X|W,Z) \n p(Y|X,W)", "p(X|Y,W)", "X -> Y \n W -> Y : X = 1")
add_case("local csi", "p(A,Y)", "p(Y|X,A=1)", "X -> Y : A = 1\nA -> Y")
add_case("local csi without cache", "p(A,Y)", "p(Y|X,A=1)", "X -> Y : A = 1\nA -> Y", cache = FALSE)
add_case("local csi 2", "p(X|A=1)", "p(X|Y,A=1)", "X -> Y : A = 1\nA -> Y")
add_case("case-by-case 5", "p(x,z=0) \n p(x,z=1)", "p(x,z)", "x -> y : w = 0 \n w -> y", rules = 5)
add_case("case-by-case -5", "p(x,z=1) \n p(x,z=0)", "p(x,z)", "x -> y : w = 0 \n w -> y", rules = -5)
add_case("general-by-case 6", "p(x|w) \n p(x,z=1)", "p(x,z=0)", "x -> y : z = 0 \n z -> y", rules = c(-3, 6))
add_case("general-by-case -7", "p(x|w) \n p(x,z=1)", "p(x,z=0)", "x -> y : z = 0 \n z -> y", rules = c(-3, -7))
add_case("general-by-case 7", "p(x|w) \n p(x,z=0)", "p(x,z=1)", "x -> y : z = 0 \n z -> y", rules = c(-3, 7))
add_case("reverse product rule", "p(a|b) \n p(b)", "p(a,b)",
         "a -> x : b = 0\nc -> x : a = 0\nb -> x : c = 0", rules = -2)
add_case("custom context variables",
         "p(rx, ry)\np(x, y, rx = 1, ry = 1)\np(x, rx = 1, ry = 0)\np(y, rx = 0, ry = 1)",
         "p(y | x, rx = 1, ry = 0)", "x -> ry\nx -> y\nrx -> ry", con_vars = c("ry", "rx"))
add_case("custom context variables absent",
         "p(rx, ry)\np(x, y, rx = 1, ry = 1)\np(x, rx = 1, ry = 0)\np(y, rx = 0, ry = 1)",
         "p(y | x, rx = 1, ry = 0)", "x -> ry\nx -> y\nrx -> ry")
add_case("example 7", "P(X,Y,Z)", "P(Y|X,I_X=1)",
         "X -> Y : Z = 1\nZ -> Y\nZ -> X : I_X = 1\nI_X -> X\nH -> X : I_X = 1\nH -> Z\nQ -> Z\nQ -> Y : Z = 0")
add_case("example 8", "P(X,Y,Z,A,W)", "P(Y|X,I_X=1)",
         "I_X -> X\nI_Z -> Z\nA -> W\nZ -> Y\nA -> Z\nX -> Z : I_Z = 1; A = 1\nX -> Y : A = 0\nW -> X : I_X = 1\nW -> Y : A = 0\nA -> Y\nU -> X : I_X = 1\nU -> Y : A = 1")

# Random cases -----------------------------------------------------------------

set.seed(20261119)
pool <- c("A", "B", "C", "D", "W", "X", "Y", "Z")

random_ldag <- function() {
  n <- sample(3:5, 1)
  vars <- sample(pool, n)
  latent <- if (runif(1) < 0.4) "L" else character(0)
  all_vars <- c(vars, latent)
  ord <- sample(all_vars)
  if (length(latent) > 0L) ord <- c("L", setdiff(ord, "L"))   # latents are roots
  edges <- list()
  for (i in seq_len(length(ord) - 1L)) {
    for (j in seq.int(i + 1L, length(ord))) {
      if (runif(1) < 0.45) edges[[length(edges) + 1L]] <- c(ord[i], ord[j])
    }
  }
  if (length(edges) == 0L) edges[[1L]] <- c(ord[1L], ord[2L])
  # intervention nodes for one or two observed variables
  targets <- sample(vars, sample(1:min(2, n), 1))
  for (t in targets) edges[[length(edges) + 1L]] <- c(paste0("I_", t), t)
  parents <- function(v) unique(vapply(Filter(function(e) e[2L] == v, edges), `[`, "", 1L))
  lines <- character(0)
  for (e in sample(edges)) {
    line <- paste(e[1L], "->", e[2L])
    others <- setdiff(parents(e[2L]), e[1L])
    if (length(others) > 0L && !startsWith(e[1L], "I_") && runif(1) < 0.45) {
      groups <- replicate(sample(1:2, 1), {
        lab <- sample(others, sample(1:min(2, length(others)), 1))
        paste(paste(lab, "=", sample(0:1, length(lab), replace = TRUE)), collapse = ", ")
      })
      line <- paste(line, ":", paste(unique(groups), collapse = "; "))
    }
    lines <- c(lines, line)
  }
  list(vars = vars, targets = targets, graph = paste(lines, collapse = "\n"))
}

fmt <- function(a, b = character(0)) {
  if (length(b) == 0L) paste0("p(", paste(a, collapse = ","), ")")
  else paste0("p(", paste(a, collapse = ","), "|", paste(b, collapse = ","), ")")
}

for (i in seq_len(150)) {
  g <- random_ldag()
  vars <- g$vars
  x <- sample(g$targets, 1)
  y <- sample(setdiff(vars, x), 1)
  query <- fmt(y, c(x, paste0("I_", x, "=1")))
  data <- fmt(sample(vars))
  if (runif(1) < 0.3) {
    sub <- sample(vars, sample(2:length(vars), 1))
    data <- paste(data, fmt(sub[1], sub[-1]), sep = "\n")
  }
  add_case(paste("random ldag", i), data, query, g$graph,
           heuristic = runif(1) < 0.5, cache = runif(1) < 0.8, verbose = runif(1) < 0.1)
}

# Run and write -----------------------------------------------------------------

con <- file(out_file, open = "w")
writeLines(c(
  "# Generated by test/fixtures/generate_ldag_fixtures.R from R dosearch 1.0.12.",
  "# Do not edit by hand.",
  ""
), con)
n_written <- 0L
n_skipped <- 0L
for (cs in cases) {
  control <- list(heuristic = cs$heuristic, cache = cs$cache, rules = cs$rules,
                  con_vars = cs$con_vars, draw_derivation = TRUE, benchmark = TRUE)
  log <- NULL
  out <- tryCatch({
    if (cs$verbose) {
      control$verbose <- TRUE
      log <- capture.output(res <- dosearch:::get_derivation_ldag(cs$data, cs$query, cs$graph, control))
    } else {
      res <- dosearch:::get_derivation_ldag(cs$data, cs$query, cs$graph, control)
    }
    list(res = res, args = parse_internal(cs))
  }, error = function(e) e)
  if (inherits(out, "error") || out$res$time > 2000) {
    n_skipped <- n_skipped + 1L
    next
  }
  res <- out$res
  args <- out$args
  lines <- c(
    "[[case]]",
    paste0("name = ", toml_str(cs$name)),
    paste0("data = ", toml_str(cs$data)),
    paste0("query = ", toml_str(cs$query)),
    paste0("graph = ", toml_str(cs$graph)),
    paste0("heuristic = ", tolower(cs$heuristic)),
    paste0("cache = ", tolower(cs$cache))
  )
  if (length(cs$rules) > 0L) lines <- c(lines, paste0("rules = ", toml_ints(cs$rules)))
  if (length(cs$con_vars) > 0L) lines <- c(lines, paste0("con_vars = ", toml_strs(cs$con_vars)))
  lines <- c(
    lines,
    paste0("identifiable = ", tolower(res$identifiable)),
    paste0("formula = ", toml_str(res$formula))
  )
  if (nchar(res$derivation) <= 20000) lines <- c(lines, paste0("derivation = ", toml_str(res$derivation)))
  if (!is.null(log)) lines <- c(lines, paste0("verbose = ", toml_str(paste(log, collapse = "\n"))))
  # internal problem
  label_map <- vapply(args$label_map, function(lm) {
    ctx <- vapply(lm$contexts, function(cx) {
      paste0("[", cx$zero, ", ", cx$one, ", ", cx$equivalence, ", ",
             toml_ints(cx$from), ", ", toml_ints(cx$to), "]")
    }, "")
    paste0("[", lm$vars, ", [", paste(ctx, collapse = ", "), "]]")
  }, "")
  local_csi <- vapply(args$local_csi, function(l) toml_ints(c(l$x, l$y, l$z, l$zero, l$one)), "")
  lines <- c(
    lines,
    paste0("vars = ", toml_strs(args$vars)),
    paste0("edges = [", paste(mapply(function(a, b) toml_ints(c(args$nums[a], args$nums[b])),
                                     args$dir_lhs, args$dir_rhs), collapse = ", "), "]"),
    paste0("con_vars_internal = ", toml_strs(args$con_vars)),
    paste0("intv_vars = ", toml_strs(args$intv_vars)),
    paste0("label_map = [", paste(label_map, collapse = ", "), "]"),
    paste0("local_csi = [", paste(local_csi, collapse = ", "), "]"),
    paste0("p_list = [", paste(vapply(args$p_list, toml_ints, ""), collapse = ", "), "]"),
    paste0("q_vec = ", toml_ints(args$q_vec))
  )
  writeLines(c(lines, ""), con)
  n_written <- n_written + 1L
}
close(con)
cat("Wrote", n_written, "cases to", out_file, "(skipped", n_skipped, "invalid or slow cases)\n")
