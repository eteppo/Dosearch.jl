# Port of reference/dosearch/tests/testthat/test-errors.R (DAG cases).
# Checks of R argument types, the `control` list and the igraph and dagitty
# inputs have no counterpart in Julia and are not ported.

# The error message of `f()`, or `nothing` if it does not throw.
function error_message(f)
    try
        f()
    catch e
        return sprint(showerror, e)
    end
    return nothing
end

macro test_error(msg, ex)
    return :(@test occursin($(esc(msg)), something(error_message(() -> $(esc(ex))), "no error")))
end

@testset "malformed alternative distribution format input fails" begin
    @test_error "all role values must be non-missing and finite" dosearch((x = NaN,), "p(y|do(x))", "x -> y")
    @test_error "all role values must be either 0, 1 or 2" dosearch((x = -1,), "p(y|do(x))", "x -> y")
    @test_error "at least one variable must have role value 0" dosearch((x = 1, y = 1), "p(y|do(x))", "x -> y")
end

@testset "transportability and selection bias nodes exist" begin
    @test_error "Transportability nodes t are not present in the graph" dosearch(
        "p(x)", "p(y)", "x -> y"; transportability = "t")
    @test_error "Selection bias nodes s are not present in the graph" dosearch(
        "p(x)", "p(y)", "x -> y"; selection_bias = "s")
end

@testset "graph errors" begin
    @test_error "Invalid graph, the graph is empty." dosearch("p(x)", "p(y)", "")
    # R allows at most 30 graph nodes (15 variables); Dosearch.jl allows 32 variables.
    big = join(["v$i -> v$(i + 1)" for i in 1:32], "\n")
    @test_error "more than 64 nodes" dosearch("p(v1)", "p(v2)", big)
    @test_error "Invalid graph, malformed lines found" dosearch("p(x)", "p(y)", "x - > y\nx z w y\nx y")
    @test_error "Invalid graph, unknown edge types found: edge" dosearch("p(x)", "p(y)", "x edge y")
    @test_error "the graph contains self-loops: x -> x" dosearch("p(x)", "p(y)", "x -> x")
    @test_error "the graph contains self-loops: x <-> x" dosearch("p(x)", "p(y)", "x <-> x")
    @test_error "the graph contains cycles" dosearch("p(x)", "p(y)", "x -> z\nz -> y\ny -> x")
end

@testset "malformed missing data mechanisms fail" begin
    @test_error "Malformed missing data mechanisms" dosearch("p(x)", "p(y)", "x -> y"; missing_data = "")
    @test_error "Malformed missing data mechanisms" dosearch("p(x)", "p(y)", "x -> y"; missing_data = "r_x")
    @test_error "A missing data mechanism cannot be a parent of a true variable" dosearch(
        "p(x)", "p(y)", "r_x -> y"; missing_data = "r_x : x, r_y : y")
end

@testset "syntactically incorrect data inputs fail" begin
    for m in ["(", "p(", "p(x", "p(x|y", "p(x|do(x", "p(x|do(x)", "p(x||y)", "p(x|do(y),do(z))",
              "p(x|z,do(y),w)", "p()", "p(x,)", "q(x)"]
        @test_error "Unable to parse input distribution" dosearch(m, "p(y)", "x -> y")
    end
    @test_error "Unable to parse target distribution" dosearch("p(x)", "p(y|)", "x -> y")
end

@testset "syntactically correct but semantically incorrect inputs fail" begin
    md = "r_x : x, r_y : y"
    @test_error "duplicated variables" dosearch("p(x,x)", "p(y)", "x -> y")
    @test_error "multiple symbols used for missing data mechanisms" dosearch(
        "p(x,r_x=2,r_y=1)", "p(y)", "x -> r_x"; missing_data = md)
    @test_error "invalid symbol used for a missing data mechanism" dosearch(
        "p(x,r_x=2)", "p(y)", "x -> r_x"; missing_data = md)
    @test_error "true and proxy versions of the same variable on the left-hand side" dosearch(
        "p(x,x*)", "p(y)", "x -> r_x"; missing_data = md)
    @test_error "true and proxy versions of the same variable on the right-hand side" dosearch(
        "p(y|x,x*)", "p(y)", "x -> r_x"; missing_data = md)
    @test_error "true variable of a proxy variable on the left-hand side" dosearch(
        "p(x|x*)", "p(y)", "x -> r_x"; missing_data = md)
    @test_error "proxy variable of a true variable on the left-hand side" dosearch(
        "p(x*|x)", "p(y)", "x -> r_x"; missing_data = md)
    @test_error "value assignment of a non-missing data mechanism" dosearch(
        "p(x = 1)", "p(y)", "x -> r_x"; missing_data = md)
    # R treats `x=1` as a variable name when there are no missing data
    # mechanisms; Dosearch.jl rejects it.
    @test_error "value assignments are only supported for missing data mechanisms" dosearch(
        "p(x = 1)", "p(y)", "x -> y")
    @test_error "same variable on the left and right-hand side" dosearch("p(x|x)", "p(y)", "x -> y")
    @test_error "same variable on the left and right-hand side" dosearch("p(x|do(x))", "p(y)", "x -> y")
    @test_error "are both conditioned on and intervened on" dosearch("p(y|x,do(x))", "p(y)", "x -> y")
    @test_error "transportability node on the left-hand side" dosearch(
        "p(t)", "p(y)", "t -> y"; transportability = "t")
    @test_error "a transportability node cannot be a child of another node" dosearch(
        "p(x)", "p(y)", "y -> t"; transportability = "t")
    @test_error "intervention on a transportability node" dosearch(
        "p(y|do(t))", "p(y)", "t -> y"; transportability = "t")
    @test_error "selection bias node on the left-hand side" dosearch(
        "p(s)", "p(y)", "y -> s"; selection_bias = "s")
    @test_error "a selection bias node cannot be a parent of another node" dosearch(
        "p(x)", "p(y)", "s -> y"; selection_bias = "s")
    @test_error "intervention on a selection bias node" dosearch(
        "p(y|do(s))", "p(y)", "y -> s"; selection_bias = "s")
    @test_error "Invalid query p(t)" dosearch("p(y)", "p(t)", "t -> y"; transportability = "t")
end

@testset "variable names" begin
    # R rewrites every "do" in a distribution, so it cannot parse names such
    # as `dose`; Dosearch.jl only treats `do(` as the do-operator.
    out = dosearch("p(dose,y)", "p(y|do(dose))", "dose -> y")
    @test out.formula == "p(y|dose)"
end

# Inputs that R accepts silently or rejects with a terse message. Dosearch.jl
# rejects them with an explanation of the syntax.
@testset "helpful errors for unusual inputs" begin
    # A backslash followed by `n` (from raw strings or text copied from code)
    # is not a line break.
    @test_error "The graph contains the characters `\\n`, which are not a line break" dosearch(
        "p(x,y)", "p(y|do(x))", raw"x -> y\nz -> y")
    @test_error "The data contains the characters `\\n`" dosearch(raw"p(x,y)\np(z)", "p(y|do(x))", "x -> y")
    @test_error "The query contains the characters `\\n`" dosearch("p(x,y)", raw"p(y|do(x))\n", "x -> y")
    @test_error "The graph contains the characters `\\n`" dosearch(
        "p(x,y)", "p(y)", raw"x -> y : z = 1\nz -> y")
    # empty inputs
    @test_error "Invalid graph, the graph is empty. Write one edge per line" dosearch("p(x)", "p(y)", " \n ")
    @test_error "No input distributions" dosearch("", "p(y)", "x -> y")
    @test_error "No input distributions" dosearch(String[], "p(y)", "x -> y")
    @test_error "No input distributions" dosearch("", "p(y)", "x -> y : z = 1\nz -> y")
    @test_error "The query is empty" dosearch("p(x)", " ", "x -> y")
    # graph syntax
    for line in ["x->y", "x -> y; y -> z", "x -> y, z -> y", "x -> y -> z"]
        @test_error "Write one edge per line, as `x -> y` or `x <-> y`, with spaces around the arrow" dosearch(
            "p(x,y,z)", "p(y)", line)
    end
    @test_error "variable names cannot contain `(`, `)`, `,`, `|` or `=`: x|z" dosearch(
        "p(x,y)", "p(y|do(x))", "x|z -> y")
    @test_error "variable names cannot contain" dosearch("p(x,y)", "p(y|do(x))", "x -> f(y)")
    # distribution syntax
    @test_error "A distribution is written `p(y|do(x),z)`" dosearch("p(x", "p(y)", "x -> y")
    @test_error "Intervention nodes such as I_x belong to LDAGs" dosearch(
        "p(X,Y)", "p(Y|X,I_X=1)", "X -> Y\nI_X -> X")
    # keyword arguments
    @test_error "Unknown keyword argument `heurstic`. Did you mean `heuristic`?" dosearch(
        "p(x,y)", "p(y|do(x))", "x -> y"; heurstic = true)
    @test !occursin("Did you mean", error_message(() -> dosearch("p(x,y)", "p(y|do(x))", "x -> y"; foo = 1)))
    @test_error "`transportability` describes the model" dosearch(
        DosearchProblem("p(x,y)", "p(y|do(x))", "x -> y"); transportability = "s")
    @test_error "Unknown keyword argument `formula`. Did you mean `formulas`?" maximal_identifiable_graphs(
        "p(x,y)", "p(y|do(x))", "x -> y"; formula = true)
    # rules and limits
    @test_error "Unknown rule code 42. The rule codes of DAGs are -8, -7" dosearch(
        "p(x,y)", "p(y|do(x))", "x -> y"; rules = [42])
    @test_error "Invalid rule \"2\": give rules as integer codes" dosearch(
        "p(x,y)", "p(y|do(x))", "x -> y"; rules = ["2"])
    @test_error "`time_limit` must be nonnegative" dosearch("p(x,y)", "p(y|do(x))", "x -> y"; time_limit = -1)
    @test_error "`term_limit` must be nonnegative" dosearch("p(x,y)", "p(y|do(x))", "x -> y"; term_limit = -1)
    @test_error "`time_limit` must be nonnegative" maximal_identifiable_graphs(
        "p(x,y)", "p(y|do(x))", "x -> y"; time_limit = -1)
end

@testset "warnings for probable typos" begin
    # A variable of the distributions that is not in the graph is unconnected.
    @test_logs (:warn, r"Variables w of the distributions are not in the graph") dosearch(
        "p(x,y)", "p(w|do(x))", "x -> y")
    @test_logs (:warn, r"Variables w of the distributions are not in the graph") dosearch(
        "p(x,y,w)", "p(y|do(x))", "x -> y")
    @test_logs dosearch("p(x,y,w)", "p(y|do(x))", "x -> y"; warn = false)
end
