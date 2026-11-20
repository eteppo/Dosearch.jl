# Ports of reference/dosearch/tests/testthat/test-dosearch.R and
# test-interface.R (DAG cases). The LDAG cases are not ported yet.

const BACKDOOR = "x -> y\nz -> x\nz -> y"
const FRONTDOOR = "x -> z\nz -> y\nx <-> y"

@testset "backdoor formula is identified" begin
    out = dosearch("p(x,y,z)", "p(y|do(x))", BACKDOOR)
    @test isidentifiable(out)
    @test out.formula == raw"\sum_{z}\left(p(z)p(y|x,z)\right)"
end

@testset "frontdoor formula is identified" begin
    out = dosearch("p(x,y,z)", "p(y|do(x))", FRONTDOOR)
    @test out.identifiable
    @test out.formula == raw"\sum_{z}\left(p(z|x)\sum_{x}\left(p(x)p(y|z,x)\right)\right)"
end

@testset "bow-arc is non-identifiable" begin
    out = dosearch("p(x,y)", "p(y|do(x))", "x -> y\nx <-> y"; heuristic = true)
    @test !out.identifiable
    @test out.formula === nothing
end

@testset "all rules are needed" begin
    data = """
      p(w|do(x_2),y,x_1)
      p(y|do(x_2),z_1,z_2,x_1)
      p(x_1|do(x_2),w)
      p(z_2,x_2|do(x_1))
      p(z_1|do(x_1,y),x_2)
    """
    query = "p(y,x_1|do(x_2),w)"
    graph = """
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
    """
    out = dosearch(data, query, graph; heuristic = true, draw_derivation = true)
    @test out.identifiable
    rules = [-2, 2, -3, 3, 4, 5, -6, 6]
    for r in rules
        @test !dosearch(data, query, graph; rules = setdiff(rules, r)).identifiable
    end
end

@testset "transportability and selection bias are checked" begin
    out = dosearch("p(x,z,y|s)\np(y,z|t,do(x))", "p(y|do(x))",
                   "x -> z\nz -> y\nx -> s\nt -> z\nx <-> y";
                   transportability = "t", selection_bias = "s", heuristic = true, improve = false)
    @test out.identifiable
    @test out.formula == raw"\sum_{z}\left(p(y|do(x),z,t)\sum_{y}p(z,y|x,s)\right)"
end

@testset "missing data mechanisms are checked" begin
    graph = "x -> y\ny -> r_y\nr_y -> r_x"
    md = "r_x : x, r_y : y"
    @test !dosearch("p(x*,y*,r_x,r_y)", "p(y|do(x))", graph; missing_data = md).identifiable
    data = "p(x*,y*,r_x,r_y)\np(y)"
    out = dosearch(data, "p(y|do(x))", graph; missing_data = md, heuristic = true, draw_derivation = true)
    @test out.identifiable
    out = dosearch(data, "p(y|do(x))", graph; missing_data = md)
    @test out.formula == raw"\frac{\left(p(y)p(x|r_x = 1,y,r_y = 1)\right)}" *
                         raw"{\sum_{y}\left(p(y)p(x|r_x = 1,y,r_y = 1)\right)}"
end

@testset "trivial non-identifiability is checked" begin
    out = dosearch("p(x)", "p(y)", "x -> y")
    @test !out.identifiable
    @test out.formula === nothing
end

@testset "trivial identifiability is checked" begin
    @test dosearch("p(y)", "p(y)", "x -> y").formula == "p(y)"
    @test dosearch("p(y)", "p(y)", "x -> y"; heuristic = true).formula == "p(y)"
end

@testset "verbose search works" begin
    io = IOBuffer()
    dosearch("p(x,y,z)", "p(y|do(x))", BACKDOOR; verbose = true, io)
    out = split(chomp(String(take!(io))), '\n')
    @test startswith(out[1], "Setting target")
    @test startswith(out[2], "Adding known distribution")
    @test startswith(out[3], "Initializing search")
    @test all(startswith("Derived"), out[4:(end - 2)])
    @test out[end - 1] == "Target found"
    @test startswith(out[end], "Index")
end

@testset "time limit works" begin
    data = "p(x*,y*,z*,r_x,r_y,r_z)\np(y)"
    graph = """
      x -> z
      z -> y
      y -> r_y
      x <-> y
      r_y -> r_x
      r_y -> r_z
      r_y <-> r_x
      r_y <-> r_z
      r_z <-> r_x
    """
    md = "r_x : x, r_y : y, r_z : z"
    for benchmark_rules in (true, false)
        out = dosearch(data, "p(y|do(x))", graph; missing_data = md, heuristic = true,
                       improve = false, time_limit = 0.001, benchmark_rules)
        @test out.time < 1.0
        @test length(out.rule_times) == (benchmark_rules ? length(out.rules) : 0)
        @test out.limit_reached && !out.identifiable
    end
    # A zero limit stops at once, but a query among the inputs is still found.
    out = dosearch("p(x,y,z)", "p(y|do(x))", "x -> z\nz -> y\nx <-> y"; time_limit = 0)
    @test out.limit_reached && !out.identifiable
    @test sprint(show, MIME"text/plain"(), out) ==
          "The query p(y|do(x)) was not identified before the search reached its limit."
    out = dosearch("p(y|do(x))", "p(y|do(x))", "x -> y"; time_limit = 0)
    @test out.identifiable && !out.limit_reached
    # trivially non-identifiable queries are decided without searching
    out = dosearch("p(x)", "p(y)", "x -> y"; time_limit = 0)
    @test !out.identifiable && !out.limit_reached
end

@testset "term limit works" begin
    frontdoor = ("p(x,y,z)", "p(y|do(x))", "x -> z\nz -> y\nx <-> y")
    out = dosearch(frontdoor...; term_limit = 10)
    @test out.limit_reached && !out.identifiable
    out = dosearch(frontdoor...; term_limit = 1000)
    @test out.identifiable && !out.limit_reached
    # The limit is checked before each term is expanded.
    @test dosearch(frontdoor...; term_limit = out.terms + 1).identifiable
    @test dosearch(frontdoor...; term_limit = out.terms ÷ 2).limit_reached
    out = dosearch("p(x,y)", "p(y|do(x))", "x -> y\nx <-> y"; term_limit = 1000)
    @test !out.identifiable && !out.limit_reached
    # LDAGs
    out = dosearch("p(A,Y)", "p(Y|X,A=1)", "X -> Y : A = 1\nA -> Y"; term_limit = 2)
    @test out.limit_reached && !out.identifiable
end

@testset "warnings" begin
    @test_logs (:warn, r"There are response indicators that are not present in any input distribution: r_x") dosearch(
        "p(x*)", "p(x)", "x -> y"; missing_data = "r_x : x")
    @test_logs (:warn, "Both lower case and upper case inputs detected.") dosearch("p(x,X)", "p(X)", "x -> X")
    @test_logs dosearch("p(x*)", "p(x)", "x -> y"; missing_data = "r_x : x", warn = false)
    @test_logs dosearch("p(x,X)", "p(X)", "x -> X"; warn = false)
end

@testset "extra variables do not influence identifiability" begin
    graph = "x -> z \n z -> y \n x <-> y"
    @test dosearch("p(x,y,z,w)", "p(y|do(x))", graph; warn = false).identifiable ==
          dosearch("p(x,y,z)", "p(y|do(x))", graph).identifiable
end

@testset "rule 1 works (redundant rule)" begin
    graph = "x -> z \n y -> z"
    out = dosearch("p(y|x)", "p(y)", graph; rules = [-1, 1], draw_derivation = true)
    @test out.formula == "p(y|x)"
    out = dosearch("p(y)", "p(y|x)", graph; rules = [-1, 1], draw_derivation = true)
    @test out.formula == "p(y)"
    # R passes `rules = c(seq(-3, 3), ...)` here, which includes the non-rule
    # code 0. Dosearch.jl rejects code 0 (see the fixture test).
    rules = [-3, -2, -1, 1, 2, 3, 4, 5, -6, 6, -7, 7, -8, 8, 9, 10]
    out = dosearch("p(x*,y*,r_x,r_y)\np(y)", "p(y|do(x))", graph;
                   missing_data = "r_x : x, r_y : y", rules, heuristic = true)
    @test out.identifiable
    @test_throws ArgumentError dosearch("p(y)", "p(y|x)", graph; rules = [0, 1])
    @test_throws ArgumentError dosearch("p(y)", "p(y|x)", graph; rules = [11])
end

@testset "division rules work" begin
    graph = "x -> y\ny -> r_x\nr_x -> r_y"
    md = "r_x : x, r_y : y"
    check(data, query, rules) = dosearch(data, query, graph; missing_data = md,
                                         draw_derivation = true, rules, warn = false)
    @test check("p(r_x=1,r_y=1) \n p(r_x=1)", "p(r_y=1|r_x=1)", 7).formula == "p(r_y = 1|r_x = 1)"
    @test check("p(r_x=1,r_y=1) \n p(r_x=1|r_y=1)", "p(r_y=1)", -7).formula == "p(r_y = 1)"
    @test check("p(r_x=1) \n p(r_x=1,r_y=1)", "p(r_y = 1|r_x = 1)", 8).formula == "p(r_y = 1|r_x = 1)"
    @test check("p(r_x=1|r_y=1) \n p(r_x=1,r_y=1)", "p(r_y = 1)", -8).formula == "p(r_y = 1)"
    @test check("p(r_x=1|z) \n p(r_x=1,r_y=1)", "p(r_y = 1|r_x = 1)", [-1, 8]).formula ==
          raw"\frac{p(r_x = 1,r_y = 1)}{p(r_x = 1|z)}"
    @test check("p(r_x=1|z) \n p(r_x=1,r_y=1)", "p(r_y = 1|r_x = 1)",
                [Dosearch.DeleteObservation, Dosearch.DenominatorMarginal]).identifiable
end

@testset "alternative distribution format works" begin
    data = "p(x)\np(y|x)\np(z|x,y)"
    expected = dosearch(data, "p(y|do(x))", BACKDOOR).formula
    data_alt = [(x = 0,), (y = 0, x = 2), (z = 0, x = 2, y = 2)]
    @test dosearch(data_alt, (y = 0, x = 1), BACKDOOR).formula == expected
    @test dosearch(["p(x)", "p(y|x)", "p(z|x,y)"], (y = 0, x = 1), BACKDOOR).formula == expected
    @test Dosearch.roles_to_distribution((y = 0, x = 1, w = 2)) == "p(y|do(x),w)"
end

@testset "results can be shown and reused" begin
    out = dosearch("p(x,y)", "p(y|do(x))", "x -> y")
    @test sprint(show, MIME"text/plain"(), out) == "p(y|x)"
    @test sprint(show, out) == "DosearchResult(identifiable = true)"
    out = dosearch("p(x,y)", "p(y|do(x))", "x -> y\nx <-> y")
    @test sprint(show, MIME"text/plain"(), out) ==
          "The query p(y|do(x)) is not identifiable by do-search."
    out = dosearch("p(x,y)", "p(y|do(x))", "x -> y"; formula = false)
    @test out.formula === nothing
    @test sprint(show, MIME"text/plain"(), out) == "The query p(y|do(x)) is identifiable."
    problem = DosearchProblem("p(x,y,z)", "p(y|do(x))", BACKDOOR)
    @test problem.labels == ["y", "x", "z"]
    @test dosearch(problem).formula == dosearch(problem; strategy = FoundOrder()).formula
    @test dosearch(problem; strategy = ScoreHeuristic()).formula ==
          dosearch(problem; heuristic = true).formula
    @test occursin("query: p(y|do(x))", sprint(show, MIME"text/plain"(), problem))
end

@testset "derivations" begin
    out = dosearch("p(x,y)", "p(y|do(x))", "x -> y"; draw_derivation = true)
    @test out.dot == """
        strict digraph InferenceTree {
        n1[shape=polygon,sides=4,label="p(y|do(x))"]
        n2[shape=polygon,sides=4,label="p(y|x)"]
        n2 -> n1[label="R2"]
        n3[shape=polygon,sides=4,label="p(y,x)"]
        n3 -> n2[label="C"]
        }
        """
    out = dosearch("p(x)", "p(y)", "x -> y"; draw_derivation = true)
    @test out.dot == ""   # trivially non-identifiable, as in R
    out = dosearch("p(x)", "p(y)", "x -> y"; draw_derivation = true, improve = false)
    @test out.dot == "strict digraph InferenceTree {\n}\n"
    out = dosearch("p(x,y)", "p(y|do(x))", "x -> y")
    @test out.dot === nothing
end
