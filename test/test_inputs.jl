# Graphs and distributions as values: they render to the text syntax and give
# the same problems as the text.

@testset "CausalGraph" begin
    g = CausalGraph("x" => "y", :z => :x, "z" => "y", "x" => "y")
    @test length(g.edges) == 3                         # repeated edges kept once
    @test string(g) == "x -> y\nz -> x\nz -> y\n"
    @test sprint(show, MIME"text/plain"(), g) == "CausalGraph with 3 edges:\n  x -> y\n  z -> x\n  z -> y"
    @test eval(Meta.parse(repr(g))) == g
    h = g ∪ CausalGraph(BidirectedPair("y", "x"))
    @test string(h) == "x -> y\nz -> x\nz -> y\nx <-> y\n"
    @test CausalGraph(g, BidirectedPair("x", "y")) == h
    @test CausalGraph("A$k" => "A$(k + 1)" for k in 1:3) == CausalGraph("A1" => "A2", "A2" => "A3", "A3" => "A4")
    @test CausalGraph(["x" => "y"], ("z" => "y",)) == CausalGraph("x" => "y", "z" => "y")
    # from text, and back
    text = """
        x -> y
        z <-> x
        z -> y
    """
    @test CausalGraph(text) == CausalGraph("x" => "y", BidirectedPair("x", "z"), "z" => "y")
    @test CausalGraph(text) != CausalGraph("x" => "y", "z" => "y", BidirectedPair("x", "z"))   # order matters
    @test CausalGraph(string(CausalGraph(text))) == CausalGraph(text)
    ldag = """
        X -> Y : A = 1; B = 0, C = 1
        A -> Y
        B -> Y
        C -> Y
    """
    lg = CausalGraph(ldag)
    @test lg.edges[1:2] == [LabeledEdge("X" => "Y", "A" => 1), LabeledEdge("X" => "Y", "B" => 0, "C" => 1)]
    @test string(lg) == "X -> Y : A = 1\nX -> Y : B = 0, C = 1\nA -> Y\nB -> Y\nC -> Y\n"
    @test sprint(show, MIME"text/plain"(), lg.edges[2]) == "X -> Y : B = 0, C = 1"
    @test eval(Meta.parse(repr(lg))) == lg
    # errors
    @test_throws "self-loop" CausalGraph("x" => "x")
    @test_throws "cycles" CausalGraph("x" => "y", "y" => "x")
    @test_throws "Invalid variable name" CausalGraph("x y" => "z")
    @test_throws "Invalid variable name" CausalGraph("x" => "y|z")
    @test_throws "string or a symbol" CausalGraph(1 => 2)
    @test_throws "no bidirected edges" CausalGraph("a" => "b", BidirectedPair("a", "c"), LabeledEdge("c" => "b", "a" => 1))
    @test_throws "0 or 1" LabeledEdge("x" => "y", "a" => 2)
    @test_throws "at least one assignment" LabeledEdge("x" => "y")
    @test_throws "twice" LabeledEdge("x" => "y", "a" => 1, "a" => 0)
    @test_throws "values 0 or 1" CausalGraph("x -> y : a = 2")
    @test_throws "not a line break" CausalGraph("x -> y\\nz -> y")
end

@testset "P" begin
    d = P("y"; given = "z", intervened = "x")
    @test string(d) == "p(y | do(x), z)"
    @test sprint(show, MIME"text/plain"(), d) == "p(y | do(x), z)"
    @test repr(d) == "P(\"y\"; given = \"z\", intervened = \"x\")"
    @test eval(Meta.parse(repr(d))) == d
    @test string(P(:y, :w; given = [:z, :v])) == "p(y, w | z, v)"
    @test string(P(["y", "w"])) == "p(y, w)"
    e = P("Y", "Z" => 0; given = ["X", "I_X" => 1])
    @test string(e) == "p(Y, Z = 0 | X, I_X = 1)"
    @test eval(Meta.parse(repr(e))) == e
    @test string(P("x*"; given = "r_x" => 1)) == "p(x* | r_x = 1)"
    @test_throws "at least one variable" P()
    @test_throws "twice" P("y"; given = "y")
    @test_throws "Invalid variable name" P("y,z")
    @test_throws "0, 1 or a string" P("y"; given = "x" => 0.5)
end

@testset "problems from values" begin
    # the same problems as from text, with the same numbering of the variables
    text = """
        L₀ -> A₀
        L₀ -> Y₁
        A₀ -> Y₁
        L₀ -> L₁
        A₀ -> L₁
        Y₁ -> L₁
        A₀ -> A₁
        L₁ -> A₁
        Y₁ -> A₁
        L₀ -> Y₂
        L₁ -> Y₂
        A₁ -> Y₂
        Y₁ -> Y₂
        L₀ <-> A₀
    """
    a = DosearchProblem("p(L₀, A₀, Y₁, L₁, A₁, Y₂)", "p(Y₂ | do(A₀, A₁))", text)
    b = DosearchProblem(P("L₀", "A₀", "Y₁", "L₁", "A₁", "Y₂"), P("Y₂"; intervened = ["A₀", "A₁"]), CausalGraph(text))
    @test (a.labels, a.directed, a.bidirected, a.data, a.query) == (b.labels, b.directed, b.bidirected, b.data, b.query)
    @test dosearch(a).formula == dosearch(b).formula
    @test DagProblem([P("x", "y")], P("y"; intervened = "x"), CausalGraph("x" => "y")) isa DagProblem
    # data as a mixed vector of values and text
    r = dosearch([P("z"; intervened = "x"), "p(y|do(z))"], P("y"; intervened = "x"),
                 CausalGraph("x" => "z", "z" => "y", BidirectedPair("x", "y")))
    @test r.formula == "\\sum_{z}\\left(p(z|do(x))p(y|do(z))\\right)"
    # special nodes as collections
    data = "p(x*,y*,r_x,r_y)\np(y)"
    graph = "x -> y\ny -> r_y\nr_y -> r_x"
    @test dosearch(data, "p(y|do(x))", CausalGraph(graph); missing_data = ["r_x" => "x", "r_y" => "y"]).formula ==
          dosearch(data, "p(y|do(x))", graph; missing_data = "r_x : x, r_y : y").formula
    sb = dosearch("p(x,y,z|s)\np(x)", "p(y|do(x))", "x -> y\nz -> y\nz -> s"; selection_bias = [:s], warn = false)
    @test sb.formula == dosearch("p(x,y,z|s)\np(x)", "p(y|do(x))", "x -> y\nz -> y\nz -> s";
                                 selection_bias = "s", warn = false).formula
    # LDAGs
    lg = CausalGraph("I_X" => "X", LabeledEdge("L" => "X", "I_X" => 1), "L" => "Y", "X" => "Y")
    @test DosearchProblem(P("L", "X", "Y"), P("Y"; given = ["X", "I_X" => 1]), lg) isa LdagProblem
    @test dosearch(P("L", "X", "Y"), P("Y"; given = ["X", "I_X" => 1]), lg).formula ==
          "\\sum_{L}\\left(p(L)p(Y|X,L)\\right)"
    # maximal identifiable graphs: x <-> y (front door) and z <-> y
    @test Set(Set.(maximal_identifiable_graphs(P("x", "y", "z"), P("y"; intervened = "x"),
                                               CausalGraph("x" => "z", "z" => "y")))) ==
          Set([Set([BidirectedPair("x", "y")]), Set([BidirectedPair("y", "z")])])
end
