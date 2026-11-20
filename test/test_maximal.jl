# Tests of the candidates and the search for maximal identifiable graphs.

using Random, TOML

isdefined(@__MODULE__, :load_fixtures) || include("fixture_utils.jl")

const CHAIN = "x -> z\nz -> y"

@testset "candidates" begin
    @test BidirectedPair("y", "x") == BidirectedPair("x", "y")
    @test hash(BidirectedPair("y", "x")) == hash(BidirectedPair("x", "y"))
    @test sprint(show, BidirectedPair("y", "x")) == "BidirectedPair(\"x\", \"y\")"
    @test sprint(show, MIME"text/plain"(), BidirectedPair("y", "x")) == "x <-> y"
    @test_throws ArgumentError BidirectedPair("x", "x")

    p = DosearchProblem("p(x,y,z)", "p(y|do(x))", CHAIN)
    @test confounding_candidates(p) ==
          [BidirectedPair("y", "z"), BidirectedPair("x", "z"), BidirectedPair("x", "y")]
    q = with_candidates(p, [BidirectedPair("x", "y")])
    @test q.bidirected == [(3, 2)]
    @test dosearch(q).formula == dosearch("p(x,y,z)", "p(y|do(x))", CHAIN * "\nx <-> y").formula
    @test confounding_candidates(q) == [BidirectedPair("y", "z"), BidirectedPair("x", "z")]
    @test_throws ArgumentError with_candidates(p, [BidirectedPair("x", "w")])

    # transportability nodes and proxies are left out, and rejected
    p = DosearchProblem("p(x,z,y|s)\np(y,z|t,do(x))", "p(y|do(x))",
                        "x -> z\nz -> y\nx -> s\nt -> z"; transportability = "t", selection_bias = "s")
    @test all(c -> c.x != "t" && c.y != "t", confounding_candidates(p))
    @test BidirectedPair("s", "y") in confounding_candidates(p)
    @test_throws ArgumentError with_candidates(p, [BidirectedPair("t", "y")])
    p = DosearchProblem("p(x*,y*,r_x,r_y)\np(y)", "p(y|do(x))", "x -> y\ny -> r_y\nr_y -> r_x";
                        missing_data = "r_x : x, r_y : y")
    @test all(c -> !endswith(c.x, "*") && !endswith(c.y, "*"), confounding_candidates(p))
    @test length(confounding_candidates(p)) == 6
    @test_throws ArgumentError with_candidates(p, [BidirectedPair("x*", "y")])
end

@testset "maximal identifiable graphs: examples" begin
    @test maximal_identifiable_graphs("p(x,y,z)", "p(y|do(x))", CHAIN) ==
          [[BidirectedPair("y", "z")], [BidirectedPair("x", "y")]]
    @test maximal_identifiable_graph("p(x,y,z)", "p(y|do(x))", CHAIN) == [BidirectedPair("y", "z")]
    @test maximal_identifiable_graph("p(x,y,z)", "p(y|do(x))", CHAIN;
                                     candidates = [BidirectedPair("x", "y"), BidirectedPair("y", "z")]) ==
          [BidirectedPair("x", "y")]
    # back-door: either confounder of z can be added, but not both
    @test Set(maximal_identifiable_graphs("p(x,y,z)", "p(y|do(x))", "x -> y\nz -> x\nz -> y")) ==
          Set([[BidirectedPair("y", "z")], [BidirectedPair("x", "z")]])
    # not identifiable to begin with
    @test isempty(maximal_identifiable_graphs("p(x,y)", "p(y|do(x))", "x -> y\nx <-> y"))
    @test maximal_identifiable_graph("p(x,y)", "p(y|do(x))", "x -> y\nx <-> y") === nothing
    # nothing can be added
    @test maximal_identifiable_graphs("p(x,y)", "p(y|do(x))", "x -> y") == [BidirectedPair[]]
    @test maximal_identifiable_graph("p(x,y)", "p(y|do(x))", "x -> y") == BidirectedPair[]
    # everything can be added
    @test maximal_identifiable_graphs("p(x,y,z)", "p(y)", CHAIN) == [confounding_candidates(
        DosearchProblem("p(x,y,z)", "p(y)", CHAIN))]
end

# Maximal identifiable and minimal non-identifiable sets of candidates, by
# checking every subset.
function brute_force_border(problem, candidates; kwargs...)
    m = length(candidates)
    has(s, i) = (s >> (i - 1)) & 1 == 1
    members(s) = Set(candidates[[i for i in 1:m if has(s, i)]])
    ident = Dict(s => dosearch(with_candidates(problem, collect(members(s))); formula = false,
                               kwargs...).identifiable for s in 0:(2^m - 1))
    maximal = [s for s in 0:(2^m - 1) if ident[s] && all(i -> has(s, i) || !ident[s | 1 << (i - 1)], 1:m)]
    minimal = [s for s in 0:(2^m - 1) if !ident[s] && all(i -> !has(s, i) || ident[s & ~(1 << (i - 1))], 1:m)]
    return Set(members.(maximal)), Set(members.(minimal))
end
brute_force_maximal(problem, candidates; kwargs...) = first(brute_force_border(problem, candidates; kwargs...))

@testset "maximal identifiable graphs agree with brute force" begin
    cases = load_fixtures()
    rng = Xoshiro(11)
    nchecked = 0
    for case in cases
        haskey(case, "rules") && continue
        problem = DosearchProblem(case["data"], case["query"], case["graph"];
                                  transportability = get(case, "transportability", nothing),
                                  selection_bias = get(case, "selection_bias", nothing),
                                  missing_data = get(case, "missing_data", nothing),
                                  md_sym = get(case, "md_sym", "1"), warn = false)
        candidates = confounding_candidates(problem)
        length(candidates) > 7 && (candidates = shuffle(rng, candidates)[1:7])
        opts = (heuristic = case["heuristic"], improve = case["improve"])
        expected, expected_minimal = brute_force_border(problem, candidates; opts...)
        found = maximal_identifiable_graphs(problem; candidates, opts...)
        @test allunique(found)
        @test Set(Set.(found)) == expected
        @test Set(Set.(found.minimal_nonidentifiable)) == expected_minimal
        one = maximal_identifiable_graph(problem; candidates, opts...)
        @test (one === nothing) == isempty(expected)
        one === nothing || @test Set(one) in expected
        nchecked += 1
        nchecked >= 80 && break
    end
    @test nchecked == 80
end

@testset "added confounding never makes a query identifiable" begin
    # The search relies on this monotonicity; check it on chains of added
    # confounders for the fixture problems.
    rng = Xoshiro(5)
    for case in load_fixtures()[1:3:end]
        rules = get(case, "rules", nothing)
        rules !== nothing && 0 in rules && continue
        problem = DosearchProblem(case["data"], case["query"], case["graph"];
                                  transportability = get(case, "transportability", nothing),
                                  selection_bias = get(case, "selection_bias", nothing),
                                  missing_data = get(case, "missing_data", nothing),
                                  md_sym = get(case, "md_sym", "1"), warn = false)
        order = shuffle(rng, confounding_candidates(problem))
        ids = [dosearch(with_candidates(problem, order[1:k]); formula = false, rules,
                        heuristic = case["heuristic"], improve = case["improve"]).identifiable
               for k in 0:length(order)]
        @test issorted(ids; rev = true)
    end
end

@testset "dualize and advance" begin
    # Berge's algorithm for minimal transversals
    T = Dosearch.add_edge_transversals([BitSet()], BitSet([1, 2]))
    T = Dosearch.add_edge_transversals(T, BitSet([2, 3]))
    @test Set(T) == Set([BitSet([2]), BitSet([1, 3])])
    @test isempty(Dosearch.add_edge_transversals(T, BitSet()))

    # Every maximal and every minimal non-identifiable set must be searched;
    # check that the search count stays within a small factor of that.
    problem = DosearchProblem("p(a,b,c,d,x,y)", "p(y|do(x))",
                              "x -> a\na -> b\nb -> y\nc -> x\nc -> y\nd -> c")
    candidates = confounding_candidates(problem)
    o = Dosearch.MonotoneOracle(problem, candidates, Dosearch.search_options())
    Dosearch.isidentifiable!(o, BitSet())
    sets, minimal_nonid = Dosearch.dualize_and_advance!(o)
    searches = o.searches
    @test length(sets) == 16
    @test allunique(sets)
    universe = BitSet(eachindex(candidates))
    @test Set(minimal_nonid) == Set(foldl(Dosearch.add_edge_transversals,
                                          [setdiff(universe, s) for s in sets]; init = [BitSet()]))
    @test length(minimal_nonid) == 65
    @test all(s -> Dosearch.isidentifiable!(o, s), sets)
    @test !any(t -> Dosearch.isidentifiable!(o, t), minimal_nonid)
    @test o.searches == searches   # all of them were searched already
    @test searches <= 3 * (length(sets) + length(minimal_nonid))
end

@testset "derivation certificates" begin
    # If every d-separation used by a derivation holds in a graph with more
    # confounding, the query is identifiable there too.
    rng = Xoshiro(9)
    nused = 0
    for case in load_fixtures()[1:2:end]
        rules = get(case, "rules", nothing)
        rules !== nothing && 0 in rules && continue
        problem = DosearchProblem(case["data"], case["query"], case["graph"];
                                  transportability = get(case, "transportability", nothing),
                                  selection_bias = get(case, "selection_bias", nothing),
                                  missing_data = get(case, "missing_data", nothing),
                                  md_sym = get(case, "md_sym", "1"), warn = false)
        opts = (; rules, heuristic = case["heuristic"], improve = case["improve"])
        st, _ = Dosearch.solve(problem; opts...)
        st.target === nothing && continue
        seps = Dosearch.derivation_separations(st)
        candidates = confounding_candidates(problem)
        for _ in 1:5
            extended = with_candidates(problem, filter(_ -> rand(rng, Bool), candidates))
            g = Dosearch.causal_graph(extended)
            if all(s -> Dosearch.dseparated(g, s), seps)
                nused += 1
                @test dosearch(extended; formula = false, opts...).identifiable
            end
        end
    end
    @test nused > 100
end

@testset "maximal graph results" begin
    r = maximal_identifiable_graphs("p(x,y,z)", "p(y|do(x))", CHAIN; formulas = true)
    @test r isa MaximalGraphs{BidirectedPair}
    @test r.maximal == collect(r)
    @test r.minimal_nonidentifiable == [[BidirectedPair("x", "z")],
                                        [BidirectedPair("y", "z"), BidirectedPair("x", "y")]]
    @test r.formulas == [dosearch(with_candidates(r.problem, s)).formula for s in r]
    @test r.formulas[2] == raw"\sum_{z}\left(p(z|x)\sum_{x}\left(p(x)p(y|z,x)\right)\right)"
    @test r.searches > 0
    @test sprint(show, MIME"text/plain"(), r) == """
        2 maximal identifiable sets of 3 candidates:
          y <-> z
          x <-> y
        2 minimal non-identifiable sets:
          x <-> z
          y <-> z; x <-> y"""
    @test isempty(maximal_identifiable_graphs("p(x,y,z)", "p(y|do(x))", CHAIN).formulas)
    # not identifiable to begin with: the empty set is the minimal blocking set
    r = maximal_identifiable_graphs("p(x,y,z)", "p(y|do(x))", CHAIN * "\nx <-> z")
    @test isempty(r)
    @test r.minimal_nonidentifiable == [BidirectedPair[]]
    # everything can be added: nothing blocks
    r = maximal_identifiable_graphs("p(x,y,z)", "p(y)", CHAIN)
    @test length(r) == 1 && isempty(r.minimal_nonidentifiable)
end

@testset "latent causes" begin
    c = LatentCause("z", "x", "y", "x")
    @test c.variables == ["x", "y", "z"]
    @test c == LatentCause(["y", "z", "x"])
    @test hash(c) == hash(LatentCause("x", "y", "z"))
    @test sprint(show, c) == "LatentCause(\"x\", \"y\", \"z\")"
    @test sprint(show, MIME"text/plain"(), c) == "latent cause of x, y, z"
    @test_throws ArgumentError LatentCause(["x", "x"])
    p = DosearchProblem("p(x,y,z,w)", "p(y|do(x))", "x -> z\nz -> y\nw -> x")
    pairs = [BidirectedPair("x", "y"), BidirectedPair("x", "z"), BidirectedPair("y", "z")]
    edges(q) = Set(minmax(e...) for e in q.bidirected)
    @test edges(with_candidates(p, [c])) == edges(with_candidates(p, pairs))
    @test with_candidates(p, [LatentCause("x", "y")]).bidirected ==
          with_candidates(p, [BidirectedPair("x", "y")]).bidirected
    # a search over latent causes, and over mixed kinds of candidates
    for candidates in ([LatentCause("x", "y", "w"), LatentCause("y", "z"), LatentCause("x", "w"),
                        LatentCause("w", "z")],
                       Candidate[LatentCause("x", "y", "w"), BidirectedPair("y", "z"), BidirectedPair("w", "z")])
        expected, expected_minimal = brute_force_border(p, candidates)
        r = maximal_identifiable_graphs(p; candidates)
        @test Set(Set.(r)) == expected
        @test Set(Set.(r.minimal_nonidentifiable)) == expected_minimal
    end
end

@testset "search limits" begin
    p = DagProblem("p(x,y,z,w)", "p(y|do(x))", "x -> z\nz -> y\nw -> x\nw -> z")
    full = maximal_identifiable_graphs(p)
    @test isempty(full.undecided)
    # At this term limit some searches stop before they decide: their sets
    # count as non-identifiable and are reported as undecided.
    r = maximal_identifiable_graphs(p; term_limit = 166)
    @test !isempty(r) && !isempty(r.minimal_nonidentifiable) && !isempty(r.undecided)
    for s in r   # identifiable, though perhaps not maximal
        @test dosearch(with_candidates(p, s)).identifiable
        @test any(t -> issubset(s, t), full)
    end
    # proven non-identifiable sets are minimal non-identifiable sets
    @test issubset(Set(Set.(r.minimal_nonidentifiable)), Set(Set.(full.minimal_nonidentifiable)))
    for s in r.undecided
        @test dosearch(with_candidates(p, s); term_limit = 166).limit_reached
    end
    @test occursin("1 minimal undecided set (search limit reached):", sprint(show, MIME"text/plain"(), r))
    # The same result in parallel: term limits are deterministic.
    o = Dosearch.MonotoneOracle(p, confounding_candidates(p), Dosearch.search_options(; term_limit = 166);
                                parallel = true)
    @test Dosearch.isidentifiable!(o, BitSet())
    maximal, minimal = Dosearch.dualize_and_advance!(o)
    @test Set(Set.(r)) == Set(Set(o.candidates[collect(s)]) for s in maximal)
    # formulas are computed without the limit
    r = maximal_identifiable_graphs(p; term_limit = 166, formulas = true)
    @test all(!isempty, r.formulas)
    @test maximal_identifiable_graph(p; term_limit = 166) in r
    # a limit large enough changes nothing
    r = maximal_identifiable_graphs(p; term_limit = 1000)
    @test r == full && r.minimal_nonidentifiable == full.minimal_nonidentifiable && isempty(r.undecided)
    # The original model is undecided.
    r = maximal_identifiable_graphs(p; time_limit = 0)
    @test isempty(r) && isempty(r.minimal_nonidentifiable) && r.undecided == [BidirectedPair[]]
    @test maximal_identifiable_graph(p; time_limit = 0) === nothing
end

# LDAGs ------------------------------------------------------------------------

const LDAG_CASES = TOML.parsefile(joinpath(@__DIR__, "fixtures", "ldag_fixtures.toml"))["case"]

# A labeled latent cause of a pair, on the edge into one of them, assigning
# another parent of it, or `nothing` if it has no other parent.
function random_labeled_cause(rng, p::LdagProblem, c::BidirectedPair)
    v = rand(rng, [c.x, c.y])
    y = findfirst(==(v), p.labels)
    pa = [p.labels[i] for (i, j) in p.edges if j == y && !(i in p.intervention_vars)]
    isempty(pa) && return nothing
    return LatentCause(c.x, c.y; labels = [v => rand(rng, pa) * " = " * string(rand(rng, 0:1))])
end

@testset "LDAG candidates" begin
    graph = "I_X -> X\nX -> Y : A = 1\nA -> Y\nA -> X : I_X = 1"
    p = LdagProblem("p(A,X,Y)", "p(Y|X,I_X=1)", graph)
    @test confounding_candidates(p) == [BidirectedPair("X", "Y"), BidirectedPair("A", "X"), BidirectedPair("A", "Y")]
    c = LatentCause("X", "Y"; labels = ["X" => "A=0"])
    @test c == LatentCause("Y", "X"; labels = ["X" => " A = 0 "])
    @test c != LatentCause("X", "Y")
    @test hash(c) == hash(LatentCause("Y", "X"; labels = ["X" => "A = 0"]))
    @test sprint(show, c) == "LatentCause(\"X\", \"Y\"; labels = [\"X\" => \"A = 0\"])"
    @test sprint(show, MIME"text/plain"(), c) == "latent cause of X, Y (into X : A = 0)"
    q = with_candidates(p, [c, BidirectedPair("A", "Y")])
    @test q.labels[1:4] == p.labels   # the numbering is kept
    @test q.graph_text == graph * """

        U_X_Y_X_A0 -> X : A = 0; I_X = 1
        U_X_Y_X_A0 -> Y
        U_A_Y -> A
        U_A_Y -> Y"""
    @test with_candidates(p, BidirectedPair[]) === p
    # The same structure twice gets a primed name.
    q = with_candidates(p, [BidirectedPair("X", "Y"), LatentCause("X", "Y")])
    @test q.labels[5:6] == ["U_X_Y", "U_X_Y'"]
    @test_throws ArgumentError with_candidates(p, [BidirectedPair("I_X", "Y")])
    @test_throws ArgumentError with_candidates(p, [BidirectedPair("X", "Z")])
    @test_throws ArgumentError with_candidates(p, [LatentCause("X", "Y"; labels = ["Y" => "I_X = 1"])])
    @test_throws ArgumentError LatentCause("X", "Y"; labels = ["Z" => "A = 1"])
    @test_throws ArgumentError LatentCause("X", "Y"; labels = ["Y" => "A = 1", "Y" => "A = 0"])
    @test_throws ArgumentError with_candidates(DagProblem("p(x,y)", "p(y)", "x -> y"),
                                               [LatentCause("x", "y"; labels = ["y" => "x = 1"])])
    # Confounding that vanishes where A = 0 keeps the effect identifiable,
    # since X has no effect on Y where A = 1.
    r = maximal_identifiable_graphs(p; candidates = [c, BidirectedPair("X", "Y"), BidirectedPair("A", "Y"),
                                                     BidirectedPair("A", "X")], formulas = true)
    @test Set(Set.(r)) == Set([Set([c, BidirectedPair("A", "Y")]), Set([c, BidirectedPair("A", "X")])])
    @test Set(Set.(r.minimal_nonidentifiable)) ==
          Set([Set([BidirectedPair("X", "Y")]), Set([BidirectedPair("A", "Y"), BidirectedPair("A", "X")])])
    @test all(==("\\sum_{A}\\left[\\left(p(Y|X,A = 0)\\left[p(A)\\right]_{A = 0}\\right) /\\ " *
                 "\\left(p(Y|A = 1)\\left[p(A)\\right]_{A = 1}\\right)\\right]"), r.formulas)
    @test r.problem === p   # A is a context variable already
    @test maximal_identifiable_graphs("p(A,X,Y)", "p(Y|X,I_X=1)", graph) == maximal_identifiable_graphs(p)
    # Variables in candidate labels become context variables throughout.
    p2 = LdagProblem("p(x,y,z)", "p(y|x)", "z -> x\nz -> y\nx -> y : z = 1")
    c2 = LatentCause("x", "w"; labels = ["x" => "z = 0"])
    @test_throws ArgumentError with_candidates(p2, [c2])
    p3 = LdagProblem("p(x,y,z,w)", "p(y|x)", "z -> x\nw -> x\nz -> y\nx -> y")
    c3 = LatentCause("x", "y"; labels = ["x" => "w = 1"])
    r = maximal_identifiable_graphs(p3; candidates = [c3])
    @test r.problem.labels[collect(r.problem.context_vars)] == ["w"]
    @test isempty(p3.context_vars)
end

@testset "LDAG maximal graphs agree with brute force" begin
    rng = Xoshiro(5)
    nchecked = 0
    certified = 0
    for case in LDAG_CASES
        p = LdagProblem(case["data"], case["query"], case["graph"];
                        context_vars = get(case, "con_vars", String[]), warn = false)
        pairs = confounding_candidates(p)
        candidates = Candidate[shuffle(rng, pairs)[1:min(end, 3)]...]
        for c in shuffle(rng, pairs)[1:min(end, 2)]
            l = random_labeled_cause(rng, p, c)
            l === nothing || push!(candidates, l)
        end
        r = maximal_identifiable_graphs(p; candidates)
        expected, expected_minimal = brute_force_border(r.problem, candidates)
        @test Set(Set.(r)) == expected
        @test Set(Set.(r.minimal_nonidentifiable)) == expected_minimal
        o = Dosearch.MonotoneOracle(r.problem, candidates, Dosearch.search_options())
        if Dosearch.isidentifiable!(o, BitSet())
            Dosearch.dualize_and_advance!(o)
            certified += o.certified
        end
        nchecked += 1
        nchecked >= 20 && break
    end
    @test nchecked == 20
    @test certified > 0   # derivations are reused across sets of candidates
end
