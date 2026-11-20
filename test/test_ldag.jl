# Tests of the search for labeled DAGs (LDAGs): R fixtures
# (fixtures/ldag_fixtures.toml, see fixtures/generate_ldag_fixtures.R) with the
# examples and tests of the R package (test-csisearch.R) and random LDAGs, and
# the LDAG cases of test-errors.R.

using TOML
isdefined(@__MODULE__, :check_derivation) || include("derivation_utils.jl")

const LDAG_FIXTURES = TOML.parsefile(joinpath(@__DIR__, "fixtures", "ldag_fixtures.toml"))["case"]

ldag_problem(case) = LdagProblem(case["data"], case["query"], case["graph"];
                                 context_vars = get(case, "con_vars", String[]), warn = false)

# The contexts and local CSIs that R passes to C++, without the edges that
# vanish in every context. Intentional difference: R means to leave those
# edges out, but it compares a named with an unnamed vector (`all.equal`) and
# keeps them; the C++ code then adds them back to the graph when it leaves a
# context that removes them. Dosearch.jl leaves them out. The search results
# agree with R on every fixture all the same.
function r_contexts(case)
    edges = Set(Tuple(e) for e in case["edges"])
    vars = case["vars"]
    intv = VarSet(Int[findfirst(==(v), vars) for v in case["intv_vars"]])
    label_map = []
    for (setvars, contexts) in case["label_map"]
        classes = Vector{Vector{Tuple{Int,Int}}}()
        out = []
        for (zero, one, _, from, to) in contexts
            removed = [(f, t) for (f, t) in zip(from, to) if (f, t) in edges]
            class = findfirst(==(removed), classes)
            class === nothing && (push!(classes, removed); class = length(classes))
            push!(out, (zero, one, class, removed))
        end
        is_intervention = VarSet(UInt64(setvars)) ⊆ intv
        (is_intervention || any(c -> !isempty(c[4]), out)) && push!(label_map, (setvars, out))
    end
    local_csi = [Tuple(l) for l in case["local_csi"]
                 if (trailing_zeros(l[1]) + 1, trailing_zeros(l[2]) + 1) in edges]
    return label_map, local_csi
end

# Intentional difference: a variable fixed to a value in the query or an input
# is a context variable, as if it were in the R package's `con_vars`. So the
# case without `con_vars` gives R's result for the same case with them.
const SAME_RESULT_AS = Dict("custom context variables absent" => "custom context variables")

@testset "LDAG fixtures" begin
    @test length(LDAG_FIXTURES) > 150
    for case in LDAG_FIXTURES
        @testset "$(case["name"])" begin
            p = ldag_problem(case)
            # the problem in internal form
            @test p.labels == case["vars"]
            @test Set(p.edges) == Set(Tuple(e) for e in case["edges"])
            assigned = setdiff(reduce(∪, (d.c ∪ d.d for d in p.data); init = p.query.c ∪ p.query.d),
                               p.intervention_vars)
            @test Set(p.labels[collect(p.context_vars)]) ==
                  Set(case["con_vars_internal"]) ∪ Set(p.labels[collect(assigned)])
            @test p.labels[collect(p.intervention_vars)] == case["intv_vars"]
            label_map, local_csi = r_contexts(case)
            @test [(Int(cs.vars.bits), [(Int(c.zero.bits), Int(c.one.bits), c.equivalence, c.removed)
                                        for c in cs.contexts]) for cs in p.context_sets] == label_map
            @test [(l.x.bits, l.y.bits, l.z.bits, l.zero.bits, l.one.bits) for l in p.local_csi] == local_csi
            @test [[d.a.bits, d.b.bits, d.c.bits, d.d.bits] for d in p.data] == case["p_list"]
            @test [p.query.a.bits, p.query.b.bits, p.query.c.bits, p.query.d.bits] == case["q_vec"]
            # the search
            kwargs = (heuristic = case["heuristic"], cache = case["cache"],
                      rules = get(case, "rules", nothing), draw_derivation = true)
            io = IOBuffer()
            result = dosearch(p; kwargs..., verbose = haskey(case, "verbose"), io)
            expected = haskey(SAME_RESULT_AS, case["name"]) ?
                       only(c for c in LDAG_FIXTURES if c["name"] == SAME_RESULT_AS[case["name"]]) : case
            @test result.identifiable == expected["identifiable"]
            @test something(result.formula, "") == expected["formula"]
            haskey(expected, "derivation") && @test result.dot == expected["derivation"]
            haskey(case, "verbose") && @test strip(String(take!(io))) == strip(case["verbose"])
            result.identifiable && check_derivation(result)
        end
    end
end

@testset "LDAG interface" begin
    graph = "X -> Y : A = 1\nA -> Y"
    p = DosearchProblem("p(A,Y)", "p(Y|X,A=1)", graph)
    @test p isa LdagProblem
    @test dosearch(p).formula == "p(Y|A = 1)"
    @test dosearch("p(A,Y)", "p(Y|X,A=1)", graph).formula == "p(Y|A = 1)"
    @test occursin("query: p(Y|X,A = 1)", sprint(show, MIME"text/plain"(), p))
    @test CsiRule(5) == Dosearch.CsiCaseByCaseZero
    @test Int.(default_csi_rules()) == [0, 1, -3, 3, -4, 4, -5, 5, 6, -7, 7, -8, 8, -2, 2]
    out = dosearch("p(x,z=0) \n p(x,z=1)", "p(x,z)", "x -> y : w = 0 \n w -> y";
                   rules = [Dosearch.CsiCaseByCaseZero], warn = false)
    @test out.formula == "p(x,z)"
    @test out.rules == [Dosearch.CsiCaseByCaseZero]
    @test_throws ArgumentError dosearch(p; rules = [99])
    @test_throws ArgumentError dosearch(p; rules = [-6])
    # context variables beyond those in labels
    data = "p(rx, ry)\np(x, y, rx = 1, ry = 1)\np(x, rx = 1, ry = 0)\np(y, rx = 0, ry = 1)"
    q = "p(y | x, rx = 1, ry = 0)"
    g = "x -> ry\nx -> y\nrx -> ry"
    @test dosearch(LdagProblem(data, q, g; context_vars = ["ry", "rx"])).formula == "p(y|ry = 1,x,rx = 1)"
    @test_throws ArgumentError LdagProblem(data, q, g; context_vars = ["nope"])
    # Intentional difference: variables fixed to a value in the query or an
    # input are context variables without `context_vars` (R: not identifiable)
    @test Set(LdagProblem(data, q, g).labels[collect(LdagProblem(data, q, g).context_vars)]) == Set(["rx", "ry"])
    @test dosearch(LdagProblem(data, q, g)).formula == "p(y|ry = 1,x,rx = 1)"
    # a query that fixes an intervened variable that is in no label (R: not
    # identifiable, although p(Y|C=0,I_C=1) = p(Y) because C affects nothing)
    g = "I_C -> C\nX -> C : I_C = 1\nX -> Y"
    @test dosearch("p(X,C,Y)", "p(Y|C=0,I_C=1)", g).formula == "p(Y)"
    @test dosearch("p(X,C)\np(Y|X,C=0)", "p(Y|C=0,I_C=1)", g).formula == "\\sum_{X}\\left(p(X)p(Y|C = 0,X)\\right)"
    @test dosearch("p(X,C)\np(Y|X,C=0)", "p(Y|C,I_C=1)", g).identifiable
end

@testset "label alternatives on separate lines" begin
    # An edge listed on several lines has the alternatives of all its labels,
    # the same as alternatives separated by `;` on one line.
    flat(p) = [(cs.vars, [(c.zero, c.one, c.equivalence, c.removed) for c in cs.contexts]) for cs in p.context_sets]
    local_csi(p) = [(l.x, l.y, l.z, l.zero, l.one) for l in p.local_csi]
    one_line = """
        I_X -> X
        X -> Y : A = 1
        L -> X : I_X = 1; A = 0
        L -> Y
        A -> X : I_X = 1
        A -> Y
    """
    lines = """
        I_X -> X
        X -> Y : A = 1
        L -> X : I_X = 1
        L -> X : A = 0
        L -> Y
        A -> X : I_X = 1
        A -> Y
    """
    a = LdagProblem("p(A,X,Y)", "p(Y|X,I_X=1)", one_line)
    b = LdagProblem("p(A,X,Y)", "p(Y|X,I_X=1)", lines)
    @test a.labels == b.labels && a.edges == b.edges
    @test flat(a) == flat(b) && local_csi(a) == local_csi(b)
    @test dosearch(a).formula == dosearch(b).formula
    # alternatives that together cover every context remove the edge
    c = LdagProblem("p(A,X,Y)", "p(Y|X)", "A -> Y\nX -> Y : A = 0\nX -> Y : A = 1")
    @test c.labels[collect(first.(c.edges))] == ["A"]
end

@testset "CSI with context variables outside labels" begin
    # Intentional difference: R reasons by cases over the label variable rx
    # in a context that also assigns ry, a context variable in no label. It
    # finds no contexts and accepts any statement, here that y is independent
    # of its parent x and of its parent U, and identifies the query by
    # \left[p(y|ry = 1,rx)\right]_{rx = 1}. Dosearch.jl does not.
    graph = "x -> ry\nx -> y\nrx -> ry\nU -> ry\nU -> y\nV -> ry : rx = 0\nV -> x"
    data = "p(rx, ry)\np(x, y, rx = 1, ry = 1)\np(x, rx = 1, ry = 0)\np(y, rx = 0, ry = 1)"
    p = LdagProblem(data, "p(y | x, rx = 1, ry = 0)", graph; context_vars = ["ry", "rx"], warn = false)
    @test !dosearch(p).identifiable
    m = Dosearch.LdagModel(p)
    v(names...) = VarSet([findfirst(==(n), p.labels) for n in names])
    @test !Dosearch.separation_holds(m, Dosearch.CsiSeparation(v("x"), v("y"), v("ry", "rx"), VarSet(), v("ry")))
    @test !Dosearch.separation_holds(m, Dosearch.CsiSeparation(v("U"), v("y"), v("ry", "x", "rx"), v("ry"), VarSet()))
    # the contexts of label variables still apply: V -> ry vanishes where rx = 0
    @test Dosearch.separation_holds(m, Dosearch.CsiSeparation(v("V"), v("ry"), v("x", "rx"), v("rx"), VarSet()))
end

@testset "LDAG errors" begin
    msg(f) = try
        f(); "no error"
    catch e
        sprint(showerror, e)
    end
    ldag = "x -> y : z = 1 \n z -> y"
    @test occursin("bidirected edges are not supported for LDAGs", msg(() -> dosearch("p(x)", "p(y)", "x <-> y : y = 1")))
    @test occursin("x cannot appear in the label", msg(() -> dosearch("p(x)", "p(y)", "x -> y : x = 1")))
    @test occursin("y cannot appear in the label", msg(() -> dosearch("p(x)", "p(y)", "x -> y : y = 1")))
    @test occursin("duplicate assignment", msg(() -> dosearch("p(x)", "p(y)", "x -> y : z = 1, z = 0 \n z -> y")))
    @test occursin("only other parents of y may be assigned", msg(() -> dosearch("p(x)", "p(y)", "x -> y : z = 0")))
    # R ignores label values other than 0 and 1; Dosearch.jl rejects them.
    @test occursin("values must be 0 or 1", msg(() -> dosearch("p(x)", "p(y)", "x -> y : z = 2 \n z -> y")))
    @test occursin("malformed lines found", msg(() -> dosearch("p(x)", "p(y)", "x z w y :: z = 1")))
    @test occursin("the graph contains self-loops", msg(() -> dosearch("p(x)", "p(y)", "x -> x : y = 1")))
    for m in ["(", "p(", "p(x", "p(x|y", "p(x|do(x)", "p(x|do(z))"]
        @test occursin("Unable to parse input distribution", msg(() -> dosearch(m, "p(y)", ldag)))
    end
    @test occursin("duplicated variables", msg(() -> dosearch("p(x,x)", "p(y)", ldag)))
    @test occursin("Invalid value assignment", msg(() -> dosearch("p(x = 2)", "p(y)", ldag)))
    @test occursin("same variable on the left and right-hand side", msg(() -> dosearch("p(x|x)", "p(y)", ldag)))
    @test occursin("not supported for LDAGs", msg(() -> dosearch("p(x)", "p(y)", ldag; transportability = "z")))
    @test occursin("only supported for LDAGs", msg(() -> dosearch("p(x)", "p(y)", "x -> y"; context_vars = ["x"])))
end
