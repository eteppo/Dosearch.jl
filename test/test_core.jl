# Unit tests of the core data structures: variable sets, the d-separation
# graph and the heuristic frontier.

using Random

@testset "VarSet" begin
    s = VarSet([1, 3, 5])
    @test collect(s) == [1, 3, 5]
    @test length(s) == 3
    @test 3 in s && !(2 in s) && !(0 in s) && !(65 in s)
    @test s ∪ VarSet([2]) == VarSet(1:3) ∪ VarSet([5])
    @test s ∩ VarSet([3, 4]) == VarSet([3])
    @test setdiff(s, VarSet([1])) == VarSet([3, 5])
    @test VarSet([1, 5]) ⊆ s && !(VarSet([2]) ⊆ s)
    @test isempty(VarSet()) && !isempty(s)
    @test VarSet([1, 2]) << 2 == VarSet([3, 4])
    @test VarSet([3, 4]) >> 2 == VarSet([1, 2])
    @test Dosearch.fullset(3) == VarSet(1:3)
    @test Dosearch.fullset(64) == VarSet(typemax(UInt64))
    @test VarSet(UInt64(5)) == VarSet([1, 3])
    @test_throws ArgumentError VarSet([0])
    @test_throws MethodError VarSet(3)
    @test sprint(show, s) == "VarSet(1, 3, 5)"
    d = Dict(VarSet([1]) => 1)
    @test d[VarSet([1])] == 1
end

@testset "banker's sequence of subsets" begin
    # Order of get_subsets in set.cpp: by size, lexicographic within a size.
    @test Dosearch.subsets_by_size(3) ==
          [VarSet([1]), VarSet([2]), VarSet([3]), VarSet([1, 2]), VarSet([1, 3]),
           VarSet([2, 3]), VarSet([1, 2, 3])]
    for n in 1:8
        subs = Dosearch.subsets_by_size(n)
        @test length(subs) == 2^n - 1
        @test allunique(subs)
        @test issorted(subs; by = length)
        @test issorted(subs; lt = (a, b) -> length(a) == length(b) && collect(a) < collect(b) ||
                                            length(a) < length(b))
    end
end

# A literal port of dcongraph::dsep from dcongraph.cpp, with boolean matrices,
# to check the bit-row implementation against.
function naive_dsep(n, B, Ce, x, y, c, j)
    hh = [Ce[i, k] for i in 1:n, k in 1:n]
    th = [B[k, i] for i in 1:n, k in 1:n]
    tt = falses(n, n)
    inset(el, s) = (s >> (el - 1)) & 1 == 1
    adj() = hh[x, y] || hh[y, x] || th[x, y] || th[y, x] || tt[x, y] || tt[y, x]
    m = ((UInt64(1) << n) - 1) & ~c & ~(UInt64(1) << (x - 1)) & ~(UInt64(1) << (y - 1))
    for el in 1:n
        inset(el, j) || continue
        for i in 1:n
            hh[i, el] = false
            th[i, el] = false
        end
    end
    adj() && return false
    for el in 1:n
        inset(el, c) || continue
        for i in 1:n, k in 1:n
            hh[i, k] |= hh[i, el] & hh[el, k]
            th[i, k] |= th[i, el] & hh[el, k]
            tt[i, k] |= th[i, el] & th[k, el]
        end
        adj() && return false
    end
    for el in 1:n
        inset(el, m) || continue
        for i in 1:n, k in 1:n
            (el == i || el == k) && continue
            hh[i, k] |= (th[el, i] & th[el, k]) | (hh[i, el] & th[el, k]) |
                        (th[el, i] & hh[el, k]) | (hh[el, i] & hh[el, k] & tt[el, el])
            tt[i, k] |= (tt[i, el] & tt[el, k]) | (th[i, el] & tt[el, k]) |
                        (tt[i, el] & th[k, el]) | (th[i, el] & th[k, el] & tt[el, el])
            th[i, k] |= (tt[i, el] & th[el, k]) | (th[i, el] & th[el, k]) |
                        (tt[i, el] & hh[el, k]) | (th[i, el] & hh[el, k] & tt[el, el])
        end
        adj() && return false
    end
    return true
end

# A random separation statement of the form that the rules make: the
# intervened nodes are conditioned on, the queried nodes are not, and `y` may
# contain intervention nodes of variables that are not intervened on.
function rule_shaped_statement(rng, nv)
    vars = randperm(rng, nv)
    nx = rand(rng, 1:max(1, nv ÷ 3))
    x = vars[1:nx]
    rest = vars[(nx + 1):end]
    ny = rand(rng, 0:min(2, length(rest)))
    y = rest[1:ny]
    rest = rest[(ny + 1):end]
    given = filter(_ -> rand(rng) < 0.5, rest)
    intervened = filter(_ -> rand(rng) < 0.5, given)
    free = filter(v -> !(v in intervened), rest)
    if (isempty(y) || rand(rng, Bool)) && !isempty(free)
        push!(y, nv + rand(rng, free))   # an intervention node
    end
    isempty(y) && push!(y, nv + x[1])
    return Dosearch.Separation(VarSet(x), VarSet(y), VarSet(given), VarSet(intervened))
end

@testset "d-separation matches a literal port of dcongraph::dsep" begin
    rng = Xoshiro(1)
    for trial in 1:400
        nv = rand(rng, 2:7)
        g = Dosearch.DsepGraph(nv)
        n = 2nv
        B = falses(n, n)     # B[to, from], as in dcongraph
        Ce = falses(n, n)
        for x in 1:nv
            B[x, nv + x] = true
        end
        order = randperm(rng, nv)
        for a in 1:nv, b in (a + 1):nv
            if rand(rng) < 0.4
                Dosearch.add_directed!(g, order[a], order[b])
                B[order[b], order[a]] = true
            end
            if rand(rng) < 0.25
                Dosearch.add_bidirected!(g, order[a], order[b])
                Ce[order[a], order[b]] = Ce[order[b], order[a]] = true
            end
        end
        for q in 1:25
            s = rule_shaped_statement(rng, nv)
            expected = all(naive_dsep(n, B, Ce, x, y, s.given.bits, s.intervened.bits)
                           for x in s.x, y in s.y)
            @test Dosearch.dseparated(g, s) == expected
        end
    end
end

@testset "the rules make only rule-shaped separation statements" begin
    # The d-separation test agrees with the C++ code on these statements.
    rng = Xoshiro(3)
    n = 7
    checked = 0
    bad = 0
    for trial in 1:300
        md = rand(rng, Bool)   # variables 1:3 are (true, indicator, proxy)
        ctx = Dosearch.RuleContext(n, md, md ? VarSet([2]) : VarSet(), md ? VarSet([3]) : VarSet(),
                                   md ? VarSet([1]) : VarSet(), VarSet())
        roles = rand(rng, 0:3, n)   # 0: absent, 1: left, 2: conditioning, 3: intervened
        any(==(1), roles) || continue
        a = VarSet(findall(==(1), roles))
        c = VarSet(findall(==(3), roles))
        b = VarSet(findall(==(2), roles)) ∪ c
        p = Dist(a, b, c, VarSet())
        for r in instances(Rule), z in Dosearch.subsets_by_size(n)
            r == Dosearch.InputDistribution && continue
            inf = Dosearch.apply_rule(r, p, z, ctx)
            (inf === nothing || isempty(inf.sep.x)) && continue
            s = inf.sep
            checked += 1
            s.intervened ⊆ s.given && isempty(s.x ∩ s.given) && isempty(s.y ∩ s.given) &&
                isempty(s.x ∩ s.y) && !isempty(s.y) || (bad += 1)
        end
    end
    @test checked > 2_000
    @test bad == 0
end

@testset "separation folded into subset masks" begin
    # For rules 1-3 the search removes the variables d-connected to the
    # left-hand side from the subset mask instead of testing each subset.
    rng = Xoshiro(4)
    n = 6
    for trial in 1:200
        g = Dosearch.DsepGraph(n)
        order = randperm(rng, n)
        for a in 1:n, b in (a + 1):n
            rand(rng) < 0.4 && Dosearch.add_directed!(g, order[a], order[b])
            rand(rng) < 0.2 && Dosearch.add_bidirected!(g, order[a], order[b])
        end
        ctx = Dosearch.RuleContext(n, false, VarSet(), VarSet(), VarSet(), VarSet())
        roles = rand(rng, 0:3, n)
        any(==(1), roles) || continue
        a = VarSet(findall(==(1), roles))
        c = VarSet(findall(==(3), roles))
        b = VarSet(findall(==(2), roles)) ∪ c
        p = Dist(a, b, c, VarSet())
        connected = Dosearch.dconnected(g, a.bits, b.bits, c.bits)
        for r in (Dosearch.InsertObservation, Dosearch.ObservationToAction, Dosearch.InsertAction)
            mask = Dosearch.subset_mask(r, p, ctx)
            reduced = setdiff(mask, Dosearch.connected_variables(r, connected, n))
            passing = [z for z in Dosearch.subsets_by_size(n)
                       if (inf = Dosearch.apply_rule(r, p, z, ctx)) !== nothing && Dosearch.dseparated(g, inf.sep)]
            kept = [z for z in Dosearch.subsets_by_size(n)
                    if z ⊆ reduced && Dosearch.apply_rule(r, p, z, ctx) !== nothing]
            @test passing == kept
        end
    end
end

@testset "d-separation basics" begin
    sep(x, y, given = Int[], intervened = Int[]) =
        Dosearch.Separation(VarSet(x), VarSet(y), VarSet(given), VarSet(intervened))
    # chain x -> z -> y
    g = Dosearch.DsepGraph(3)
    Dosearch.add_directed!(g, 1, 3)
    Dosearch.add_directed!(g, 3, 2)
    @test !Dosearch.dseparated(g, sep([1], [2]))
    @test Dosearch.dseparated(g, sep([1], [2], [3]))
    # collider x -> z <- y, and a descendant w of the collider
    g = Dosearch.DsepGraph(4)
    Dosearch.add_directed!(g, 1, 3)
    Dosearch.add_directed!(g, 2, 3)
    Dosearch.add_directed!(g, 3, 4)
    @test Dosearch.dseparated(g, sep([1], [2]))
    @test !Dosearch.dseparated(g, sep([1], [2], [3]))
    @test !Dosearch.dseparated(g, sep([1], [2], [4]))
    # x <-> y <-> z: conditioning on the collider y connects x and z, and
    # intervening on y removes its bidirected edges
    g = Dosearch.DsepGraph(3)
    Dosearch.add_bidirected!(g, 1, 2)
    Dosearch.add_bidirected!(g, 2, 3)
    @test Dosearch.dseparated(g, sep([1], [3]))
    @test !Dosearch.dseparated(g, sep([1], [3], [2]))
    @test Dosearch.dseparated(g, sep([1], [3], [2], [2]))
    # intervention node I(y) of y in x <- y
    g = Dosearch.DsepGraph(2)
    Dosearch.add_directed!(g, 2, 1)
    @test !Dosearch.dseparated(g, sep([1], [4]))
    @test Dosearch.dseparated(g, sep([1], [4], [2]))
end

@testset "heuristic frontier matches libstdc++ priority_queue" begin
    # heap_ops.txt is the output of fixtures/heap_reference.cpp compiled with
    # GCC: pushes with many tied scores and the ids that pop() returned.
    h = Dosearch.HeapFrontier()
    pops = 0
    for line in eachline(joinpath(@__DIR__, "fixtures", "heap_ops.txt"))
        parts = split(line)
        if parts[1] == "push"
            Dosearch.heap_push!(h, parse(Int, parts[2]), parse(Int, parts[3]))
        else
            @test Dosearch.heap_pop!(h) == parse(Int, parts[2])
            pops += 1
        end
    end
    @test pops == 200
end

@testset "deposit" begin
    mask = VarSet([2, 5, 7])
    @test Dosearch.deposit(VarSet([1]), mask) == VarSet([2])
    @test Dosearch.deposit(VarSet([2, 3]), mask) == VarSet([5, 7])
    @test Dosearch.deposit(VarSet(1:3), mask) == mask
    @test Dosearch.deposit(VarSet(), mask) == VarSet()
    # Subsets of a mask in banker's order are the subsets of the mask in the
    # banker's order of all variables, which the search relies on.
    every = Dosearch.subsets_by_size(7)
    @test [Dosearch.deposit(r, mask) for r in Dosearch.subsets_by_size(3)] == filter(z -> z ⊆ mask, every)
end

@testset "subset masks contain every subset a rule accepts" begin
    rng = Xoshiro(2)
    rules = setdiff(instances(Rule), [Dosearch.InputDistribution])
    n = 7
    subsets = Dosearch.subsets_by_size(n)
    for trial in 1:400
        # variables 1:3 are (true, indicator, proxy) of a missing data variable
        md = rand(rng, Bool)
        ctx = Dosearch.RuleContext(n, md, md ? VarSet([2]) : VarSet(), md ? VarSet([3]) : VarSet(),
                                   md ? VarSet([1]) : VarSet(), rand(rng) < 0.3 ? VarSet([n]) : VarSet())
        roles = rand(rng, 0:3, n)   # 0: absent, 1: left, 2: conditioning, 3: intervened
        a = VarSet(findall(==(1), roles))
        c = VarSet(findall(==(3), roles))
        b = VarSet(findall(==(2), roles)) ∪ c
        d = VarSet(filter(_ -> rand(rng, Bool), collect((a ∪ b) ∩ ctx.md_s)))
        p = Dist(a, b, c, d)
        for r in rules
            mask = Dosearch.subset_mask(r, p, ctx)
            @test all(z -> z ⊆ mask || Dosearch.apply_rule(r, p, z, ctx) === nothing, subsets)
        end
    end
end

@testset "rules" begin
    @test Rule(-2) == Dosearch.ActionToObservation
    @test Int(Dosearch.ProductMarginal) == -6
    @test Dosearch.isdocalculus(Dosearch.DeleteAction)
    @test !Dosearch.isdocalculus(Dosearch.Marginalize)
    @test !Dosearch.isdocalculus(Dosearch.InputDistribution)
    @test Int.(default_rules(missing_data = false, improve = true)) == [4, 5, -2, -3, 2, 3, 6, -6]
    @test Int.(default_rules(missing_data = false, improve = false)) == [4, 5, -1, -2, -3, 1, 2, 3, 6, -6]
    @test Int.(default_rules(missing_data = true, improve = true)) ==
          [4, 5, 9, 10, -2, -3, 2, 3, 6, -6, 7, -7, 8, -8]
    @test Int.(default_rules(missing_data = true, improve = false)) ==
          [4, 5, 9, 10, -1, -2, -3, 1, 2, 3, 6, -6, 7, -7, 8, -8]
    @test Dosearch.rule_label.(Rule.([1, -2, 3, 4, 5, -6, 7, -8, 9, 10])) ==
          ["R1", "R2", "R3", "M", "C", "P", "D", "D", "A", "EX"]
end
