# Checks of the step-by-step derivation of an identified query, run on every
# R fixture: the last step is the identifying formula of the R package, every
# do-calculus statement of the derivation holds in its mutilated graph by an
# independent d-separation test, and the derivation renders in every format.

# The terms that a formula consists of.
formula_terms(f::Dosearch.TermFormula) = [f.node]
formula_terms(f::Dosearch.Formula) =
    reduce(vcat, (formula_terms(getfield(f, k)) for k in fieldnames(typeof(f))
                  if getfield(f, k) isa Dosearch.Formula); init = Int[])

# d-separation of `x` and `y` given `c` in the DAG with the given parents, by
# moralizing the ancestral graph (Lauritzen et al., 1990).
function moral_dseparated(parents::Vector{Vector{Int}}, x, y, c)
    anc = Set{Int}(union(x, y, c))
    stack = collect(anc)
    while !isempty(stack)
        for u in parents[pop!(stack)]
            u in anc || (push!(anc, u); push!(stack, u))
        end
    end
    adj = Dict(v => Set{Int}() for v in anc)
    for v in anc
        pa = parents[v]
        for u in pa
            push!(adj[u], v); push!(adj[v], u)
            for w in pa
                u == w || push!(adj[u], w)
            end
        end
    end
    seen = Set{Int}(x)
    stack = collect(x)
    while !isempty(stack)
        for w in adj[pop!(stack)]
            (w in c || w in seen) && continue
            w in y && return false
            push!(seen, w); push!(stack, w)
        end
    end
    return true
end

# Pearl's statement of a do-calculus step in the explicitly mutilated graph,
# with a latent node for every bidirected edge.
function pearl_statement_holds(p::DagProblem, j::Dosearch.DoCalculusStep)
    n = length(p.labels)
    parents = [Int[] for _ in 1:n]
    for (from, to) in p.directed
        push!(parents[to], from)
    end
    for (a, b) in p.bidirected
        push!(parents, Int[])
        push!(parents[a], length(parents))
        push!(parents[b], length(parents))
    end
    for v in j.overline ∪ j.overline_z
        empty!(parents[v])
    end
    for v in 1:n
        filter!(u -> !(u in j.underline), parents[v])
    end
    return moral_dseparated(parents, collect(j.y), collect(j.z), collect(j.given))
end

function check_derivation(r::DosearchResult)
    d = r.derivation
    p = r.problem
    term(i) = Dosearch.term_string(p, d.nodes[i].dist)
    final(steps) = isempty(steps) ? Dosearch.TermFormula(1) : last(steps)[2]
    steps = Dosearch.rewrite_steps(d)
    @test Dosearch.r_formula(p, final(steps), term) == r.formula
    @test all(i -> d.nodes[i].primitive, formula_terms(final(steps)))
    full = Dosearch.rewrite_steps(d; inputs = true)
    @test all(i -> d.nodes[i].parent1 == 0, formula_terms(final(full)))
    @test length(full) >= length(steps)
    if p isa DagProblem
        g = Dosearch.causal_graph(p)
        for (j, _) in full
            just = Dosearch.justification(p, g, d.nodes, d.nodes[j])
            just isa Dosearch.DoCalculusStep && @test pearl_statement_holds(p, just)
        end
    end
    tex = latex(d)
    @test startswith(tex, "\\begin{align*}\n") && endswith(tex, "\\end{align*}\n")
    @test count("&=", tex) == length(steps)
    typ = typst(d; inputs = true)
    @test startswith(typ, "\$ ") && count("&=", typ) == length(full)
    @test count('\n', sprint(show, MIME"text/plain"(), d)) == length(steps)
end
