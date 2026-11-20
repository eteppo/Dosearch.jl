# Benchmarks of single searches, for optimizing the search core. Run in the dev
# session (BenchmarkTools comes from the global environment):
#
#     include("<repo>/benchmark/core.jl")
#     run_core_benchmarks(core_cases())
using Dosearch, Random, BenchmarkTools
include(joinpath(@__DIR__, "maximal.jl"))

# Add random confounders to `problem` until the query just becomes
# non-identifiable: a typical expensive search of the maximal graph search.
function just_nonidentifiable(problem, rng)
    chosen = BidirectedPair[]
    for c in shuffle(rng, confounding_candidates(problem))
        trial = vcat(chosen, c)
        p = with_candidates(problem, trial)
        dosearch(p; formula = false).identifiable || return p
        push!(chosen, c)
    end
    error("identifiable with every candidate")
end

function core_cases()
    probs10 = let rng = Xoshiro(10)
        [random_dag_problem(rng, 10) for _ in 1:3]
    end
    nonid10 = just_nonidentifiable(probs10[2], Xoshiro(1))
    example3 = DosearchProblem(
        "p(z_2,x_2|do(x_1))\np(z_1|x_2,do(x_1,y))\np(x_1|w_1,do(x_2))\np(y|z_1,z_2,x_1,do(x_2))\np(w|y,x_1,do(x_2))",
        "p(y,x_1|w,do(x_2))",
        "x_1 -> z_2\nx_1 -> z_1\nx_2 -> z_1\nx_2 -> z_2\nz_1 -> y\nz_2 -> y\nx_1 -> w\nx_2 -> w\nz_1 -> w\nz_2 -> w")
    example6 = DosearchProblem("p(x*,y*,z*,m_x,m_y,m_z)", "p(x,y,z)",
        "z -> x\nx -> y\nx -> m_z\ny -> m_z\ny -> m_x\nz <-> y"; missing_data = "m_x : x, m_y : y, m_z : z")
    chain8 = DosearchProblem("p(a,b,c,d,e,f,g,h)", "p(h|do(a))",
        "a -> b\nb -> c\nc -> d\nd -> e\ne -> f\nf -> g\ng -> h\na <-> h\nb <-> g\nc <-> f")
    return [("example3", example3), ("example6", example6), ("chain8", chain8), ("nonid10", nonid10)]
end

function run_core_benchmarks(cases)
    for (name, p) in cases
        b = @benchmark Dosearch.solve($p) samples = 5 evals = 1 seconds = 30
        println(rpad(name, 10), "min ", round(minimum(b).time / 1e6, digits = 1), " ms, median ",
                round(median(b).time / 1e6, digits = 1), " ms, ", minimum(b).allocs, " allocs")
    end
end
