# Benchmarks of the maximal identifiable graph search on random problems with
# several observational distributions. Run in the dev session:
#
#     include("<repo>/benchmark/maximal.jl")
#     problem = random_dag_problem(Xoshiro(9), 9)
#     measure_maximal(problem)
#
# Runs longer than a minute: `start_bench(path, problems, variants)` runs in a
# background task of the session and appends one line per problem to `path`.
using Dosearch, Random

function random_dag_problem(rng, n; p_edge = 0.35, nsmall = 2)
    while true
        vars = ["v$i" for i in 1:n]
        order = shuffle(rng, vars)
        edges = String[]
        for i in 1:n, j in (i + 1):n
            rand(rng) < p_edge && push!(edges, "$(order[i]) -> $(order[j])")
        end
        isempty(edges) && push!(edges, "$(order[1]) -> $(order[2])")
        x = order[rand(rng, 1:(n ÷ 2))]
        y = order[rand(rng, (n ÷ 2 + 1):n)]
        reach = Set([x])   # require a directed path from x to y
        for e in edges     # edges are listed in topological order
            a, _, b = split(e)
            a in reach && push!(reach, b)
        end
        y in reach || continue
        data = String[]
        big = shuffle(rng, vars)[1:rand(rng, (n - 1):n)]
        push!(data, "p(" * join(big, ",") * ")")
        for _ in 1:nsmall
            sub = shuffle(rng, vars)[1:rand(rng, 2:max(2, n - 2))]
            if length(sub) >= 3 && rand(rng) < 0.5
                c = rand(rng, 1:(length(sub) - 1))
                push!(data, "p(" * join(sub[1:c], ",") * "|" * join(sub[(c + 1):end], ",") * ")")
            else
                push!(data, "p(" * join(sub, ",") * ")")
            end
        end
        problem = DosearchProblem(join(data, "\n"), "p($y|do($x))", join(edges, "\n"); warn = false)
        dosearch(problem; formula = false).identifiable && return problem
    end
end

function measure_maximal(problem; use_certificates = true, parallel = Threads.nthreads() > 1, kwargs...)
    cands = confounding_candidates(problem)
    o = Dosearch.MonotoneOracle(problem, cands, Dosearch.search_options(; kwargs...); use_certificates, parallel)
    t = @elapsed begin
        sets = Dosearch.isidentifiable!(o, BitSet()) ? Dosearch.dualize_and_advance!(o) :
               (BitSet[], [BitSet()])
    end
    return (m = length(cands), maximal = length(sets[1]), minimal_nonid = length(sets[2]),
            undecided = count(s -> s in o.undecided, sets[2]),
            searches = o.searches, certified = o.certified, time = round(t, digits = 2),
            id_time = round(o.identifiable_time, digits = 2),
            nonid_time = round(o.nonidentifiable_time, digits = 2)), sets
end

# Run `measure_maximal` on `problems` for each keyword set in `variants`,
# appending results to `path`.
function start_bench(path, problems, variants)
    return Threads.@spawn try
        open(path, "a") do io
            for (k, p) in enumerate(problems), (name, kw) in variants
                res, _ = measure_maximal(p; kw...)
                println(io, "#$k $name n=$(length(p.labels)) ", res)
                flush(io)
            end
            println(io, "done")
        end
    catch e
        open(io -> println(io, "ERROR ", sprint(showerror, e)), path, "a")
    end
end
