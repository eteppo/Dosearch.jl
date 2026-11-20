# Parsing of the text syntax of the R package: graphs, distributions and the
# lists of special nodes. Port of the string handling in `internal_dag.R` and
# `internal.R`. The results are still in terms of variable names; numbering
# and validation happen in problem.jl.

# Comma-and-space separated list, as in R's internal `cs`.
cs(xs) = join(xs, ", ")

# Text inputs list one item per line. A backslash followed by `n` is not a line
# break: it comes from raw strings, or from text copied from code, where the
# escape sequence was not interpreted.
function check_line_breaks(text::AbstractString, what::AbstractString, item::AbstractString)
    if occursin(r"\\[nr]", text)
        throw(ArgumentError("The $what contains the characters `\\n`, which are not a line break. " *
                            "Write each $item on its own line, for example in a multi-line string."))
    end
    return nothing
end

const GRAPH_SYNTAX = "Write one edge per line, as `x -> y` or `x <-> y`, with spaces around the arrow."
const DISTRIBUTION_SYNTAX = "A distribution is written `p(y|do(x),z)`: comma separated variables, " *
                            "the conditioning variables after `|` and the intervened ones in `do(...)`."

"""
    ParsedGraph

Edges of a graph given in the text syntax, in the order of the R parser:
`directed` holds `(from, to)` pairs and `bidirected` holds `(x, y)` pairs with
`x` sorted before `y`.
"""
struct ParsedGraph
    directed::Vector{Tuple{String,String}}
    bidirected::Vector{Tuple{String,String}}
end

"""
    parse_graph(graph::AbstractString) -> ParsedGraph

Parse a graph with one edge per line: `x -> y` for a directed edge and
`x <-> y` (or `x -- y`) for a bidirected edge. The three parts of a line are
separated by whitespace. Throws an `ArgumentError` for malformed lines, unknown
edge types, self-loops and directed cycles.
"""
function parse_graph(graph::AbstractString)
    text = replace(graph, "<->" => "--")
    lines = String[strip(l) for l in unique(split(text, r"\r|\n"))]
    tokens = [split(l) for l in lines]
    malformed = [!isempty(t) && length(t) != 3 for t in tokens]
    if any(malformed)
        throw(ArgumentError("Invalid graph, malformed lines found: " * cs(lines[malformed]) * ". " * GRAPH_SYNTAX))
    end
    keep = [length(t) == 3 for t in tokens]
    lines = lines[keep]
    tokens = tokens[keep]
    bad = unique([v for t in tokens for v in (t[1], t[3]) if occursin(r"[(),|=]", v)])
    if !isempty(bad)
        throw(ArgumentError("Invalid graph, variable names cannot contain `(`, `)`, `,`, `|` or `=`: " * cs(bad)))
    end
    edges = [t[2] for t in tokens]
    unknown = [e != "->" && e != "--" for e in edges]
    if any(unknown)
        throw(ArgumentError("Invalid graph, unknown edge types found: " * cs(edges[unknown])))
    end
    isdir = edges .== "->"
    directed = Tuple{String,String}[(String(t[1]), String(t[3])) for t in tokens[isdir]]
    loops = [from == to for (from, to) in directed]
    if any(loops)
        throw(ArgumentError("Invalid graph, the graph contains self-loops: " * cs(lines[isdir][loops])))
    end
    is_acyclic(directed) || throw(ArgumentError("Invalid graph, the graph contains cycles."))
    # R sorts the endpoints of a bidirected edge with `sort`, whose order
    # depends on the collation locale. Here the order is by code point, which
    # equals R in the C locale.
    bidirected = unique(Tuple{String,String}[
        Tuple(sort([String(t[1]), String(t[3])])) for t in tokens[.!isdir]])
    bloops = [x == y for (x, y) in bidirected]
    if any(bloops)
        bad = [replace(l, "--" => "<->") for l in lines[.!isdir] if split(l)[1] == split(l)[3]]
        throw(ArgumentError("Invalid graph, the graph contains self-loops: " * cs(bad)))
    end
    return ParsedGraph(directed, bidirected)
end

# Kahn's algorithm on the directed edges.
function is_acyclic(edges::Vector{Tuple{String,String}})
    indegree = Dict{String,Int}()
    out = Dict{String,Vector{String}}()
    for (from, to) in edges
        indegree[to] = get(indegree, to, 0) + 1
        get!(indegree, from, 0)
        push!(get!(out, from, String[]), to)
    end
    queue = [v for (v, d) in indegree if d == 0]
    seen = 0
    while !isempty(queue)
        v = pop!(queue)
        seen += 1
        for w in get(out, v, String[])
            indegree[w] -= 1
            indegree[w] == 0 && push!(queue, w)
        end
    end
    return seen == length(indegree)
end

"""
    ParsedDist

A distribution `p(lhs | do(intervened), conditioning)` in terms of variable
names, with the names of the enabled missing data mechanisms in `enabled`.
`text` is the distribution as given, for error messages.
"""
struct ParsedDist
    lhs::Vector{String}
    conditioning::Vector{String}
    intervened::Vector{String}
    enabled::Vector{String}
    text::String
end

const NAME_FORBIDDEN = ('|', ',', '(', ')', '$')

function parse_names(s::AbstractString)
    names = split(s, ',')
    all(n -> !isempty(n) && !any(in(NAME_FORBIDDEN), n), names) || return nothing
    return String.(names)
end

# Split the right-hand side into conditioning and intervened variables. The
# do-operator may appear at most once, before or after the conditioning
# variables.
function parse_rhs(rhs::AbstractString)
    if startswith(rhs, "do(")
        close = findfirst(')', rhs)
        close === nothing && return nothing
        intervened = parse_names(rhs[4:prevind(rhs, close)])
        rest = rhs[nextind(rhs, close):end]
        if isempty(rest)
            conditioning = String[]
        elseif startswith(rest, ",")
            conditioning = parse_names(rest[2:end])
        else
            return nothing
        end
    elseif endswith(rhs, ")") && occursin(",do(", rhs)
        k = findlast(",do(", rhs)
        conditioning = parse_names(rhs[1:prevind(rhs, first(k))])
        intervened = parse_names(rhs[(last(k) + 1):prevind(rhs, lastindex(rhs))])
    else
        conditioning = parse_names(rhs)
        intervened = String[]
    end
    (conditioning === nothing || intervened === nothing) && return nothing
    return conditioning, intervened
end

"""
    parse_distribution(d, kind; missing_data, md_sym) -> ParsedDist

Parse a distribution of the form `p(y)`, `p(y|z)`, `p(y|do(x))`,
`p(y|z,do(x))` or `p(y|do(x),z)`, where `y`, `z` and `x` are comma separated
variable names. Whitespace is ignored. With `missing_data = true`, a missing
data mechanism can be enabled by an assignment `r_x = md_sym`. `kind` names the
distribution in error messages.
"""
function parse_distribution(d::AbstractString, kind::AbstractString; missing_data::Bool, md_sym::AbstractString)
    s = replace(d, r"\s+" => "")
    parse_error() = throw(ArgumentError("Unable to parse $kind: $d. " * DISTRIBUTION_SYNTAX))
    m = match(r"^[Pp]\((.*)\)$", s)
    m === nothing && parse_error()
    inner = something(m.captures[1])   # the group always matches
    count('|', inner) > 1 && parse_error()
    bar = findfirst('|', inner)
    if bar === nothing
        lhs = parse_names(inner)
        conditioning, intervened = String[], String[]
    else
        lhs = parse_names(inner[1:prevind(inner, bar)])
        rhs = parse_rhs(inner[nextind(inner, bar):end])
        rhs === nothing && parse_error()
        conditioning, intervened = rhs
    end
    lhs === nothing && parse_error()
    enabled = String[]
    for part in (lhs, conditioning, intervened)
        dups = unique(part[nonunique_mask(part)])
        if !isempty(dups)
            throw(ArgumentError("Invalid $kind: $d, duplicated variables $(cs(dups))."))
        end
        equals = findall(contains('='), part)
        isempty(equals) && continue
        if !missing_data
            hint = any(i -> startswith(part[i], "I_"), equals) ?
                   " Intervention nodes such as I_x belong to LDAGs, whose edges have labels; in a DAG, write do(x)." : ""
            throw(ArgumentError("Invalid $kind $d: value assignments are only supported for missing data mechanisms." * hint))
        end
        split_eq = [split(part[i], '=') for i in equals]
        symbols = unique([length(e) >= 2 ? e[2] : "" for e in split_eq])
        if length(symbols) > 1
            throw(ArgumentError("Invalid $kind $d: multiple symbols used for missing data mechanisms."))
        end
        if symbols[1] != md_sym
            throw(ArgumentError("Invalid $kind $d: invalid symbol used for a missing data mechanism."))
        end
        for (i, e) in zip(equals, split_eq)
            part[i] = String(e[1])
            push!(enabled, part[i])
        end
    end
    return ParsedDist(lhs, conditioning, intervened, enabled, String(d))
end

# Elements that equal an earlier element (R's `duplicated`).
function nonunique_mask(xs::Vector{String})
    seen = Set{String}()
    mask = falses(length(xs))
    for (i, x) in enumerate(xs)
        x in seen ? (mask[i] = true) : push!(seen, x)
    end
    return mask
end

"""
    roles_to_distribution(roles::NamedTuple) -> String

Convert the alternative distribution format of the R package to text. Each
variable has a role: `0` for the left-hand side, `1` for an intervened
variable and `2` for a conditioning variable, so `(y = 0, x = 1)` is
`p(y|do(x))`.
"""
function roles_to_distribution(roles::NamedTuple)
    vals = collect(values(roles))
    names = string.(collect(keys(roles)))
    if !all(v -> v isa Real && isfinite(v), vals)
        throw(ArgumentError("Invalid distribution format $roles: all role values must be non-missing and finite."))
    end
    if !all(in((0, 1, 2)), vals)
        throw(ArgumentError("Invalid variable roles in distribution format $roles: all role values must be either 0, 1 or 2."))
    end
    if all(>(0), vals)
        throw(ArgumentError("Invalid variable roles in distribution format $roles: at least one variable must have role value 0."))
    end
    roles = Int[v for v in vals]
    a = join(names[roles .== 0], ",")
    b = join(names[roles .== 1], ",")
    c = join(names[roles .== 2], ",")
    s = "p(" * a
    (!isempty(b) || !isempty(c)) && (s *= "|")
    !isempty(b) && (s *= "do(" * b * ")")
    !isempty(b) && !isempty(c) && (s *= ",")
    return s * c * ")"
end

"""
    parse_missing_data(spec) -> Vector{Tuple{String,String}}

Parse missing data mechanisms `"r_x : x, r_y : y"` into `(indicator, variable)`
pairs.
"""
function parse_missing_data(spec::AbstractString)
    pairs = filter(!isempty, [replace(s, r"\s+" => "") for s in split(spec, ',')])
    malformed = ArgumentError("Malformed missing data mechanisms.")
    isempty(pairs) && throw(malformed)
    out = Tuple{String,String}[]
    for p in pairs
        parts = split(p, ':')
        length(parts) == 2 && !isempty(parts[1]) && !isempty(parts[2]) || throw(malformed)
        push!(out, (String(parts[1]), String(parts[2])))
    end
    allunique(last.(out)) && allunique(first.(out)) ||
        throw(ArgumentError("Malformed missing data mechanisms: duplicated mechanisms."))
    return out
end

# Comma separated node list of transportability or selection bias nodes.
parse_node_list(spec::AbstractString) =
    unique(filter(!isempty, [replace(s, r"\s+" => "") for s in split(spec, ',')]))
