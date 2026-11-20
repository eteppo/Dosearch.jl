"""
    Dist(a, b, c, d)

A distribution term ``P(A \\mid do(C), B \\setminus C)`` over variable indices.

- `a`: the variables on the left-hand side.
- `b`: every variable on the right-hand side, both conditioning and
  intervened variables, so `c ⊆ b`.
- `c`: the intervened variables, inside `do(⋅)`.
- `d`: the missing data mechanisms (response indicators) that are enabled,
  i.e. fixed to the observed value, in this term.

`Dist` is an immutable value type and can be used as a `Dict` key.
"""
struct Dist
    a::VarSet
    b::VarSet
    c::VarSet
    d::VarSet
end

Base.hash(p::Dist, h::UInt) = hash(p.a.bits, hash(p.b.bits, hash(p.c.bits, hash(p.d.bits, h))))

"""
    Separation(x, y, given, intervened)

A d-separation statement checked by the do-calculus rules: every node of `x` is
d-separated from every node of `y` given `given` in the graph where the
incoming edges of `intervened` are removed. Node sets range over the variables
and their intervention nodes.
"""
struct Separation
    x::VarSet
    y::VarSet
    given::VarSet
    intervened::VarSet
end

const NO_SEPARATION = Separation(VarSet(), VarSet(), VarSet(), VarSet())
const NO_DIST = Dist(VarSet(), VarSet(), VarSet(), VarSet())

"""
    Inference{S}

The result of applying a rule to a term: the derived term `to`, the term
`required` that must already be known (`required.a` is empty if none), the
separation statement `sep::S` that must hold (`sep.x` is empty if none), and
whether the required term must be looked up with every candidate assignment
(`enumerate`).
"""
struct Inference{S}
    to::Dist
    required::Dist
    sep::S
    enumerate::Bool
end
