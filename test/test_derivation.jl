# The step-by-step derivation (`Derivation`, `latex`, `typst`) on examples. The
# derivations of all R fixtures are checked in test_fixtures.jl and
# test_ldag.jl with `check_derivation` (derivation_utils.jl).

const BACKDOOR_PROBLEM = ("p(x,y,z)", "p(y|do(x))", "x -> y\nz -> x\nz -> y")

@testset "back-door" begin
    r = dosearch(BACKDOOR_PROBLEM...)
    @test r.derivation isa Derivation
    @test latex(r) == latex(r.derivation) == """
        \\begin{align*}
        p(y \\mid \\mathrm{do}(x)) &= \\sum_{z} p(y, z \\mid \\mathrm{do}(x)) && \\text{marginalization} \\\\
        &= \\sum_{z} p(z \\mid \\mathrm{do}(x)) p(y \\mid \\mathrm{do}(x), z) && \\text{chain rule} \\\\
        &= \\sum_{z} p(z) p(y \\mid \\mathrm{do}(x), z) && \\text{rule 3: } (z \\perp\\!\\!\\!\\perp x)_{\\mathcal{G}_{\\overline{x}}} \\\\
        &= \\sum_{z} p(z) p(y \\mid x, z) && \\text{rule 2: } (y \\perp\\!\\!\\!\\perp x \\mid z)_{\\mathcal{G}_{\\underline{x}}}
        \\end{align*}
        """
    @test typst(r) == """
        \$ p(y | "do"(x)) &= sum_(z) p(y, z | "do"(x)) && quad "marginalization" \\
          &= sum_(z) p(z | "do"(x)) p(y | "do"(x), z) && quad "chain rule" \\
          &= sum_(z) p(z) p(y | "do"(x), z) && quad "rule 3:" (z ⫫ x)_(cal(G)_(overline(x))) \\
          &= sum_(z) p(z) p(y | x, z) && quad "rule 2:" (y ⫫ x | z)_(cal(G)_(underline(x))) \$
        """
    @test sprint(show, MIME"text/plain"(), r.derivation) == """
        p(y|do(x))
          = Σ_z p(y,z|do(x))             (marginalization)
          = Σ_z p(z|do(x)) p(y|do(x),z)  (chain rule)
          = Σ_z p(z) p(y|do(x),z)        (rule 3: (z ⊥⊥ x) in G_{x\u0305})
          = Σ_z p(z) p(y|x,z)            (rule 2: (y ⊥⊥ x | z) in G_{x\u0332})"""
    @test sprint(show, MIME"text/latex"(), r.derivation) == latex(r)
    @test sprint(show, r.derivation) == "Derivation(p(y|do(x)) = \\sum_{z}\\left(p(z)p(y|x,z)\\right))"
    # custom labels
    tex = latex(r; labels = Dict("x" => "X", "z" => "Z_1"))
    @test occursin("\\sum_{Z_1} p(Z_1) p(y \\mid X, Z_1)", tex)
    @test occursin("(y \\perp\\!\\!\\!\\perp X \\mid Z_1)_{\\mathcal{G}_{\\underline{X}}}", tex)
end

@testset "front-door" begin
    # The inner sum is over a second copy of x, which is renamed x'.
    r = dosearch("p(x,y,z)", "p(y|do(x))", "x -> z\nz -> y\nx <-> y")
    @test r.formula == "\\sum_{z}\\left(p(z|x)\\sum_{x}\\left(p(x)p(y|z,x)\\right)\\right)"
    @test sprint(show, MIME"text/plain"(), r.derivation) == """
        p(y|do(x))
          = Σ_z p(z,y|do(x))                           (marginalization)
          = Σ_z p(z|do(x)) p(y|do(x),z)                (chain rule)
          = Σ_z p(z|x) p(y|do(x),z)                    (rule 2: (z ⊥⊥ x) in G_{x\u0332})
          = Σ_z p(z|x) p(y|do(z,x))                    (rule 2: (y ⊥⊥ z | x) in G_{x\u0305,z\u0332})
          = Σ_z p(z|x) p(y|do(z))                      (rule 3: (y ⊥⊥ x | z) in G_{z\u0305,x\u0305})
          = Σ_z p(z|x) Σ_x' p(y,x'|do(z))              (marginalization)
          = Σ_z p(z|x) Σ_x' p(x'|do(z)) p(y|do(z),x')  (chain rule)
          = Σ_z p(z|x) Σ_x' p(x') p(y|do(z),x')        (rule 3: (x ⊥⊥ z) in G_{z\u0305})
          = Σ_z p(z|x) Σ_x' p(x') p(y|z,x')            (rule 2: (y ⊥⊥ z | x) in G_{z\u0332})"""
    tex = latex(r)
    @test endswith(tex, "&= \\sum_{z} p(z \\mid x) \\sum_{x'} p(x') p(y \\mid z, x') && " *
                        "\\text{rule 2: } (y \\perp\\!\\!\\!\\perp z \\mid x)_{\\mathcal{G}_{\\underline{z}}}\n\\end{align*}\n")
    @test occursin("(y \\perp\\!\\!\\!\\perp x \\mid z)_{\\mathcal{G}_{\\overline{z}\\,\\overline{x}}}", tex)
end

@testset "rule 3 with observations" begin
    # Rule 3 removes the edges into the nodes of Z that are not ancestors of
    # the observed W, Z(W) in Pearl's notation: here none, as x -> w.
    r = dosearch("p(y|w)", "p(y|do(x),w)", "x -> w\nw -> y")
    @test r.formula == "p(y|w)"
    j = Dosearch.justification(r.problem, Dosearch.causal_graph(r.problem), r.derivation.nodes,
                               r.derivation.nodes[1])
    @test j isa Dosearch.DoCalculusStep && j.rule == 3 && isempty(j.overline_z)
    @test occursin("(y ⊥⊥ x | w) in G)", sprint(show, MIME"text/plain"(), r.derivation))
end

@testset "known terms and inputs" begin
    # The formula keeps terms computed from an input by probability
    # manipulations; `inputs = true` rewrites them down to the inputs.
    r = dosearch("p(x,y)", "p(y|x)", "x -> y")
    @test r.formula == "p(y|x)"
    @test sprint(show, MIME"text/plain"(), r.derivation) == "p(y|x)"
    @test latex(r) == "\\begin{align*}\np(y \\mid x)\n\\end{align*}\n"
    @test typst(r) == "\$ p(y | x) \$\n"
    @test Dosearch.plain_text(r.derivation; inputs = true) ==
          "p(y|x)\n  = p(y,x) / (Σ_y' p(y',x))  (conditioning)"
    @test typst(r; inputs = true) ==
          "\$ p(y | x) &= frac(p(y, x), sum_(y') p(y', x)) && quad \"conditioning\" \$\n"
    r = dosearch("p(x,y)", "p(y|do(x))", "x -> y")
    @test Dosearch.plain_text(r.derivation; inputs = true) == """
        p(y|do(x))
          = p(y|x)                   (rule 2: (y ⊥⊥ x) in G_{x\u0332})
          = p(y,x) / (Σ_y' p(y',x))  (conditioning)"""
    # an input query
    r = dosearch("p(y|do(x))", "p(y|do(x))", "x -> y")
    @test sprint(show, MIME"text/plain"(), r.derivation) == "p(y|do(x))"
end

@testset "missing data" begin
    r = dosearch("p(x*,y*,r_x,r_y)\np(y)", "p(y|do(x))", "x -> y\ny -> r_y\nr_y -> r_x";
                 missing_data = "r_x : x, r_y : y")
    text = sprint(show, MIME"text/plain"(), r.derivation)
    @test startswith(text, "p(y|do(x))\n  = p(y|do(x,r_x = 1,r_y = 1))")
    @test endswith(text, "= (p(y) p(x|r_x = 1,y,r_y = 1)) / (Σ_y' p(y') p(x|r_x = 1,y',r_y = 1))" *
                         "  (rule 2: (x ⊥⊥ r_y | r_x,y) in G_{r\u0332_\u0332y\u0332})")
    tex = latex(r)
    @test occursin("(y \\perp\\!\\!\\!\\perp \\{r_{x}, r_{y}\\} \\mid x)_{\\mathcal{G}_{\\overline{x}\\,\\overline{r_{x}, r_{y}}}}", tex)
    # proxies are rewritten to their true variables where observed
    r = dosearch("p(x*,r_x)", "p(x|r_x = 1)", "x -> y"; missing_data = "r_x : x")
    @test r.identifiable
    text = Dosearch.plain_text(r.derivation; inputs = true)
    @test occursin("(x* = x if r_x = 1)", text)
    @test occursin("x^{*} = x \\text{ if } r_{x} = 1", latex(r; inputs = true))
    # rule 2 on an intervention that enables r_y: the rewritten term is
    # evaluated at r_y = 1, which the R formula leaves implicit
    r = dosearch("p(x,y*,r_y)", "p(y|do(x))", "x -> y\nx -> r_y"; missing_data = "r_y : y")
    @test r.formula == "p(y*|r_y,x)"
    text = sprint(show, MIME"text/plain"(), r.derivation)
    @test endswith(text, "= p(y*|r_y = 1,x)      (rule 2: (y* ⊥⊥ {r_y,x}) in G_{r̲_̲y̲,x̲})")
    @test occursin("&= p(y^{*} \\mid r_{y} = 1, x) &&", latex(r))
end

@testset "LDAGs" begin
    r = dosearch("p(A,Y)", "p(Y|X,A=1)", "X -> Y : A = 1\nA -> Y")
    @test sprint(show, MIME"text/plain"(), r.derivation) ==
          "p(Y|X,A = 1)\n  = p(Y|A = 1)  (CSI: (Y ⊥⊥ X | A = 1))"
    @test latex(r) == """
        \\begin{align*}
        p(Y \\mid X, A = 1) &= p(Y \\mid A = 1) && \\text{CSI: } (Y \\perp\\!\\!\\!\\perp X \\mid A = 1)
        \\end{align*}
        """
    r = dosearch("p(x|w) \n p(x,z=1)", "p(x,z=0)", "x -> y : z = 0 \n z -> y"; warn = false)
    @test sprint(show, MIME"text/plain"(), r.derivation) == """
        p(x,z = 0)
          = p(x) - p(x,z = 1)    (law of total probability)
          = p(x|w) - p(x,z = 1)  (CSI: (x ⊥⊥ w))"""
    graph = "I_X -> X\nX -> Y : A = 1\nL -> X : I_X = 1; A = 0\nL -> Y\nA -> X : I_X = 1\nA -> Y"
    r = dosearch("p(A,X,Y)", "p(Y|X,I_X=1)", graph)
    text = sprint(show, MIME"text/plain"(), r.derivation)
    @test occursin("  = Σ_A {p(Y,A = 0|X,I_X = 1) if A = 0; p(Y,A = 1|X,I_X = 1) if A = 1}", text)
    @test endswith(text, "= Σ_A {p(Y|X,A = 0) [p(A)]_{A = 0} if A = 0; p(Y|A = 1) [p(A)]_{A = 1} if A = 1}" *
                         "  (CSI: (A ⊥⊥ I_X))")
    tex = latex(r)
    @test occursin("\\sum_{A} \\begin{cases} p(Y \\mid X, A = 0) \\left[p(A)\\right]_{A = 0} & \\text{if } A = 0 \\\\ " *
                   "p(Y \\mid A = 1) \\left[p(A)\\right]_{A = 1} & \\text{if } A = 1 \\end{cases}", tex)
    @test occursin("(Y \\perp\\!\\!\\!\\perp I_{X} \\mid X, A = 0)", tex)
end

@testset "parentheses" begin
    p = DagProblem("p(x,y,z)", "p(y)", "x -> y\nz -> y")
    nodes = [Dosearch.Node(Dist(VarSet([i]), VarSet(), VarSet(), VarSet()), Dosearch.InputDistribution,
                           0, 0, true, Dosearch.NO_SEPARATION) for i in 1:3]
    d = Dosearch.Derivation(p, nodes, Dosearch.Formula[Dosearch.TermFormula(i) for i in 1:3], "")
    t(i) = Dosearch.TermFormula(i)
    S(f) = Dosearch.SumFormula(VarSet([2]), f)
    P(f, g) = Dosearch.ProductFormula(f, g)
    D(f, g) = Dosearch.DifferenceFormula(f, g)
    F(f, g) = Dosearch.FractionFormula(f, g)
    text(f) = Dosearch.markup(Dosearch.TextNotation(), d, f)
    labels = p.labels
    @test text(P(t(1), S(t(2)))) == "p($(labels[1])) Σ_$(labels[2]) p($(labels[2]))"
    @test text(P(S(t(2)), t(1))) == "(Σ_$(labels[2]) p($(labels[2]))) p($(labels[1]))"
    @test text(P(P(t(1), S(t(2))), t(3))) == "p($(labels[1])) (Σ_$(labels[2]) p($(labels[2]))) p($(labels[3]))"
    @test text(S(D(t(1), t(2)))) == "Σ_$(labels[2]) (p($(labels[1])) - p($(labels[2])))"
    @test text(D(t(1), D(t(2), t(3)))) == "p($(labels[1])) - (p($(labels[2])) - p($(labels[3])))"
    @test text(D(S(t(2)), t(1))) == "(Σ_$(labels[2]) p($(labels[2]))) - p($(labels[1]))"
    @test text(P(F(t(1), t(2)), t(3))) == "(p($(labels[1])) / p($(labels[2]))) p($(labels[3]))"
    @test text(F(P(t(1), t(2)), t(3))) == "(p($(labels[1])) p($(labels[2]))) / p($(labels[3]))"
    latex_(f) = Dosearch.markup(Dosearch.LatexNotation(Dict{String,String}()), d, f)
    @test latex_(P(D(t(1), t(2)), t(3))) == "\\left(p($(labels[1])) - p($(labels[2]))\\right) p($(labels[3]))"
end

@testset "labels" begin
    L = Dosearch.LatexNotation(Dict{String,String}())
    T = Dosearch.TypstNotation(Dict{String,String}())
    @test Dosearch.label(L, "x") == "x"
    @test Dosearch.label(L, "x", 1) == "x'"
    @test Dosearch.label(L, "age") == "\\mathit{age}"
    @test Dosearch.label(L, "x_1") == "x_{1}"
    @test Dosearch.label(L, "I_age") == "I_{\\mathit{age}}"
    @test Dosearch.label(L, "x*") == "x^{*}"
    @test Dosearch.label(L, "x*", 2) == "x^{*\\prime\\prime}"
    @test Dosearch.label(L, "x'", 1) == "x''"
    @test Dosearch.label(L, "a&b") == "\\mathit{a\\&b}"
    @test Dosearch.label(L, "12") == "12"
    # trailing Unicode subscript digits
    @test Dosearch.label(L, "A₁") == "A_{1}"
    @test Dosearch.label(L, "L₁₂") == "L_{12}"
    @test Dosearch.label(L, "A₁", 1) == "A'_{1}"
    @test Dosearch.label(L, "age₀") == "\\mathit{age}_{0}"
    @test Dosearch.label(L, "₁") == "\\mathit{₁}"
    @test Dosearch.label(L, "R_Y₁") == "R_{Y_{1}}"     # not a double subscript
    @test Dosearch.label(L, "R_Y₁", 1) == "R'_{Y_{1}}"
    @test Dosearch.label(L, "Y₁*") == "Y_{1}^{*}"
    @test Dosearch.label(T, "R_Y₁") == "R_(Y_(1))"
    @test Dosearch.label(T, "A₁") == "A_(1)"
    @test Dosearch.label(Dosearch.TextNotation(), "A₁") == "A₁"
    @test Dosearch.label(T, "age") == "italic(\"age\")"
    @test Dosearch.label(T, "x_1", 1) == "x'_(1)"
    @test Dosearch.label(T, "x*", 1) == "x^(*')"
    @test Dosearch.label(T, "a\"b") == "italic(\"a\\\"b\")"
    U = Dosearch.LatexNotation(Dict("x" => "X^{a}"))
    @test Dosearch.label(U, "x", 1) == "{X^{a}}'"
    # A renamed bound variable avoids the names of other variables.
    r = dosearch("p(x,y,x')", "p(y|x)", "x -> y\nx' -> y")
    @test occursin("Σ_y' p(y',x)", Dosearch.plain_text(r.derivation; inputs = true))
    r = dosearch("p(x,y,y')", "p(y|x)", "x -> y\ny' -> y")
    @test occursin("Σ_y'' p(y'',x)", Dosearch.plain_text(r.derivation; inputs = true))
end

@testset "errors" begin
    r = dosearch("p(x)", "p(y)", "x -> y")
    @test r.derivation === nothing
    @test_throws ArgumentError latex(r)
    @test_throws ArgumentError typst(r)
    r = dosearch(BACKDOOR_PROBLEM...; formula = false)
    @test r.identifiable && r.derivation === nothing && r.formula === nothing
    @test_throws ArgumentError latex(r)
end
