# Orzel's higher-characteristic Sym sufficiency question

Run `python3 check.py` with Python 3.12 or newer. The program uses only the
standard library and exact arithmetic in F3 and F9 = F3[i]/(i² + 1).
The recorded output is `check.out`.

Source: Ian Orzel, *Computing the Elementary Symmetric Polynomials in Positive
Characteristics*, arXiv:2509.05009v3, Section 4, third bullet, with Claim 3.7
and Corollary 3.4. Source fidelity checks reproduce the degree-four Newton
identity in Appendix A and Lemma 3.1's free addition of fourth powers.

For A = Σ xᵢ² and B = Σ yⱼ², the checks verify symbolic Hessian formulas,
the decomposition into two nonzero reducible quartics, and explicit expressions
AB = e₄(L) − e₄(M) over F9. Exhaustive F3 point checks cover (a,b) = (4,4),
(4,5), (5,4), with 45,927 points in total. Another 600 seeded F9 points check
the extension-field rank formulas. Explicit isotropic vectors span the
quadratic nullcones in dimensions four through six.

These computations support the proposed counterfamily. They do not establish
nonmembership in Sym over an algebraic closure: that requires the Hessian-locus
dimension and irreducibility argument in the preregistration. No Lean proof
is included. The predicted minimum Sym top fan-in over an algebraic closure is
two; this program verifies the upper bound in the listed instances.
