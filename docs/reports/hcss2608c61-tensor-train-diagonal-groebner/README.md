# Tensor-train diagonal Gröbner obstruction

This report independently checks the proposed refutation of Conjecture 6.1 of
Hoşten, Chakrabarty Paul, Schmidt and Skurt, *Equations of Tree Tensor Network
Varieties*, arXiv:2608.19071v1 (19 August 2026), for the page arrangement and
lexicographic order in Definition 6.3. It contains computational evidence, not
a Lean proof or a claim about arbitrary other term orders.

Run with Python 3 and SymPy (verified with SymPy 1.14.0):

    python3 check.py

The stored output is `output.txt`; success ends with `ALL_CHECKS_OK`. All
mathematical calculations use integers, exact SymPy polynomials over ZZ, and
Fraction. Every flattening is reconstructed from coordinate tuples ordered as
in the source: rows (j₂,…,jᵢ,j₁), columns (jᵢ₊₁,…,jₙ₋₁,jₙ), with variables ordered
by (j₂,…,jₙ₋₁,j₁,jₙ).

The fidelity checks enumerate minors and verify the diagonal leading-term
property (Proposition 6.4). They reproduce Theorem 6.5 in four order-three
instances: dimensions/ranks (2,2,2)/(1,1), (3,2,3)/(2,2), (3,2,2)/(1,2), and
(2,2,3)/(2,1). The complete selected generator set passes every non-coprime
S-pair; the coprime pairs satisfy Buchberger's product criterion. Every
original generator reduces to zero by the selected set.

For dimensions (3,2,2,3), all bond ranks two, the program constructs a cut-two
minor F and a cut-three minor G from those flattenings and checks the exact
identity R = β₂₃ F − δG. It verifies the 16 distinct terms, their coefficients
±1, degree five, and the leading monomial. It enumerates all 840 cubic minors,
checks their diagonal initial terms, and verifies that none divides any term
of R. A separate longest-chain computation records all 48 term/cut bounds.

Padding is tested for ten dimension tuples with 4 ≤ N ≤ 8, including enlarged
end and interior dimensions. These are finite tests of the proposed general
padding argument, not a computational proof for unbounded N or dimensions.
The rank-two example is checked through its two independent nonzero rows and
columns. The general argument and the distinction between the generated ideal
and the vanishing ideal belong to the preregistration draft.
