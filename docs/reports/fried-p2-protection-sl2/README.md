# Fried's CA-group protection target: exact SL2 checks

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored output is
`output.txt`; a successful run exits 0 and ends with `ALL_CHECKS_OK`.
All evidence uses integers, exact Smith normal form, and finite-field arithmetic.

The fidelity check reproduces the SL(2,3) row of Table 2 in Matthew Fried,
*Game Conductors of Finite Groups: Determinantal Torsion from Structured Payoff
Probes*, arXiv:2607.05698v2. It enumerates the 24 matrices, their actual
centralizers, and the full integer commuting row-difference matrix. Its nonzero
Smith factors are six copies of 1 and one 6, and its seven noncentral
centralizers give the stated conductor {2,3}.

The characteristic-two checks construct F4, F8 and F64 from the binary
polynomials X²+X+1, X³+X+1 and X⁶+X+1. Irreducibility, field operations and
inverses are checked exactly. All SL2 matrices are enumerated. Nonidentity
matrices are grouped by their projective coordinates (b,c,a+d) in Mat2(F)/FI;
every projective signature occurs. Each representative's full matrix commutant
has dimension two, and it commutes with every matrix in its signature block.
All ordered commuting pairs are additionally checked for q=4 and q=8.

For q=64 there are 262080 matrices and 4161 signatures: 65 centralizers of order
64, 2080 of order 63, and 2016 of order 65. The program checks displayed diagonal
commutators for every upper and lower transvection, and an explicit transvection
factorization of every enumerated matrix. The expected arithmetic refutation
uses 9 dividing 262080 and 3 not dividing 4160. The q=64 full commuting matrix and
all its ordered pairs are not constructed.

The additional ten samples check the arithmetic choice q=2^(p(p−1)k), for
p in {3,5,7,11,13} and k in {1,2}. These samples support the separate general
argument; they do not establish a statement for every prime or every k.
The repository settlement target is the arithmetic implication in Section 10,
with b counting actual maximal abelian subgroups. No Lean proof is included.
