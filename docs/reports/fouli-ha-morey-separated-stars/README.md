# Separated stars and an extra associated prime

`check.py` uses exact SymPy polynomial arithmetic and Groebner elimination. Run
`python3 docs/reports/fouli-ha-morey-separated-stars/check.py` from the repository
root. Its saved standard output is `output.txt`; success ends in `ALL_CHECKS_OK`.

The fidelity check reproduces Example 3.6 of Louiza Fouli, Tài Huy Hà and Susan
Morey, *Regular sequences of linear forms on monomial ideals*,
[arXiv:2606.15125v1](https://arxiv.org/abs/2606.15125v1): for
`I=(ab,bc,be,de,ef)` and `J=I+(a+b)`, the ideal quotient `J:b` is
`(a,b,c,e)`, an embedded associated prime extending `(a,c,e)`. The smaller
prime `(a,b,e)` is minimal.

For each `2 <= r,s <= 4`, over the rationals and the fields of 2, 3 and 5
elements, the program constructs the edge ideal of `K_{r,s}` and the star
`f=x1+y1+...+ys`. It checks the defining star condition, the two minimal
vertex covers, `I=A intersect B`, a nonzero dual-number image of `x1`, and the
exact ideal quotient `(I+(f)):x1=m`. Evaluation maps show that neither `A+(f)`
nor `B+(f)` equals `m`. It also checks `I:f=I` and the ideal identity obtained
after eliminating `x1`: `YT=Y intersect T intersect n^2`, including all three
irredundancy witnesses. Ideal intersections are computed by elimination;
`J:z` is computed as `(J intersect (z))/z`, with exact division.

The product-ideal primary-decomposition formula is prior art:
Aldo Conca and Jürgen Herzog, *Castelnuovo-Mumford regularity of products of
ideals*, [arXiv:math/0210065](https://arxiv.org/abs/math/0210065), Lemma 3.2.
The proposed application is the edge-ideal star quotient and the resulting
failure of Conjecture 3.7. The finite checks do not prove the unbounded family
or substitute for a Lean proof. The characteristic-free family argument uses
the socle element `[x1]`, the dual-number map, and the two original minimal
primes. Python and SymPy versions appear in the saved output.
