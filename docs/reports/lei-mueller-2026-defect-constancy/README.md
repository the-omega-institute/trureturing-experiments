# Lei–Müller defect constancy: exact checks

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored output is
`check.out`; success requires exit 0 and the final line `ALL_CHECKS_OK`.
All matrix entries, characteristic polynomials, polynomial divisions and ranks
are exact. Integer matrix rank is calculated over its fraction field, the
rationals. Smith normal form independently checks ranks for order at most 12.
The free integer rank of the cokernel of a square matrix M is its order minus
its rank. No floating point eigenvalue approximations are used.

The source is Antonio Lei and Katharina Müller, *Iwasawa theory for abelian
towers of digraphs*, arXiv:2601.19571v1, Conjecture 7.11. The spectral definitions
of a and b are also explicit in Kundu–Müller–Newton–Stewart–Wang,
arXiv:2609.16702v2, Definition 2.13. Adjacency entries count arrows from the
column vertex to the row vertex.

Fidelity reproduces the connected case r=3, p=2 of the original Example 7.10
at levels 0 through 3: the characteristic polynomial is x^(3*2^m)-1,
a=b=1, and defect zero. The original Example 7.2 has a discrepancy: its printed
matrix has characteristic polynomial (x-1)(x^2-4x+2), giving a=b=1 and defect
zero, whereas the text states a=2 and defect one. That example is not used as
evidence for the proposed counterexample.

The proposed family is A=cJ+(I+uv^T) tensor C_q, with q=p^h,
u=(1,-1,0)^T, v=(1,1,-2)^T and c=3. For every combination p in {2,3},
h in {1,2}, levels 0 through h+1 give

    a = 3*gcd(q,p^m)-1
    b = 2*gcd(q,p^m)-1
    defect = gcd(q,p^m).

The base characteristic polynomial is checked independently. For each
combination the full first lift is constructed and its characteristic
polynomial and integer matrix rank are calculated directly. Other levels use
det(xI-A tensor C_s)=det(x^s I-A^s) and
nullity(I-A tensor C_s)=nullity(I-A^s), consequences of the cyclic block
recurrence. Reachability searches inspect every tested lift directly.
The first lift therefore supplies a direct finite refutation independent of
the formulas used to check higher levels.

The six-vertex-class construction produces zero-diagonal 0/1 matrices on
18q vertices. Exact matrix equalities check the quotient A and its invariant
sum-zero complement. Reachability is checked at the same levels as above,
and length-2 and length-3 closed walks are checked at every base vertex.
For q=2 the full 36-vertex graph and its 72-vertex first lift give respectively
(a,b,defect)=(4,3,1) and (14,12,2).

These are finite computational checks. They do not constitute a Lean proof,
certify the family at all parameters, or establish novelty. A mismatch in
the expected family law, a disconnected required lift or a failed direct
counterexample check rejects this proposed route.
