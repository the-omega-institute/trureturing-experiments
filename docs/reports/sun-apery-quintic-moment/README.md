# Sun's quintic squared moment

The source is Zhi-Hong Sun, *Generalizations of the Christoffel-Darboux
formula and congruences involving Apéry-like numbers*,
[arXiv:2608.13192v1](https://arxiv.org/abs/2608.13192v1), Conjecture 4.1.
The lowercase sequence has initial values v₀(x)=1 and v₁(x)=x and satisfies

\[
(n+1)^3v_{n+1}(x)=(2n+1)(n(n+1)+x)v_n(x)-n^3v_{n-1}(x),\qquad n\ge1.
\]

Run from the repository root with Python 3; there are no external dependencies:

```sh
python3 docs/reports/sun-apery-quintic-moment/check.py
```

`check.out` is the program's complete output. All evidence uses integers and
exact rational arithmetic. The program exits zero and ends with `ALL_CHECKS_OK`.

The fidelity checks reproduce Theorem 3.1's product formula in 64 exact
instances, Corollary 4.3's first linear congruence at every prime 3<p≤97,
and Theorem 4.4's linear congruence at 264 rational inputs for primes 5≤p≤43.

For m=2, the 168 inputs x=5/2+pt, with every prime 7<p≤97 and
t∈{0,−2,−1,1,2,1/2,−1/3,2/3}, all give

\[
\sum_{n=0}^{p-1}(2n+1)^5v_n(x)^2\equiv p/64\pmod{p^2},
\]

a nonzero residue under the source's hypotheses. Full exact sequences,
sums, and term residues are printed at p=11 and p=13 for both x=5/2 and
the positive integer x=(p+5)/2. At p=11 the residue is 55 modulo 121;
at p=13 it is 156 modulo 169.

Independent polynomial arithmetic checks the proposed degree and leading
coefficient of F modulo each prime, the odd-polynomial reflection rule in
four cases per prime, and the coefficient identity

\[
S(X)\equiv-44(X-5)(X-3)(X-2)(X-1)(X+4)(X^2-2X+2)\pmod{121}.
\]

Here S(X) is the quintic moment at p=11 computed by polynomial recurrence.
All 121 integer residues are also compared against the rational recurrence.
These are finite checks supporting a proposed general argument, not a proof
of the all-prime family and not a Lean settlement. The source's stronger
Conjecture 4.3 modulo p⁴ is outside this experiment's settlement claim.
