# Reflection complexity: Conjecture 28

Exact finite checks for Allouche–Campbell–Li–Shallit–Stipulanti,
*The reflection complexity of sequences over finite alphabets*,
arXiv:2406.09302v4, Definitions 1 and 3, Example 4 and Conjecture 28.
The target is the unrestricted bounded-difference statement, restated in
Section 6 of Dvořáková–Pelantová, DOI 10.1007/s00224-026-10278-7.

Run `python3 check.py`. Python 3's standard library suffices. All arithmetic
is integral; there is no random sampling or floating point. `output.txt`
is the complete deterministic output.

The fidelity check reproduces the fifteen source values r(0),...,r(14)
for Thue–Morse. To obtain its exact length-n language, choose k with
2^k >= n and enumerate factors of the four words mu^k(ab), where mu(0)=01,
mu(1)=10 and a,b range over the binary alphabet. All four pairs occur in
Thue–Morse. Every length-n window lies in at most two aligned mu^k blocks,
so this enumerates all factors rather than relying on prefix stabilization.

The remaining checks exhaust every presentation x = prefix.period^omega
in these ranges, allowing repetitions of the same infinite word:

| Alphabet size | Prefix length | Nonempty period length | Presentations |
| --- | --- | --- | --- |
| 1 | 0–3 | 1–5 | 20 |
| 2 | 0–3 | 1–5 | 930 |
| 3 | 0–2 | 1–3 | 507 |

For each presentation and n=0,...,17, all factors are obtained by starting
at positions 0,...,length(prefix)+length(period)-1. Later starts repeat
modulo the period. The two-sided marked mirror is z(0)=-1 and
z(i)=x(abs(i)-1) for i != 0; -1 is outside the nonnegative-letter alphabet.
The enumerated start interval contains every left-tail, marker-crossing,
and right-tail factor, using the same periodic reduction in both tails.

The program checks the exact mirror factor description and p_K(n)=u_n+n,
2r(n)=u_n+P_n, reversal invariance of graph endpoints, palindromic edges
being loops, connectivity after their deletion, and the resulting
P_n+P_(n+1)<=u_(n+1)-u_n+2 and
-1<=r(n+1)-r(n)<=u_(n+1)-u_n+1, for n=0,...,16.
There are 24,769 comparisons, including 2,160 negative reflection steps;
1,090 presentations lack reversal closure in the inspected range.

These results support the preregistered bridge lemmas on the declared
finite families. They do not prove the universal conjecture or the
Cassaigne theorem about factorial prolongable languages. No Lean proof
of the conjecture or its negation is included.
