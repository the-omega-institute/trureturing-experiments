# Quaternionic Deutsch equality: exact finite checks

Source: K. Mahesh Krishna, *Quaternionic Buzano-Dragomir Inequality and
Deutsch Entropic Uncertainty Principle*, arXiv:2610.35881v1, Theorem 2.1,
Theorem 2.9 equation (11), and Question 2.12.

Run `python3 check.py` with Python 3. No third-party packages are needed.
The stored output is `output.txt`. Successful completion prints
`ALL_CHECKS_OK` and exits zero. All arithmetic uses integers and
`fractions.Fraction`; no floating point or numerical logarithms occur.

The program implements Hamilton multiplication, quaternionic conjugation,
and the right Hilbert inner product sum(conjugate(a_i) b_i). It checks the
source's quaternionic Buzano inequality on every generated triple in pools
of size 5, 7, and 9 in dimensions 1, 2, and 3 (1,197 triples, including
repetitions). Radical inequalities are decided by sign-aware exact squaring.

For each dimension 1, 2, 3 it generates eight frame parameter combinations:
the coordinate basis or the rational rotation with entries 3/5 and 4/5,
quaternionic row and column phases, scale 1 or 1/2, and one or two copies.
Each frame atom has weight 1/(copies*scale^2). The operator identity
sum(weight*a*a^*) = Id is checked entry by entry, hence Parseval is verified
for every state, rather than only the sampled states. Some parameter
combinations generate the same frame. All states in {1,i,j,1+i}^d are used.

The source's entropy bounds and the proposed lower-equality classification
are checked on all 5,376 resulting frame-pair/state tuples. The program
decides the entropy comparisons without approximating log: with N clearing
the rational probability denominators, it forms exactly
Z = product(p^(N*weight*p))*product(q^(N*weight*q)) = exp(-N*(S_tau+S_omega)).
The Deutsch threshold M^N is evaluated in the rational quadratic algebra
with M=(1+sqrt(D))^2/4 and D the exact maximum squared cross modulus.
Comparing an element a+b*sqrt(D) to a rational uses exact signs and squaring.
Equality occurs on 64 generated tuples, all in dimension one with unit
frame vectors at positive-mass atoms. No equality occurs in dimensions 2
and 3. The upper bound is checked by comparing Z with (total_mass)^(-N).

Two further weighted-atom examples distinguish the source's everywhere
nonzero coefficient domain from almost-everywhere entropy properties:
a zero frame vector at a null atom destroys the domain, whereas a nonzero
subunit vector at a null atom preserves equality. Exact coefficient
multiplication checks the saturation polynomial factorization and its two
roots, 1 and -1/2.

These checks reproduce source results on the specified finite families.
They do not prove the classification for arbitrary quaternionic Hilbert
spaces or continuous frames. The planned settlement uses a general
entropy-rigidity and Parseval argument; finite enumeration is supporting
evidence only.
