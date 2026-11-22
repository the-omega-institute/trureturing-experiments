# Mąkowski–Schinzel inequality for integers with smooth prime predecessors

The Mąkowski–Schinzel conjecture (Colloq. Math. 13 (1964/65), 95–99) asserts σ(φ(n)) ≥ n/2 for every positive integer n. C. Pomerance, *A bound for the Mąkowski–Schinzel conjecture concerning the arithmetic functions σ and φ* (preprint dated 2026-11-13, https://math.dartmouth.edu/~carlp/msconj2.pdf), proves the cases n ≤ 10^60 and P(n) < 10^6. On p. 7 he remarks that the case P(φ(n)) ≤ 5 follows "with a little extra work" and that "this type of result can be further extended"; no proof of that case is given there.

`smooth_totient_certificate.py B X` checks the following statement exactly: **if every prime factor q of n has B-smooth q − 1, then σ(φ(n)) ≥ n/2, and σ(φ(n)) ≥ n when n is odd.** This hypothesis is weaker than P(φ(n)) ≤ B; for example 101² satisfies it for B = 97, while φ(101²) = 10100 has the prime factor 101.

    python3 experiments/makowski-schinzel-smooth-predecessor/smooth_totient_certificate.py 97 1e6   # needs sympy and gmpy2

## Reduction and certificate

By α(rad n) ≤ α(n) and α(2m) = α(m)/2 for odd m, where α(n) = σ(φ(n))/n, it suffices that every odd squarefree m with the hypothesis satisfies α(m) ≥ 1. Writing s_p(e) = σ(p^e)/p^e and e_p = v_p(φ(m)), one has α(m) = s_2(A) ∏_{odd p} s_p(e_p) ∏_{q|m} (1 − 1/q). The Fermat part s_2(A) ∏_{Fermat q|m} (1 − 1/q) is at least 1 (Cohen, Colloq. Math. 74 (1997) 1–8). So it suffices to prove ∏_{q|m, q not Fermat} q/(q − 1) ≤ ∏_{odd p ≤ B} s_p(e_p).

The certificate follows the grouping method of Cohen (1997, Theorem 3).
- **Candidates.** Every q = m + 1 with m ≤ X even, B-smooth and not a power of two is a candidate unless gmpy2's Miller–Rabin test proves it composite. A composite verdict is a proof, and survivors need no primality proof.
- **Labels.** Each candidate gets a label p, an odd prime dividing q − 1. Labels are assigned greedily (forced labels first, then local search). A prime q > X is labelled by the largest odd prime factor of q − 1.
- **Prefix tests.** For each label p, with members c_1 < c_2 < … and L_j = ∏_{i ≤ j} (1 − 1/c_i), the program checks s_p(j) L_j ≥ 1 and s_p(j + 1) L_j (1 − T_p) ≥ 1 for all j ≥ 0.
- **Tail bound.** T_p = (1/(2p)) (∏_{r ≤ p} r/(r − 1) − Σ_{m′ ≤ X/(2p), m′ p-smooth} 1/m′) bounds the reciprocal mass of labelled primes above X.
- **Why it suffices.** Every member contributes at least one factor p to φ(m), so e_p is at least the number of members. The j smallest chosen members are at least c_1, …, c_j. Each group therefore contributes a factor at least 1.
- **Exact arithmetic.** L_j is a floor-rounded fixed-point lower bound with scale 10^60, and the finite sum in T_p is floor-rounded, so every comparison is exact.

## Results

| B | X | candidates | result | minimal ratio | time |
|---:|---:|---:|---|---:|---:|
| 31 | 10^5 | 1 159 | pass | 1.00273 | 0.4 s |
| 47 | 10^5 | 1 686 | pass | 1.00011 | 0.3 s |
| 61 | 10^5 | 2 027 | fail | 0.99922 | 0.3 s |
| 97 | 10^6 | 9 992 | pass | 1.00004 | 0.9 s |
| 127 | 10^7 | 45 770 | pass | 1.00006 | 13 s |
| 157 | 10^7 | 56 651 | pass | 1.00007 | 18 s |
| 199 | 10^7 | 71 603 | fail | 0.99986 | 26 s |
| 199 | 10^8 | 274 998 | pass | 1.00019 | 696 s |

A failure only means that this certificate does not close at the given X; with a larger X the tail bounds shrink and the same B passes. The passing row B = 199 proves the statement for every n whose prime factors q all have 199-smooth q − 1. Times are wall-clock seconds of one run each, measured with `/usr/bin/time` on the author's laptop (Apple silicon, Python 3.15, gmpy2 2.3.2, sympy 1.14.0).
