# Prefix reversal generating-triple count

The target is Conjecture 1 in Saúl A. Blanco, Mikhail P. Golubyatnikov,
Elena V. Konstantinova, Natalia V. Maslova and Luka A. Nikiforov,
*Generating the symmetric group by three prefix reversals*,
[arXiv:2511.16959v3](https://arxiv.org/abs/2511.16959v3), Section 4.1;
[DMTCS 28:1, article 7](https://doi.org/10.46298/dmtcs.16975).

For each n, f(n) counts pairs 2 ≤ k < m < n whose three prefix reversals
of lengths n, m and k generate the full symmetric group. The conjectured
leading coefficients are 1/13, 1/9, 1/10 and 1/7 for the four residues
of n modulo four. The error term is O(n) on each residue class.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py > output.txt
```

The program uses integer permutation images and deterministic Schreier–Sims
subgroup orders, without pruning candidate pairs. An independent tuple-based
breadth-first enumeration verifies all subgroup orders through n=7.
Rational comparisons use fractions.Fraction; no floating-point values enter
the checks or output.

All 1,771 pairs for 3 ≤ n ≤ 24 pass. Fidelity checks reproduce the exact
eight n=7 triples listed in Section 4.2 and all 210 instances of Theorem 4
with 5 ≤ n ≤ 24. There are 20 independent BFS cross-checks and 275
involution checks. The parity-correct triangular exclusions and the two
forward-shift properties also hold throughout the computed range.

The finite consistency test fixes C=2 and N=6 before execution and checks
|f(n) − c_q n²| ≤ 2n for every 6 ≤ n ≤ 24. The maximum error divided by n
is 4/65, 20/117, 1/10 and 80/161 in residues 0, 1, 2 and 3 respectively.
These are finite measurements. Neither the unbounded O(n) error nor any
global classification of generating triples follows from them.

The scripts and output are independent verification material for a
preregistration draft. They contain no Lean proof and do not settle the
conjecture.
