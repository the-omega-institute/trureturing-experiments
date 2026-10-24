# Generalized Glasby–Paseman eventual log-concavity

Source: Seok Hyun Byun and Svetlana Poznanović, *Unimodality and log-concavity of Generalized Glasby-Paseman sequences*, arXiv:2604.14639v1, equation (1.4), Theorem 1.2(a),(b), and Table 2. The target is Conjecture 1.1(b).

Run `python3 check.py`; `output.txt` is its complete standard output. The script uses Python's standard library, arbitrary-precision integers and `fractions.Fraction`. It uses no floating-point evidence.

For a reduced weight `a=u/v`, both literal sums at index `r` are multiplied by `v^(r*l)`. Their ratio is unchanged; integer cross-products decide every inequality. An independent literal Fraction evaluation checks this representation on 1,080 entries.

Source fidelity reproduces the central-binomial denominator identity through row 100, Theorem 1.2(a),(b) through `m=100`, and all 200 cells of Table 2 from 20,000 complete finite sequences. These known cases earn no new resolution credit.

The proposed proof route is checked at 2,232 triples: `a` in `{1/7,1/2,1,3/2,3,10}`, `l=1..12`, `r=1..31`. The program verifies the all-index denominator inequality

`D[r-1]*D[r+1]/D[r]^2 >= ((r+1)/(r+2))^l`

and the strict fixed-index limiting margin

`((r+1)/r)^l * D[r-1]*D[r+1]/D[r]^2 > 1`.

A separate grid scans `m=1..100` for the same weights and powers, checks numerator-prefix log-concavity, and reports finite exceptions. A largest observed exception equal to 100 is censored by the scan ceiling. In particular `a=1/7,l=10,m=100,r=18` is an exact failure of global log-concavity, also recomputed using literal Fractions. No universal finite threshold such as 80 is claimed.

These bounded checks establish neither eventual log-concavity for all parameters nor the analytic expansion used in the ordinary proof outline. The universal settlement and its Lean proof remain future work.
