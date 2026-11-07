# Centered anisotropic fidelity composition

`check.py` uses exact rational arithmetic and SymPy algebra to reproduce Example 1
of Jesús Lacalle, Luis Miguel Pozo Coronado and Rafael Martín-Cuevas, “Fidelity of
the sum of independent quantum computing errors”, Quantum Information Processing
25, 38 (2026), DOI [10.1007/s11128-025-05037-5](https://doi.org/10.1007/s11128-025-05037-5).
It obtains squared fidelities 4/9 and 13/28 and their composed value 127/252.

For normalized surface measure on S³, it checks the proposed density family

\[
f_{k,a}(q)=\frac{(1+q_0)^k[1+a q_2(q_0-k/(k+5))]}{Z_k},\qquad
Z_k=\frac{2(2k+1)!!}{(k+2)!},\quad k\ge2,\quad 0<a<1/20.
\]

The first moment is k/(k+3) times e₀. Its positive scalar makes e₀ the unique
extrinsic least-squares center. The nonzero mixed second moment is

\[
m_{k,a}=\frac{5a(2k+3)(2k+5)}{(k+3)(k+4)(k+5)^2(k+6)}.
\]

Under independent sequential SU(2) errors with first column
(q₀+iq₁, q₂+iq₃), its squared fidelity is
p=(k²+3k+6)/((k+3)(k+4)); composition gives p²+(1−p)²−2m².
The scalar identity in Conjecture 1 would give p²+(1−p)².
At k=2 and a=1/40, the exact defect is −9/10035200.

The program independently integrates polynomial densities using the normalized
sphere monomial formula for 39 rational cases (k=2 through 14, and a=1/100,
1/40, 1/21). Symbolic Beta-moment identities verify the parameter formulas,
composition, sufficient interior bounds and limiting expressions. There are no
floating point inputs or approximate comparisons.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/anisotropic-fidelity-composition/check.py
```

The stored output is `check.out`, ending with `ALL_CHECKS_OK`.

The calculation verifies this explicit SU(2) interpretation. The journal defines
the sum of errors only for isotropic variables; identifying the operation for
two anisotropic laws with the proposed SU(2) composition requires a separate
source-faithfulness decision. No Lean proof or completed settlement is claimed.
The SU(2) parametrization is also given in Saha Roy and Lacalle,
[arXiv:2606.04817v1](https://arxiv.org/abs/2606.04817v1), equation (6).
