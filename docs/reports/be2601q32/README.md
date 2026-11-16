# Core quandle ring idempotents: exact checks

The target is Bardakov–Elhamdadi, *Idempotents and Powers of Ideals in Quandle
Rings*, arXiv:2601.07057v3, Theorem 3.1 and Question 3.2.

Run `python3 docs/reports/be2601q32/check.py` from the repository root.
The standard-library program uses only integer and rational arithmetic, and
integer residue arithmetic for prime coefficient fields. The output is in
[output.txt](output.txt); success ends with `ALL_CHECKS_OK` and exit code 0.

The fidelity test checks Theorem 3.1 on every integer coefficient vector in
{-2,-1,0,1,2} on support window {-2,-1,0,1,2}. The extension checks all vectors
in the displayed coefficient windows on Z, and all supports of size at most
three in the displayed windows of lexicographically ordered Z² and the
noncommutative integral Heisenberg group. The latter uses coordinates (a,b,c)
with product (a,b,c)(d,e,f)=(a+d,b+e,c+f+ae).

Every square is computed by the source's right core product y x⁻¹ y, with
indices in the infinite group. No indices are reduced modulo a prime.
The checks compare the idempotent equation to the zero-or-basis classification,
verify the unique extreme contributing pairs and their nonzero coefficients,
and test the additive gap identity and Heisenberg group/order/quandle laws.

These finite checks provide fidelity and error detection. They do not prove
the universal assertion. The planned settlement uses a general finite-support
extremum argument, with no bound on the group or its support.

The program is independently authored; no third-party program or Lean source
was copied. It requires Python 3.12 or newer under the project execution rule.
