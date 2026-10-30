# MPS cumulant coefficient nonvanishing

`check.py` checks the proposed refutation of Remark 3.31 of Paul Mansanarez,
Guillaume Poly and Yvik Swan, *Derivatives of entropy and the MMSE conjecture*,
arXiv:2311.04831v2. It uses exact SymPy integer and rational algebra throughout.
Run `python3 check.py` with Python 3 and SymPy 1.14.0; `output.txt` is its output.

The source recurrence independently reproduces the full displayed polynomials
R₂ through R₆ and the 42-term R₉ example, including its coefficient a₈.
The program checks cumulants by both the moment recurrence and formal logarithm,
the quartic cancellation at a unique root in (1/100,1/20), and the nonvanishing
of the second, fourth, sixth and eighth cumulants at that root. It also checks
the algebra for Gaussian addition, scaling, the posterior-variance Hessian bound,
and the strict variation of the normalized fourth cumulant.

Remark 1.5 prints the mixed term in a₆ as −210 T₃T₆. The defining recurrence
gives −210 T₃T₅. This discrepancy does not affect the symmetric specialization,
where the mixed term vanishes.

The exact algebra supplies a three-point counterexample and certificates for the
outlined continuum of Gaussian-smoothed examples. Continuity and the intermediate
value theorem, zero integral of a nonnegative function, differentiation under the
finite Gaussian sum, and the classical finite-MGF moment-determinacy criterion are
mathematical steps in the accompanying preregistration, rather than conclusions
certified by SymPy. No Lean proof is included. This experiment addresses only the
universal nonvanishing assertion, not MMSE injectivity or entropy monotonicity.
