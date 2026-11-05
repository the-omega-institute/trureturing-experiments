# Periodic convergents and monochromatic progression lengths

Source: Gandhar Joshi, *Automatic sequences from algebraic, dynamical, and
computational perspectives*, Open University PhD thesis, Conjecture 6.1.2,
DOI [10.21954/ou.ro.00106631](https://doi.org/10.21954/ou.ro.00106631).

Run `python3 docs/reports/joshi-periodic-convergent-multiples/check.py` from the
repository root. Python 3.8 or later suffices; no third-party packages are used.
The stored output was produced with Python 3.14.8. Success ends with
`ALL_CHECKS_OK`; any failed check raises an exception and exits nonzero.

The source word uses phase {(i+1) alpha}, with the first letter on
(0,1-alpha) and the second on (1-alpha,1). Convergents use q_0=1.
The script represents every real quantity as a+b sqrt(D), with rational a,b.
It compares signs by rational squaring and certifies every computed floor
by k <= x < k+1. No floating-point arithmetic enters the evidence.

The circle computation successively intersects the intervals of starting
phases permitting each letter at each progression position. It stops only
when that set is empty, giving an upper bound for every start, and records
nonempty open intervals at the preceding length. Irrational-orbit density
connects these open intervals to actual word starts. That density statement
is mathematical background, not a theorem proved by the script. A finite
prefix is never used to assert a global upper bound.

Results in `check.out`:

- Fidelity: for alpha=sqrt(5)-2, it reproduces the thesis's explicit
  A(17)=59, A(34)=30, A(51)=20 and A(94316)=61; it checks q_0,...,q_3
  = 1,4,17,72 and five constant-period convergent-error identities.
- Family: for a=3,...,8 and n=1,...,8, it checks 48 exact identities
  q_n alpha-p_n=(-1)^n a^ceil(n/2) alpha^(n+1), where
  alpha=(sqrt(a^2+4a)-a)/(2a)=[0;overline(a,1)]. The corrected ceiling is
  strictly below the conjectured ceiling in all 48 cases.
- Independent maxima: for the same six a values and n=1,...,4, the circle
  computation gives 24 global phase maxima, agreeing with the corrected
  formula. It also checks 54 multiplier cases y=2,...,4 satisfying
  y a^ceil(n/2) alpha^(n+1) <= alpha.
- Witness: a=3,n=1 gives A(3)=4 versus the conjectured 11. Source positions
  13,16,19,22 all carry the first letter.

These are finite exact checks. They corroborate the proposed infinite-family
refutation; they do not prove that family, the conjecture's negation, or any
Lean theorem. The preregistration's universal argument remains a proof
obligation for the later formalization stage.
