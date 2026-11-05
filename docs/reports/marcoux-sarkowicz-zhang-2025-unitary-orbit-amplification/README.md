# Exact unitary-orbit amplification gap

Target: Marcoux–Sarkowicz–Zhang, arXiv:2508.13834v1, Question 4.16,
noncoprime clause (85). Run with Python 3 and SymPy:

    python3 check.py

`output.txt` contains the output. Acceptance uses integers, fractions and exact
SymPy algebra. Decimal displays are directed rational enclosures, not numerical
optimization. `ALL_CHECKS_OK` covers the executed identities, bounds and data;
the analytic argument below supplies the universal unitary lower bound.

The program reproduces the source's Section 3.1 spectral-repetition formula
(404 cases), checks amplified padding identities and polar-compression examples,
and constructs explicit Gaussian-rational seeds of dimensions 3072 and 4096.
For the 4096 seed it gives delta_1 > 0.000005840674603 and
delta_3 < 0.000004567088664; their gap exceeds 1/3200000. Padding and independent
irreducibility perturbations give a contraction family in every dimension at
least 4096, with uniform gap at least 1/102425600000.

## Construction and lower-bound argument

Take h(t)=t I_3 and five Hermitian matrices E11,E22,E33,E12+E21,E13+E31.
The eleven functions are h and (cos(2 pi t)-1)K, sin(2 pi t)K for these K.
Evaluate at j/86, j=1,...,85, with one scalar endpoint at 0 on the left and
one at 1 on the right. Each tuple acts on dimension 256. Their h-traces differ
by -1. No common unitary can intertwine all eleven entries with error <=1/8:
for u_a=I+(exp(2 pi i t)-1)E_aa, the three determinant-one relative unitaries
have principal-log trace integers. Rank-one projection rotations make the
three integers equal. Their product remains within 3/4 of identity, so its
trace-log integer is their sum, a multiple of three. The path formed from the
two h-exponentials remains within pi/4 of identity and forces that same integer
to be -1. These paths never encounter the principal-log cut.

Encode the tuple with real diagonal labels spaced by 1000 and a common
imaginary star of weight 32, using q=12 or 16 labels. A single-matrix residual
epsilon<=1 gives off-diagonal unitary blocks <=epsilon/1000. Polar correction
and the star connections recover one common tuple unitary with error less
than (1.24525+0.00018 epsilon)epsilon < (5/4)epsilon. Hence the encoded base
orbit distance exceeds 1/10. Triple amplification combines three scalar
endpoints into one 3-dimensional endpoint and matches adjacent mesh points;
its explicit permutation witness costs at most pi/43 < 7/86.

Replace cos and sin by Taylor polynomials through degrees 24 and 23, at the
rational arguments 22j/301. Keep endpoint values exact. The integral identity
22/7-pi = integral_0^1 x^4(1-x)^4/(1+x^2) dx and the integer inequality
25! > 10000*7^25 give entry-block error less than eta=2063/630000.
The base distance is then greater than 1/10-2 eta >93/1000, and the amplified
witness is less than 7/86+2 eta <88/1000. Divide the q=12 matrices by 12000
and the q=16 matrices by 16000 to obtain contractions. All matrix formulas,
permutation coordinates, constants and directed enclosures are in `check.py`.

The finite construction is a proposed derivation supplied through a Nyx broker
reasoning call and checked by the codex-cli research seat. It is not copied
from published program source and has not been formalized in Lean. The source's
Proposition 4.14 and Marcoux–Zhang, JFA 280 (2021), 108778, Theorem 5.1 provide
prior existence results; they do not specify these finite seeds.

Sources: https://arxiv.org/abs/2508.13834v1 and
https://doi.org/10.1016/j.jfa.2020.108778 . The accepted 34-page JOT version
remains unread; the forthcoming-list date is a submission date.
