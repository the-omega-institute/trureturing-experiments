# Seidel energy of chain graphs

Exact finite checks for Santanu Mandal, Ranjit Mehatari and Kinkar Chandra Das,
*On the spectrum and energy of Seidel matrix for chain graphs*,
arXiv:2205.00310v1, Conjecture 4.14.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py > output.txt
```

The check builds the graph from the binary creation string and uses
`S = J - I - 2A`. It reproduces Lemma 4.12 on all 66 members of
`Gamma_{n,2}` with `4 <= n <= 14`, verifying the characteristic polynomial
and its stated spectrum and energy consequence exactly. It verifies the
proposed signed-cycle inverse and signed two-step rotation for
`2 <= k <= 12`.

For all 511 positive run-length tuples of even length and total order
`2 <= n <= 10`, it verifies the proposed compression matrix and
characteristic-polynomial identities, and compares energies using certified
rational intervals isolating every real characteristic root. The comparison
is strict in 416 cases and equality holds in exactly 95 cases, precisely
those with every intermediate run of length one. In the 466 cases with
`k >= 2`, exact root-sign certificates also establish the finite instances
of the proposed energy identity.

`output.txt` contains the output of the successful run. All evidence uses
integers, exact radicals or rational intervals; no floating point is used.
These finite computations support the proposed universal proof and do not
prove the unbounded conjecture. The source definitions, binding Lean
statement and preregistration outline belong to the main-repository issue.
