# Exact certificates for the Ilyas–Khan–Loi–Serra-Capizzano GLT conjecture

The target is Conjecture 1 in arXiv:2505.03256v1, submitted 6 May 2025, with the candidate symbol in Definition 3, equation (4.5). The expected settlement is Refuted. These certificates support a proposed settlement; they are not a Lean proof.

Run with Python 3 and SymPy:

```sh
python3 -m venv /tmp/ilyas-glt-venv
/tmp/ilyas-glt-venv/bin/pip install sympy==1.14.0
/tmp/ilyas-glt-venv/bin/python check.py > check.output.txt
```

`check.py` uses exact integers, rational numbers, symbolic radicals and symbolic limits. It reproduces the source's explicit overlap matrix in Section 4.4.2, Case 1, Example 2, by checking its native geometric-mean formula and Riccati identity.

For the proposed family it checks positivity, determinants, Riccati equations, the singular-boundary limits for every positive t and nonnegative s, fixed noncommuting rank-three symbols, uniform bounds, norm-null perturbation certificates, and the two exact normalized-trace formulas. Direct-sum and Kronecker identities are independently checked at block multiplicities 1 through 8 for each shear. The two subsequential limits are 19/20 and 9/10, separated by 1/20. Output eigenvalues lie in [0,10/3], so a continuous compactly supported test equal to the identity on [0,4] gives the same nonconvergent averages.

The two-by-two mean formula and discontinuity at singular boundaries are classical: Ando–Li–Mathias, *Geometric Means* (2004), Proposition 2.1 and Corollary 2.2; Bhatia, *Positive Definite Matrices* (2007), Exercises 4.1.6 and 6.2.4. The proposed addition is their GLT lifting with fixed noncommuting symbols and alternating boundary paths. The GLT generator/distribution implications and the final refutation remain obligations for the later Lean probe.

`check.output.txt` is the complete successful output; its final line is `ALL_CHECKS_OK`.
