# Doust–Weston circle negative type

Source: Ian Doust and Anthony Weston, *Squaring the circle: embedding S¹ in ℓ₁*, arXiv:2608.10422v1, Section 4.3, Question 4.9. Source URL: https://arxiv.org/html/2608.10422v1 .

Run with Python 3 and SymPy:

```sh
python3 docs/reports/doust-weston-circle-negative-type/check.py
```

`output.txt` is the complete successful output. All arithmetic uses integers, rational numbers, or symbolic exact algebra. No floating-point estimate supplies evidence.

The fidelity checks reproduce both odd-cycle identities in Section 4.3 for n = 3, 5, …, 23, and the exact C₅ threshold polynomial. The candidate checks establish the symbolic zero-sum decomposition of the powered-distance matrix with entries 1, A, B, its negative-square certificate at B = (A−1)², and its positive obstruction vector when B < (A−1)². They also check every distance for 81 rational geometric gap parameters, 1,875 endpoint vectors, and 45 obstruction cases. A sign-reversal negative control is rejected.

The planned universal realization is b = (2ᵖ−1)^(2/p), a = 2π/(b+3), X = {0,a,2a,3a}, for every real 1 < p < 2. The symbolic calculations apply to the full parameterized matrix; the finite checks are additional controls. The program checks the algebraic derivative rewrite but does not verify the real-analytic positivity argument, the real-power identities, or the universal claim in Lean. Those obligations belong to the proof stage.
