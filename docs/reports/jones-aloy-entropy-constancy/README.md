# Jones–Aloy–Higgins–Müller entropy constancy

Source: [Certified randomness from quantum speed limits](https://arxiv.org/abs/2506.14526v3), Appendix E(ii), conjecture after Lemma 3; PRX Quantum 7, 010349, DOI [10.1103/bnrq-sp5w](https://doi.org/10.1103/bnrq-sp5w). The inputs have probability one half; time separation and Planck's constant are one. Branch energy labels are nonnegative and have exact prescribed mean.

Run with Python 3 and SymPy 1.14:

```sh
python3 docs/reports/jones-aloy-entropy-constancy/check.py
```

The stored output is `check.out`. The program uses integers, rational arithmetic and exact symbolic algebra. It encloses logarithms by the series ln(x) = 2 sum z^(2k+1)/(2k+1), z=(x−1)/(x+1), after reducing x to [1,2]. The omitted positive tail is bounded by 2z^(2N+1)/((2N+1)(1−z²)). No floating point or numerical optimizer supplies evidence.

Fidelity checks reproduce the reflection invariance and midpoint identity in source Lemmas 2–3 for arbitrary symbolic correlations. Family checks verify the three-branch ensemble, exact mean energy, derivatives of asin(sqrt(p)), the entropy-gap algebra and Pauli realizations. A separate rational logarithm enclosure proves the entropy gap at s=24/25 is greater than 4330451/100000000 nats.

The analytic interpretation is as follows. Let b be binary entropy in bits, 24/25≤s<1, q=s², c=sqrt(1−s²), E=asin(s), B=(1,1−2q), and M=(q,−q). For f(p)=asin(sqrt(p)), its tangent at q lies strictly below f except at q: f is strictly convex above one half; below one half its difference from the tangent is concave with positive endpoints. At zero, the difference is s/(2c)−E≥12/7−π/2>0. An ensemble with mean B has every active first coordinate equal to one. The exact mean energy then forces every active second coordinate to equal 1−2q. Consequently H*(B,E)=b(q)/2.

The points (s,−s), (1,1), (−1,−1) with weights s, (1−s)/2, (1−s)/2, all labeled by energy E, have mean M and entropy s b((1+s)/2). Set a=(1−s)/2 and L=log₂(1/a)>11/2. The inequalities −(1−a)ln(1−a)<a and ln(2)>2/3 give

b(q)/2 − s b((1+s)/2) > a[L−11/2+7a] > 7a² > 0.

Thus the proposed counterexamples form an interval of energy values, with equal coordinate differences at B and M. Compactness of the capped-energy convex roof gives lower semicontinuity. It extends the strict inequality to B_t=(1−t)B+tM for sufficiently small positive t. The check verifies that the observable −tsc σ_x+(1−tc²)σ_z realizes B_t on states (1,0) and (c,s), and both effects have trace one and determinant c²t(2−t)/4>0 for 0<t<1.

The program checks algebraic and inequality certificates used by this argument. It does not implement a proof checker for boundary rigidity, compactness or lower semicontinuity; no Lean proof of the conjecture or its negation is supplied. The work is an exact reproducible check accompanying a preregistration draft.
