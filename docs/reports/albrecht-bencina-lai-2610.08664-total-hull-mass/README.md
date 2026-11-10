# Total Galois hull mass: exact checks

The proposed uniform formula counts actual subspaces C of F_(p^d)^n of dimension k whose intersection of all right Galois hulls has dimension h. It answers the counting direction in Albrecht, Benčina and Lai, *Algorithms for Sampling Self-Orthogonal and Totally Self-Orthogonal Codes in Odd Characteristic*, arXiv:2610.08664v1, §I-B. These finite checks do not prove the uniform formula.

Run with Python 3.12 or newer, using only the standard library:

```sh
python3 check.py
```

`output.txt` contains the complete output. All subspaces are enumerated once by their unique reduced-row-echelon matrices. Right-dual equations are stacked after applying the inverse Frobenius to make the equations linear in the codeword coefficients. The implementation independently checks the literal vector-membership definition for every code of length at most two and the first twelve RREF codes for each larger parameter pair. It also compares exact Frobenius-closure counts with subspace-lattice Möbius inversion and checks the hull-incidence moments for every parameter tuple.

The formula uses integer Gaussian coefficients, the classical odd-characteristic dot-product isotropic-space products, exact-closure Möbius inversion over F_p, and Gaussian incidence inversion over F_(p^d). All divisions use `fractions.Fraction` and are checked to have integral result; no floating-point arithmetic is used as evidence.

The exhaustive scope is p=3, d=1,2, 0≤n≤4 and p=5, d=1, 0≤n≤3, with every 0≤h≤k≤n. The 9,627 enumerated subspaces verify 90 mass parameter tuples. The largest Grassmannian has 7,462 two-dimensional subspaces of F_9^4; its observed hull-dimension counts are [7326,128,8]. All boundary cases in this scope, including n=k=h=0, pass.

For source fidelity, Lemma 45 is checked on all ordered vector pairs in F_q^2 for q=3,9,5 (81, 6561, 625 pairs). For every enumerated code the analogous subspace-level equivalence between total self-orthogonality and Euclidean isotropy of its Frobenius closure is checked. The lemma and the closure mechanism are credited to the source. The check program is independently written from the mathematical definitions and does not use the source's bundled programs.

The field F_9 is represented as F_3[u]/(u²+1); all arithmetic tables are exact. The field's Frobenius addition and multiplication laws and x^q=x are checked. This representation is a computational choice and is not a restriction on the proposed uniform theorem.
