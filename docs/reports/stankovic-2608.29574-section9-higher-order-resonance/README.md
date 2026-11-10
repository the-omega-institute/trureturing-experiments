# Higher-order resonance: exact finite checks

Source: Hranislav Stanković, *Spectral Rigidity of Commutators: Dynamics, Resonance, and Nilpotency*, arXiv:2608.29574v1 (30 August 2026), Section 9, first open direction. The relation is P(ad_A)(T)=0, where ad_A(T)=AT−TA.

The proposed complete characterization is: in dimension n≥1 a non-nilpotent solution exists if and only if the roots of a nonconstant complex polynomial P admit a positive zero-sum word of length at most n, with repetitions permitted. The universal argument uses primary-block support and a repeated vertex in a block walk, with a cyclic construction for the converse. This report contains finite checks only; no universal proof or Lean proof is claimed.

Run with Python ≥3.12 and SymPy 1.14.0:

```sh
python3 docs/reports/stankovic-2608.29574-section9-higher-order-resonance/check.py
```

All arithmetic uses integers, rationals and Gaussian rationals. The random seed is 260829574. The output is stored in `check.out`.

- Fidelity: the source's displayed 3×3 example in Section 6 has roots 2 and −1, satisfies ad_A²(T)−ad_A(T)−2T=0, and has det(T)=1 and det(ad_A(T))=2.
- Third order: roots 1, i, −1−i have shortest positive zero-sum word length 3. For the displayed cyclic pair, P(ad_A)(T)=0, T³=I and ad_A(T)³=(1−i)I.
- Small cases: 15 adversarial root multisets and 12 seeded random root multisets, each in dimensions 1 through 4. For 108 diagonal A and 12 integer conjugates of diagonal A, the program solves the exact polynomial-operator kernel. It computes the characteristic polynomial of a symbolic generic element of the entire kernel; all coefficients equal those of Xⁿ precisely when every element of that fixed solution space is nilpotent. It also checks the independent directed-support cycle criterion. A fixed A need not realize a root word: only the global existence statement allows A to vary.
- Jordan cases: 32 exact polynomial kernels for non-diagonalizable A, including repeated eigenvalues and repeated polynomial roots. Repeated roots affect kernel dimension: for the 2×2 nilpotent Jordan block, dim ker(ad_A)=2 and dim ker(ad_A²)=3. The proposed global root-word criterion ignores multiplicity.
- Construction: 45 exact padded cyclic witnesses in dimensions at most 4. When every root on the cycle is nonzero, the weighted-cycle identity is also checked for iterates 0 through 4.
- Primary-block algebra: 16 nilpotent Sylvester perturbations, for both block sizes 1 through 4, satisfy the predicted nilpotence bound and the explicit two-sided finite geometric inverse whenever P(a)≠0.
- Boundaries: n=1 reduces to P(0)=0; coprime two-root arithmetic thresholds p+q are checked for 1≤p,q,n≤4.

Across the 152 exact kernel spaces, 109 consist entirely of nilpotent matrices and 43 contain a non-nilpotent matrix. A nonzero coefficient polynomial of the generic characteristic polynomial gives existence over the infinite field ℂ. These are bounded test results, not a classification for arbitrary A or arbitrary dimension.

The program exits 0 and ends with `ALL_CHECKS_OK`. The Lean conventions and preregistration draft are separate artifacts; this report does not certify a kernel-checked settlement.
