# Three-point lattice sumset diameter

Target: Melvyn B. Nathanson, *A problem on sumset sizes of sets of lattice points*, arXiv:2607.20571v1, Problem 3. Intervals and integer cubes include both endpoints. The proposed inequality is `N_n(h,k) ≤ N(h,k)^(1/n)`.

Run with Python 3, without third-party packages:

```sh
python3 docs/reports/nathanson-2607.20571-problem3-three-point-diameter-family/check.py
```

The stored output is `check.out`. The program uses only integer arithmetic and exits successfully with final line `ALL_CHECKS_OK`.

The fidelity check reproduces finite instances of the source's Theorem 1: its mixed-radix map preserves the actual h-fold sumsets and is injective on them. It tests all triples in `[0,2]^2` and `[0,1]^3`, for `1 ≤ h ≤ 6` (840 cases). A separate check compares 270 coprime integer triples against the published three-point counting lemma in Nathanson's Theorem 9, arXiv:2411.02365v5. That lemma is prior art.

The exhaustive check computes actual three-point sumsets for every translation class fitting the tested cubes. It examines all 56 triples of `[0,7]`, all 18,424 triples of `[0,6]^2`, and all 6,666,891 triples of `[0,6]^3`; coordinatewise translation leaves 21, 2,604, and 340,299 representatives respectively. Inclusion-exclusion independently checks the representative counts. Mixed-radix encoding replaces vector sums by exact integers; the enumeration uses neither the affine-rank lemma nor the coprime counting formula. Taking cumulative spectra by cube side also covers every smaller nonnegative side. The target full integer spectrum used for comparison is the published Theorem 9 spectrum.

For `h = 2,...,6`, the minima are `N(h,3) = h+1` and `N_2(h,3) = N_3(h,3) = h`. Each tested inequality fails by the exact integer comparison `h^n > h+1`. The extra `h=1` control gives minima 2 in one dimension and 1 in dimensions two and three; it is outside Problem 3's inherited `h ≥ 2` domain.

These are exhaustive finite checks. They do not establish the infinite family for all `h ≥ 2` and `n ≥ 2`. No Lean proof is included.

The candidate is withdrawn from preregistration because the universal inequality is already refuted in Henry Shin, *Iterated-sumset spectra: The complete exponent law and its rank geometry*, [arXiv:2609.01690v1](https://arxiv.org/html/2609.01690v1), Corollary 10.8, equation (199). That corollary explicitly identifies Nathanson’s Problem 3 and proves its failure for `k=4` and every `h,n ≥ 2`, with the same inclusive endpoint convention. The data here concerns `k=3` and is retained as an independent finite check, not a new settlement of the already-refuted universal assertion.
