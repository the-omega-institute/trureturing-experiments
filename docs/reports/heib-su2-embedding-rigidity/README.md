# Compact su(2) embedding rigidity: exact checks

Source: Tim Heib, *Finite-dimensional Lie subalgebras of the universal enveloping algebra of sl(2,C)*, arXiv:2610.04001v1 (2 October 2026), Conjecture 33; Example 26 and Lemma 31.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored output is `check.out`. Arithmetic uses integers, the exact Gaussian-rational domain QQ(i), and symbolic rational matrices; no floating point is used as evidence.

The program reproduces Example 26's nonlinear complex embedding using PBW reduction and shows that it fails the fixed compact dagger. It reproduces Lemma 31's polynomial spin representations for n = 0 through 8, including the factorial Gram matrix, cyclic brackets and positive Casimir n(n+2)/4. PBW first-reduction choices and dagger involutivity are cross-checked on all 1,093 words of length at most six.

The expected positive settlement is checked on all 531,441 affine triples whose generator coefficients lie in {-1,0,1} and whose imaginary constant coefficients lie in {-i,0,i}. Bracket constraints leave exactly 25 maps: the zero map and 24 embeddings, all with the canonical image. All 96 partitions of representation dimensions 1 through 9 satisfy the Casimir bound. For all 125 central-polynomial rescalings with polynomial degree at most two and coefficients -2 through 2, exact PBW bracket defects leave only zero and the canonical embedding.

These finite checks do not prove Conjecture 33, the PBW theorem, the center calculation, the abelian-centralizer theorem, or the unbounded polynomial-growth argument. They supply source fidelity, bounded checks and negative controls for a separate general proof. No Lean proof is part of this report.
