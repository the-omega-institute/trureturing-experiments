# Exact checks for Fillman–Lim Question 2

Source: Jake Fillman and Yi-Sheng Lim, *On the identification of certain classes of gap-rich sets*, arXiv:2609.21029v1, Section 6, Question 2. This report checks a proposed negative answer to all three universal closure assertions. It contains no Lean proof.

Run with Python 3.12 or newer and SymPy 1.14:

```sh
python3 check.py
```

`check.out` is the complete successful output. All arithmetic is integer, rational, or symbolic polynomial algebra. The final line is `ALL_CHECKS_OK`.

The fidelity check reproduces Appendix A, Lemma A.1, for lengths 1–6: the transfer trace has leading term ∏(E−xᵢ), with residual degree at most n−1 in the potential coordinates. It also checks determinant one.

The candidate uses qₙ=(2j+1)/2^(m+1) for n=2^m+j, rₙ=3^(−n), Iₙ=[10+rₙ,10+11rₙ/10], and P=([0,1]×{10})∪⋃ₙ({qₙ}×Iₙ), where n≥1. For any h>2 put Lₕ(q,t)=(t,q,2q−h−t), Uₕ(q,t)=(t,q,2q−h−t,q,q), and Vₕ(q,t)=(q,t,q,2q+h−t,q). Their images are nonempty compact perfect inputs.

For S(a)=[[a,−1],[1,0]], H=S(h), and J=S(0), the symbolic checks establish:

- Each of Lₕ, Uₕ, Vₕ has relative trace 2+(E−q)²(t−s)² when varying t within a component.
- At E=q the respective transfers are −H, H, and −H⁻¹.
- The coordinatewise sum Lₕ(q,t)+L₋ₕ(q,s) has transfer −J at E=2q.
- The concatenated block (Uₕ(q,t),Vₕ(q,s)) has transfer −I at E=q.

The exact interval inequalities isolate component n with radius rₙ/4. The tests cross-check these bounds for n=1–12, summand indices 0–16 (index 0 denotes the baseline), and the entire tails are covered by the displayed rational inequalities in the program. Short-word tests check invariance under fixed prefixes and suffixes and periodic alignment, including the impossibility of odd-length perturbations of a separated alternating union word. The semigroup check covers all 44 nonhyperbolic integer determinant-one matrices with entries in [−2,2], including all 592 ordered pairs with hyperbolic relative matrix.

Mathematical scope: the program does not decide infinite gap-richness or certify compactness/perfectness. The full argument additionally uses elliptic positive-power recurrence; the trace-±2 identity tr((I+N)^p(I+M)^q)=2+pq·tr(NM) for square-zero N,M; density of the dyadic midpoint tail; and the source's closure-invariance proposition. Together these supply input gap-richness at every energy and dense sets of obstructed energies for each output. A discrete exceptional set cannot contain either dense set. These analytic and topological steps remain obligations for the post-preregistration Lean probe.

Literature boundary: the source and prior spectral papers supply definitions and semigroup ingredients. These exact checks are an independent check of the proposed construction, not a claim of a kernel-verified settlement or an exhaustive prior-art search.
