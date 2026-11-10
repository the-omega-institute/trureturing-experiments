# Three quintuple products modulo 31

The source is Daniels–Huber–McLaughlin–Ye, *Vanishing Coefficients in Products of Quintuple Products*, arXiv:2606.06701v1, Section 5. The companion arXiv:2603.04666v1, Section 6, gives the related modulo-13 and modulo-19 experiments.

For each required specialization, the program expands

\[
Q(q^a,q^m)=(q^a,q^{m-a},q^m;q^m)_\infty
(q^{m+2a},q^{m-2a};q^{2m})_\infty.
\]

All five initial exponents are positive. Omitted factors have degree greater than the requested cutoff and cannot affect any reported coefficient. A separate implementation uses the bilateral quintuple identity, with a bound obtained by completing the square and taking an integer square root. Both implementations use Python integers only. Matrix checks use exact rational numbers from `fractions.Fraction`.

Run `python3 check.py`. The output is stored in `output.txt`.

- Fidelity: reproduce Example 4.1(1) at b=1, including both modulo-37 zero progressions through t=20, and compare every coefficient with the source's quintuple identity.
- Target: compare all 639 coefficients through degree 638 and check all 63 asserted zeros for t=0,…,20.
- Related checks: the companion's modulo-13 and modulo-19 examples, the composite modulus 49, and 16 admissible pairs with 1≤A<B≤8 through t=5.
- Mechanism: test the reflection and signed cyclic maps on all nonempty target shells through degree 638. Their matrices are orthogonal and square to the identity; each map preserves its shell and reverses its weight.

These are finite checks. They establish neither the universally quantified vanishing assertion nor a Lean proof. The planned settlement requires the full formal product-to-lattice bridge and a cancellation argument valid at every degree. Code in this directory is independently authored and licensed under Apache-2.0.
