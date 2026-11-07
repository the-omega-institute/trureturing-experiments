# Membrane signature dimension checks

Source: Carlos Améndola, Gabriel Riffo and Leonard Schmitz, *Universal Varieties
for Signatures of Membranes*, [arXiv:2610.16724v1](https://arxiv.org/abs/2610.16724v1),
Conjecture 5.10, equation (26). The target is the simultaneous pair of
second-level signature matrices.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The recorded output is
`check.out`; a successful run exits zero and ends with `ALL_CHECKS_OK`.

The fidelity checks integrate the moment-curve formula exactly, reproduce the
displayed polynomial in the example labelled `ex:momentMembrane`, and compare
the Jacobian implementation with symbolic differentiation. The rank checks use
three reproducible integer parameter matrices for every triple with
1 ≤ d ≤ 4 and 1 ≤ m,n ≤ 5, retaining both signature matrices. All 300 ranks are
computed over the rational field. For each triple the maximum sampled rank
equals the conjectured dimension. The relevant source dimension-table entries
are also checked. The 25 cases d = 0 and 1 ≤ m,n ≤ 5 are zero-dimensional.

With B_r having diagonal 1, upper entries 2 and lower entries 0, the checked map
is A ↦ (A(B_m ⊗ B_n)Aᵀ, A(B_m ⊗ B_nᵀ)Aᵀ). It scales each integral coordinate
by 4, which preserves Jacobian rank. No floating point or finite-field ranks
are used.

These sampled ranks give lower bounds on generic rank. They do not certify the
generic upper bound or prove the universal conjecture. The uniform-grid
definition divides by m and n, so zero grid sizes are outside these checks;
s_d(0) = 0 is checked where it occurs in the dimension formula.
