# Pemmasani zero facet residual converse

Source: A. Pemmasani, *Occupation counts of rarely switching Markov chains: stratified local laws, support selection and the boundary layer*, arXiv:2610.09258v2, Proposition 6.2 and the open question in Section 9. The v2 version date is 26 October 2026.

Run `python3 check.py` with Python 3 and SymPy. The adjacent output records the Python and SymPy versions. All arithmetic used as evidence is rational or symbolic exact algebra.

The fidelity check reproduces Proposition 6.2(e): the exact value of omega, the third weighted-rate derivative, and the cubic residual coefficient. An independent expansion of the square-root variational objective gives the same coefficient. The Laplacian pseudoinverse satisfies both reflexive Penrose identities on this example.

The finite screen covers all 256 three-state generators with internal facet rates in {1,2} and the four entry/exit rates in {0,1,2}, with each entry/exit vector nonzero. For every generator it uses optimizer ratios in {1/2,1,2}; the corresponding stationary facet profiles are defined exactly in the code. There are 768 profiles: 156 have equal laws, 612 unequal laws, and 576 include zero entry or exit rates. Nonzero cubic coefficients detect all 612 mismatches. For the 156 equal cases, exact nonnegative-square identities show the fixed optimizer gives the affine value for every nonnegative square-root occupation parameter. Nine singleton-facet examples are also checked.

These finite checks support the proposed general analytic argument. They do not enumerate all real profiles, all larger chains, or certify the universal converse. No Lean proof of either the claim or its negation is included.
