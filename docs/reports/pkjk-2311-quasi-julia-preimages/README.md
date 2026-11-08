# Quasi-Julia preimage checks

Source: Portik, Kálmán, Jex and Kiss, *Robustness of chaotic behavior in iterated quantum protocols*, arXiv:2311.13280v2, Sections II–IV; Physical Review A 109, 042410, DOI 10.1103/PhysRevA.109.042410.

Run `python3 check.py` with Python 3 and SymPy. The program uses integers, exact symbolic algebra and rational arithmetic. Its output is in `check.out`.

The exact checks derive equation (8) from the normalized entrywise square and the imperfect Hadamard gate, reproduce the purity and ideal two-cycle, and verify the product-coordinate identity for orders 2 through 6. They check 375 rational interior inputs, including the critical axis. Exhaustive checks over the 1, 4 and 29 topologies on respectively 1, 2 and 3 points verify 1,386 frontier identities for 179 continuous open self-maps and both common-boundary and boundary-union invariance for 790 cyclic pairs of subsets.

The Appendix’s printed inverse is inconsistent with equation (8): at ε = 0 it sends F(0, 0, 1/2) = (4/5, 0, 0) to (√(8/5), 0, 0), outside the physical ball. The check records this exact discrepancy and derives physical preimages directly from equation (8). This does not refute the boundary-preimage assertion.

The finite-resolution diagnostic samples six brackets between pure phases for coherent errors ε = 2 arctan(r), with r = −1/25, 0, 1/25. Each of the two physical inverse branches is within distance 1/10⁶ of points with both finite-time phase labels; squared forward residuals are below 1/10⁶⁰. Square roots and orbit steps are rounded to dyadic rationals with denominator 2¹²⁸; no binary floating-point operations are used. A negative control records a transverse arc that enters the maximally mixed basin, so it cannot serve as a two-phase boundary bracket.

These diagnostics certify finite algebra and finite trajectories. They do not certify infinite-time basin membership, exact frontier membership, density of inverse orbits, fractal dimension, a critical purity, or the universal preimage assertion. The proposed settlement uses a general topological argument; these finite checks are supporting evidence.
