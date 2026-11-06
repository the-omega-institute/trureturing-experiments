# Measurement-induced mixed-state basin

`check.py` checks the proposed rational Lyapunov route for the all-interior-state
convergence assertion in Kálmán, Kiss and Jex, *Sensitivity to initial noise in
measurement-induced nonlinear quantum dynamics*, arXiv:1808.08151v1,
J. Russ. Laser Res. 39, 382–388, DOI 10.1007/s10946-018-9732-2.

Run with Python 3 and SymPy:

    python3 check.py > check.out

The fidelity check derives the source's Eq. (10) from normalized entrywise
squaring followed by the source's unitary, and reproduces its purity formula.
The algebra checks establish zero symbolic residuals for the purity defect,
the defect-normalized recurrence and the Lyapunov identity. Nonnegative
polynomial coefficients certify the coercive estimates used in the outline.
Finite rational grids check the same identities and bounds, including zero
coordinates and the maximally mixed state; exact boundary examples check the
source's fixed points and transient purification.

All arithmetic used as evidence is exact: integers, Fraction and symbolic
rational functions. The program never uses floating-point convergence as
evidence. `check.out` contains the program's output and ends in `ALL_CHECKS_OK`.
These checks support a mathematical proof outline; they are not a Lean proof
of convergence or a numerical exclusion of exceptional orbits. The proposed
settlement uses the general Lyapunov and summability argument, not finite
positive enumeration.

The program is original code under Apache-2.0. Source mathematics is attributed
above; no third-party program or Lean proof has been copied.
