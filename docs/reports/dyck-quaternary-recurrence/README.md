# OEIS A379464: exact recurrence checks

Run `python3 check.py` with Python 3 and SymPy. `output.txt` is the output of a successful run. The last line is `ALL_CHECKS_OK`.

The source is [OEIS A379464, revision 20](https://oeis.org/A379464/internal), with [revision history](https://oeis.org/history?seq=A379464). Emely Hanna Li Lobnig supplied the path definition; Thomas Scheuerle supplied the PARI binomial sum; R. J. Mathar supplied the conjectured eight-term recurrence. Source excerpts and initial values are attributed to these contributors and the OEIS Foundation under the [OEIS license agreements](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement) (CC BY-SA 4.0 for OEIS content).

The program uses integer height-state dynamic programming, independent exhaustive word enumeration through length 10, exhaustive forest codes through four internal vertices, exact `Fraction` factorial/binomial expressions, and symbolic polynomial identities over `ZZ`. It reproduces the source's initial 32 terms and its four paths of length four. It compares the original path count with Scheuerle's sum through length 512 and checks every conjectured recurrence instance from 7 through 512. It also checks the differential certificate's eight coefficient polynomials as identities and retains the inhomogeneous constant 6 and forcing polynomial `384*(1-z)^3`.

These are finite counting checks and exact operator-algebra checks. They provide no Lean proof of the all-index statement and do not certify the literature status. The proposed all-index bridge remains a proof obligation for the later formalization probe.
