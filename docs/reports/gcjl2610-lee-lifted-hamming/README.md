# Lifted-Hamming toric-code checks

Exact-arithmetic evidence for the proposed all-exponents settlement of Example 7.10 in Guo–Han–Chakraborty–Jain–Liu–Albert–Lucas, *Quantum codes in the Lee metric*, arXiv:2610.04834v1, Section 7.3. This directory contains a checker and its output; it contains no Lean proof.

Run `python3 check.py`. Python 3.12 or later and SymPy are required; the recorded run used Python 3.14.8 and SymPy 1.14.0. The last line must be `ALL_CHECKS_OK` and the exit code must be zero. All calculations use integers, symbolic rational expressions with exact cancellation, or explicit finite quotient rings. No floating point is used.

The symbolic checks establish `h*k = 1+t+...+t^6`, `h*(t-1)*k = t^7-1`, and the proposed seven-site matrix identity modulo the integral relation `lambda^2-lambda+2`. Both basis determinants are one, and their inverse entries belong to `Z[lambda]`. A separate exact Laurent-adjoint calculation checks the symplectic lifts of the unit-block shears.

For q = 2, 4, 8, 16, 32, 64 and both roots of the source parameter equation, the checker reproduces the source's cyclic-check annihilator cardinality q^3 at L = 7 and 14. It checks the complete 98-coordinate tensor complex after evaluation at two pairs of coarse-variable units, its reduction to nine toric summands and forty two-qudit product summands, and the physical symplectic transformation at coarse variables one. The source check is v = -h; changing the overall check sign leaves the stabilizer module unchanged. These results give eighteen encoded q-dimensional qudits at L = 7, agreeing with the source's reported q = 4 case.

The finite-ring checks also verify the cubic roots, seventh-order unit, invertible CRT evaluation matrix on the nine-dimensional charge module, and the fine translation fixed-charge directions (1,1), (2,1), (4,1). The q = 4 cubic and ring cardinality reproduce the source's h = t^3+3t^2+2t+3 and |GR(4,3)| = 64. The standard toric string crossing is the coordinate dot product modulo q. Finite specializations support the route and do not establish its uniform local-equivalence or physical braiding conclusion; those remain proof obligations.
