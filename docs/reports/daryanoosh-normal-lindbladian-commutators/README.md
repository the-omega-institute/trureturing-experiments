# Normal Lindbladian with three nonzero commutator sectors

Source: Shakib Daryanoosh, *From normal Lindbladians to non-normal quantum trajectories*, [arXiv:2608.04775v1](https://arxiv.org/abs/2608.04775v1), Section VIII and the second construction challenge in Section IX.

Run `python3 check.py` with Python 3 and SymPy. The recorded output is in `output.txt`. The program uses exact symbolic algebra, with no floating-point evidence.

The fidelity check reproduces the source's driven thermal qubit: its displayed Pauli-coordinate generator, normality, and three separately zero commutators, with gamma_th = q².

The proposed construction takes 0 < r < t and a,b,h > 0, jumps rX+aI, rY+ibI, tZ, and Hamiltonian hZ-rbY. It checks the GKSL reduction and all three nonzero commutator formulas symbolically. Their sum is zero. It also checks the matrix-unit representation independently through Kronecker products, the coordinate change, the characteristic polynomial, and the stationary density I/2.

These are mathematical checks of the proposed family and a source example. They do not establish Lean kernel verification, independent peer review, literature novelty, a new gauge-equivalence class, or a simulation advantage. The family uses the identity-shift freedom explicitly allowed by the source.
