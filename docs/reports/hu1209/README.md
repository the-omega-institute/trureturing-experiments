# Exact checks for Hu–Fan–Zhou–Liu QCP candidate

The source is arXiv:1209.0058v2, *Super-additivity of quantum-correlating power*, and Physical Review A 88, 012315, DOI 10.1103/PhysRevA.88.012315. The proposed target is the contextual strict reading of its unnumbered conjecture at the end of Section III. The printed conjecture omits the word “strict”; the already-proved non-strict inequality is not a refutation target.

Run `python3 check.py` with Python 3.12 or newer and SymPy 1.14.0. The stored output is `check.out`. The program uses integer, rational and symbolic algebra, with no floating-point entropy estimates.

The fidelity check reproduces the source's claim that its two phase-damped orthogonal inputs have noncommuting outputs. For computational-basis damping with off-diagonal multiplier r and p=1-r², direct matrix multiplication gives the Pauli factor Y⊗I+I⊗Y, whereas the source displays Y⊗Z+Z⊗Y. Both give a nonzero operator for 0<r<1. The program checks the directly computed identity and the source's nonzero conclusion, not the mismatching printed Pauli factor.

For the proposed two-parameter five-dimensional family, the program checks POVM completeness, characteristic polynomials, exact lower eigenvalue endpoints, selector outputs, the output noise floor, noncommuting outputs and noncommuting adjoint range. At a=t=1/2 it independently assembles the Choi matrix in two ways, verifies positive definiteness, and checks its output partial trace equals I₅.

The final checks verify the symbolic entropy/information algebra used in the correlated-prior upper bound and independent-prior lower bound. They also check the conditional-independence probability identity for all 80 nonzero binary arrays with weights in {0,1,2}. These are structural certificates supporting the mathematical outline. They do not optimize quantum discord, prove data processing, search all POVMs, or establish the full QCP equality by computation. That equality requires the general analytic argument in the preregistration, and remains without a Lean proof at this stage.

The channel class itself is not new: essentially classical-quantum channels appear in Fukuda–Nechita–Wolf, arXiv:1408.2340, Definition 2.1(ii). Their image-additivity theorem is distinguished from the proposed optimized-discord conclusion.
