# Arbitrary-spin matrix-product scar checks

`check.py` checks the periodic Eq. (23) construction in Sascha Gehrmann and Fabian H. L. Essler, *Exact Quantum Many-Body Scars by a generalized Matrix-Product Ansatz*, arXiv:2605.03020v2, Section IV. Run `python3 check.py` with Python 3 and SymPy. The stored stdout is `check.out`; successful execution exits zero and ends with `ALL_CHECKS_OK`.

All arithmetic uses exact integers, rationals, complex units and algebraic square roots. The physical basis index k denotes |S−k>, with n=2S. The full spin Hamiltonian uses the oriented periodic bonds j→j+1, including the seam, and all three components of 4 D·(S_j×S_{j+1}). Its local matrix is constructed independently of the proposed cancellation tensors.

The fidelity tests reproduce the source's zero-energy, translation-invariant, constrained nonzero state for n=1,2,3,4 at N=4. Higher-spin checks cover (n,N)=(5,2),(5,4),(6,4),(7,4),(8,4),(5,6). These use a=3/5+i/7 and D=(2√n,3√n,5). Additional n=5,N=4 checks use D=(0,√5,1),(√5,0,−2),(0,0,2), including the collapsed vacuum family. The physical action checks every output reachable from the supported input words using the complete local physical matrix; this includes leakage to level 2. All other output coefficients vanish by the explicit local matrix support.

The independent symbolic calculation constructs the three-level Hamiltonian directly from ladder matrices with entries r,q and S_z=diag(S,S−1,S−2). It verifies all 36 scalar entries of hAA=EA−AE modulo r²=2S, leaving q symbolic. Omitting E₂ fails the leakage equation. The source gap formula is compared with all trace words at N=2,4,6. At N=2,4,6,8 its palindromic coefficient ranks, ranks of the first N/2+1 Taylor coefficients in 1/a, and ranks of N+1 distinct admissible parameter samples are 2,3,4,5. The N=8 rank check uses the gap formula and is not an independent trace-formula comparison.

These checks support a preregistered universal algebraic argument. They do not prove the all-spin, all-length theorem in Lean, determine the full zero-energy eigenspace, or establish ETH violation or chaotic dynamics.
