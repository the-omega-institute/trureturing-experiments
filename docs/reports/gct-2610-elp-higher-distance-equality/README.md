# Higher-distance nonlinear ELP equality checks

Source: Songping Ge, Han Cai and Xiaohu Tang, *Entropy/Length Profile Bounds for the Conversion of (Non)linear Codes in the Merge Regime*, arXiv:2610.15252v1, Section VI-C, Construction C, Proposition 3, Corollary 5 and the following open question.

Run `python3 check.py` with Python 3.12 or newer and NumPy. The stored output is `check-output.txt`; successful completion prints `ALL_CHECKS_OK`.

The fidelity check enumerates the source's quasigroup construction with m=3 and t=2. It verifies parameters (4,27,2) and (7,729,2) over alphabet 5, all 144 coordinate subsets, the stated profile value and inverse, and both conversion resource equalities. It uses the induced projection probabilities, rather than a uniform law on the projection image.

The proposed family uses a Cauchy parity matrix and systematic alphabet {0,...,q-2}. The exact checks cover (r,t,m,q)=(2,2,3,11) and (2,2,3,13). For each case they check every nonempty square minor, every coordinate subset in both initial codes and the final code, the exact minimum distance 3, both resource equalities, and strict nontriviality. All 1,000,000 and 2,985,984 conversion input tuples are enumerated.

For each coordinate subset, enumerate every assignment of its missing systematic symbols and group the observed parity values. Fixed observed systematic symbols only translate those parity values, so this checks every conditional fibre size without repeating identical translated histograms. The program verifies the upper bound (q-1)^max(0,j-r) for every fibre and constant-size fibres attaining it when all parity coordinates are selected. When j<=r all fibres are singletons. These exact size checks yield the profile h_j=max(0,j-r) log_q(q-1), including the induced conditional entropy weights. Singleton fibres for all selected sets of size r exclude all pairs at distance at most r; a unit systematic difference attains distance r+1.

Only int64 arrays, integers and rational coefficients of symbolic logarithms enter the evidence. All logarithmic comparisons needed for rounding and nontriviality are reduced to integer power comparisons. The arrays used here stay below int64 overflow. These finite checks support a preregistration; they are not a proof of the unbounded family and carry no Lean verification claim.
