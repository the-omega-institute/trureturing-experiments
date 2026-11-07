# Mixed distance-law classification: exact helper checks

Source: Christoph Thäle and Philipp Tuchel, *Realising atomless laws as distance distributions on metric measure spaces*, arXiv:2608.28330v1, Section 6(i). The proposed answer is that a Borel probability Θ on the nonnegative half-line is feasible exactly when zero belongs to its support and either Θ is concentrated on a countable set or Θ((0,ε)) is positive for every positive ε. Feasibility uses one Polish metric probability space and two independent samples. These checks do not prove that answer.

Run `python3 check.py`. The program uses rational arithmetic and exact SymPy expressions. The stored output ends in `ALL_CHECKS_OK`.

The fidelity check reproduces the distance-law formula in the source's lemma labelled `lem:distance-identity(ii)`: internal coefficients are p_i² and the coefficient for a later component is 2p_i times the total earlier sampling mass. It directly enumerates all point pairs in two-point components, for every positive integer weight composition with total 2 through 9 and between 2 and 5 components. The same examples check the square telescoping identity and every ordered triple's triangle inequality. The sufficient small-weight bound is tested only when its required lower bound on p_0 holds.

The local models are explicit. On X=[0,1], take uniform sampling and d(t,u)=|t−u|/4. Its distance law ν has density 2(L−r)/L² on [0,L], with L=1/4. This is derived from the area of the diagonal strip; eight independently integrated moments cross-check the coefficients. For any R, the map φ(t)=R+max(0,2t−1)/8 is 1-Lipschitz and has pushforward (δ_R+Uniform[R,R+1/8])/2. Its atom has mass 1/2 and the slope on the other half is 1/4.

For ζ=(δ_1+Uniform[0,1])/2, choosing β=1/32 gives βν density at most 1/4, below ζ's continuous density 1/2. For ζ=(δ_0/4+3Exp(1)/4) restricted to [0,1], the same β works: exp(−r)≥1−r gives continuous density at least 9/16 on [0,1/4]. The program verifies the exact analytic certificate for that inequality: the difference and its first derivative vanish at zero and its second derivative is exp(−r)>0. The model has no distance atom at zero, so subtraction preserves ζ's zero atom. Thus both requested mixed donors pass the local domination check, including a genuinely mixed image law.

An additional donor is purely atomic: the uniform four-ary product space, with distance 4^(−n) at the first differing coordinate n, has law Σ_(n≥1) 3/4^n δ_(4^(−n)). All pair coefficients are enumerated at depths 1 through 4 with the exact unresolved prefix mass 4^(−depth). A scaled base-four coordinate followed by the same plateau map gives a Lipschitz mixed image. The geometric normalization and the tail oscillation bound are checked symbolically.

Repeated normalized bins split outer atoms without splitting the atoms as measurable sets. The program checks coefficient conservation and the ordered-bin inequality used by the unbounded gluing argument.

The general nonatomic lift, measurable compact completion, countable assembly and universal iff remain paper-level arguments for subsequent verification. No Lean proof or refutation is included. The checks support helper identities; the planned positive settlement uses a general construction, with utility `none`.
