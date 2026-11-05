# Wang–Qin–Xia loopless stability: exact finite checks

Run `python3 docs/reports/wang-qin-xia-loopless-stability/check.py` from the repository root with Python 3.10 or later. The standard-library program uses Boolean adjacency matrices, integers, finite sets, permutations and exact bit-mask recursion. Its captured output is `output.txt`; a successful run exits 0 and ends with `ALL_CHECKS_OK`.

The source is Xiaomeng Wang, Yan-Li Qin and Binzhou Xia, “The existence of unexpected automorphisms in direct product graphs”, arXiv:2509.26170v2, Conjecture 1.3. The fidelity check reproduces the Section 2 equality S(C_n)=B(C_n), with no dispensable edges, for n=3,5,6,7,8,9,10. A second fidelity check reproduces the neighborhoods, three loops and diagonal two-fold automorphisms of Wang–Gao's four-vertex graph Δ from arXiv:2606.02203v1, Section 5.

For m=3,…,8 the proposed graph Γ_m has vertices (t,a,b) in {0,…,m−1} × {0,1}² and adjacency t≠u, R(a,d), R(c,b), where R={(0,1),(1,0),(1,1)}. The program checks its precise neighborhoods, looplessness, connectedness, absence of twins, a triangle, connected canonical double cover, and independence number m+2. It checks the proposed noncomponentwise involution on Γ_m × P₄ against all ordered pairs.

For Γ₃ it exhausts all degree-preserving row permutations. In a twin-free matrix, each column map is uniquely determined by the row map, so this finds all two-fold automorphisms, including any nondiagonal ones. There are 12, all diagonal. It also enumerates all loop-allowed factors of P₄: the only connected nonbipartite candidates are the two labelled copies of the edge with exactly one loop. Their unlooped slices force a half-size independent set in any direct product, conflicting with the checked independence numbers.

These results support a preregistered refutation route. They do not prove the assertions for every m≥3 and do not constitute a Lean certificate. No Lean proof is included.

The checker is independently written for these definitions; no source-paper code or external data is copied. The cited articles supply mathematical definitions and comparison results.
