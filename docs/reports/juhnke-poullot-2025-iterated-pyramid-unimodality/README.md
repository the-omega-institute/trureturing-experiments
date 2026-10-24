# Juhnke–Poullot Problem 3.12 checks (arXiv:2504.20739v3)

Each script is run from this directory.
- `check_fidelity.py`: computes edges by exact rational LP (sympy); counts monotone paths exactly and coherent paths by sampled projections. Reproduces Example 3.5 (3-cube: six monotone paths of length 3, all coherent) and Example 4.1 (Lop₃, c = (1, 1, 1): two paths of length 2 and four of length 4, all coherent). Checks N(Pyr^k) = (1+X)^(k−1) N(Pyr¹) for k = 2, 3. Run with `/opt/homebrew/bin/python3`.
- `check_coherent_exact.py`: decides edges and coherence exactly by LP feasibility (HiGHS). A path is coherent iff some ω puts every other projected vertex strictly below each projected step. It covers three random polytopes whose coherent paths form a strict subset of the monotone paths. Run with a Python that has numpy and scipy.
- `check_edges_exact.py`: computes the edges of the simplicial random 4-polytope exactly from its facets (rational arithmetic) and cross-checks the HiGHS edge test.
- `check_unimodal_bound.py`: for random and adversarial nonnegative Q of degree d ≤ 7, checks that (1+X)^m Q is unimodal when m = 4d².
