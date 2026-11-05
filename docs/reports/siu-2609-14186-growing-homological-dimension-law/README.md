# Growing homological dimension for Siu tensor landscapes

Source: Chunyin Siu, *The Global Topology of Orthogonally Decomposable Tensor Landscapes*, arXiv:2609.14186v1. The target is the explicit growing-dimension question in Section 1.2, using Theorems 3.1, 3.3, Corollary 3.2 and Definition 3.7.

Run the exact checks with Python 3 and SymPy:

    python3 check.py

The program uses integers, fractions and exact SymPy algebra only. It checks the source's singleton critical-value identity and the q=0 and q=D−2 multiplicities, then the candidate's subset counts, exact finite-population conditional mean and variance, geometric tails, odd-family mixture masses, derivative coefficient, and half-dimension scaling. Every D from 2 through 10 and all permissible subset sizes are enumerated. Constant-array proportional checks also use D=200,400. The final output is ALL_CHECKS_OK. These checks verify finite algebra and local expansions; they do not prove the random-matrix limit or the persistence identification.

The separately labelled exploratory simulation requires NumPy:

    OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 simulate.py

It samples one real Gaussian Wishart matrix at each D=200,400, with n=2D; uses k=3,4 and q=floor(theta D) for theta=1/3,1/2; and samples 8000 bars per family from the exact multiplicity distribution. The MP integral uses 256-point Gauss–Legendre quadrature after a cosine substitution. Its floating-point results in simulation.jsonl are exploratory readings, never exact evidence or a proof of convergence. Odd superlevel readings apply (x,y) ↦ (−x−y/D,y) to the sublevel readings.

For theta=1/2 and k=4 the target lifetime is 0.0857864376. The sampled sublevel lifetime means are 0.0928002017 at D=200 and 0.0890102698 at D=400; superlevel means are 0.0894587457 and 0.0870206553. For k=3 the target is 0.00245310429; the negative-family means are 0.00274866953 and 0.00257056607, and positive-family means are 0.00283743797 and 0.00260677211. Full endpoint means and standard deviations, both theta values and both filtrations are in the output. One realization per size does not quantify sampling uncertainty across matrices.

The general proof route uses weighted-subset concentration, a uniform lifetime expansion, exact odd-family binomial masses, and pathwise Wishart spectral convergence. Formalization remains a research line: almost-sure Marchenko–Pastur convergence, both spectral edges (with the coupling and quantifiers needed for triangular arrays), and the bridge from these facts to weak empirical-measure limits are obligations. The deterministic barcode identification is supplied by the source, not by these experiments.
