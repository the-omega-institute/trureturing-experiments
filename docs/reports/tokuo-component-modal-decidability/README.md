# Component-modal decidability: exact finite checks

`check.py` uses Python 3 standard-library integer bit masks and exhaustive deterministic enumeration. Run `python3 check.py`; the expected output is stored in `output.txt`, ending in `ALL_CHECKS_OK`.

The reference is Kenji Tokuo, *Component Modalities of Quantum Logic*, arXiv:2607.14757v1, Sections 2 and 4, https://arxiv.org/html/2607.14757v1. The program is original experimental code; it does not embed the article or third-party implementation.

Fidelity checks reproduce the two-component modal-alternatives example from Section 2 and Theorem 4.2's component interior and cover identities on every reflexive symmetric graph with 1–5 labelled vertices and every stable subset. There are 1,099 graphs and 9,074 stable subsets.

The filtration check enumerates every connected reflexive symmetric graph with 1–4 labelled vertices, every pair of stable atomic valuations, and all 122 formulas of depth at most 2 over two named atoms. For each formula it adjoins the negations of occurring atoms, identifies equal truth types, and takes the existential image of compatibility. It checks reflexivity, symmetry, connectedness, stability of the image valuations, the bound of 2 to the power of the subformula count plus atom count, and preservation of every subformula's truth. Atoms absent from the input are assigned the empty set. This covers 1,148 models and 140,056 quotients; 120,224 quotients strictly reduce the number of worlds.

These finite checks support the proposed general filtration argument. They do not prove the unbounded finite-model property, source completeness, correctness of the Lean decision procedure, or decidability of the calculus. Those remain proof obligations after preregistration.
