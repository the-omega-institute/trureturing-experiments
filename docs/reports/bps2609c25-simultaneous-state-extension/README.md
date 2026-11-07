# Simultaneous state extension: exact finite checks

Source: Dominika Burešová, Pavel Pták and Jan Ševic, *Simultaneous extensions of states on quantum logics*, arXiv:2609.10770v1, §2, conjecture following Proposition 2.7.

Run `python3 check.py` with Python 3. The script uses only the standard library, integers and `fractions.Fraction`; the complete deterministic output is in `check.out`.

The source-fidelity check reconstructs the 38 atoms and 21 maximal Boolean blocks of the diagram in the source's `simplified_js_2.tex`. It verifies the critical state with a = b = 1 and all ten named internal atoms assigned 1/2, and exhausts the 32 possible critical two-valued assignments. It also tests transportation for the two end blocks, including the forbidden m = n = 1 pair.

The transportation check compares an independent integer augmenting-path solver with every weighted Hall inequality for all bipartite graphs with one to three vertices in each part and all marginals with denominator 4, including zero masses: 125,482 cases, 32,651 feasible. Each resulting matrix is checked for nonnegativity, support and both marginals.

The cut and lifting check tests all pairs of Boolean subalgebras represented by partitions of sets with one to four points: 1,534 subset cuts and 8,271 rational common-source marginal cases. It verifies ambient subset inequalities and lifts couplings through Dirac states.

These are finite checks of the proposed proof mechanism. They do not prove the conjecture for arbitrary orthomodular posets, arbitrary real states, or infinite Boolean subalgebras. The planned settlement uses a general transportation argument and product compactness.
