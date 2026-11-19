# Aboulker–Aubian–Charbit–Lopes Question 3.5

Run `python3 check.py` with assertions enabled. The checker uses only the Python standard library and integer arithmetic. Its output is stored in `check.out`.

The fidelity check reproduces the example following Theorem 3.3 of arXiv:2310.04265v2 for transitive tournaments of orders 1 through 6. The separator check tests every two-color assignment and every reversed-endpoint order for both orientations of the endpoint arc. For the Mycielski constructions at q = 5 and q = 6, the program verifies the exact chromatic number of the base by exhaustive coloring search, the displayed backedge graph, and a dicoloring with q − 1 colors. The endpoint check and the forced base order then certify the clique-optimal chromatic minimum at each tested q.
