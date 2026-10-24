# Strongly polynomial graph families with constant coefficient depth

Source: T. Bogart and S. Jiménez-Salazar, *Strongly regular graph, strongly polynomial sequences, and two-level polynomials*, arXiv:2609.32077v1, Question 3.19 and Proposition 3.7(1).

Run `python3 check.py` with Python 3. The program uses the standard library and exact integer arithmetic. `output.txt` contains its output.

The fidelity check reproduces the source identity hom(S_q,S_n)=n^q+n by enumerating every vertex map in 30 cases. The proposed affirmative construction is G_q=K_{q+3,2}, H_n=K_{1,n,n}. Nine full-map enumerations and 42 common-neighbor counts check its proposed homomorphism polynomial. Proper-coloring fibre counts are checked for all 76 labelled simple graphs on at most four vertices. Exact coefficient expansions test the first three codegrees and finite differences test the exponential residual.

These finite checks support a preregistration proof outline. They do not prove the unbounded assertion, establish literature novelty, or replace a Lean proof. All program code is independently authored for these checks; no source-paper code is copied.
