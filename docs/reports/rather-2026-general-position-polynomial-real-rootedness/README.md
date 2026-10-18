# General position polynomials: which graphs are real-rooted (Rather, arXiv:2603.06930v3, Section 8 item 1)

`check.py N` enumerates every labelled simple graph on at most N vertices, computes its general position polynomial from the definition (a set is in general position if no vertex of it lies on a shortest path between two others), counts real roots exactly with sympy, and checks: connected graphs are real-rooted iff complete; all graphs are real-rooted iff every component is complete; alpha_0 = 1, alpha_1 = n, alpha_2 = n(n-1)/2; real-rooted polynomials have degree n.

Run: `python3 check.py 6` with sympy 1.14.0 (exit 0, `RESULT PASS`, 33867 graphs; output in `out.txt`).
